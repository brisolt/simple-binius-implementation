# Author: Brian Soltani
# Date: 08/06/2026
# Description:
#   The simple Binius prover and verifier. Ties together the field (binary_fields.py), the
#   polynomial and Reed-Solomon utilities (utils.py), and the Merkle commitment (merkle.py)
#   into a full commit-prove-verify protocol.

from typing import Any

from utils import *
from merkle import *

EXPANSION_FACTOR = 8
NUM_CHALLENGES = 32
BYTES_PER_ELEMENT = 2


def choose_row_length_and_count(log_evaluation_count: int) -> list[int]:
    """
    Splits a log2 evaluation count into row/column dimensions for the commitment grid, factoring
    it as close to a square as possible: row length is rounded down, row count rounded up.
    """
    log_row_length = log_evaluation_count // 2
    log_row_count = (log_evaluation_count + 1) // 2
    row_length = 2 ** log_row_length
    row_count = 2 ** log_row_count
    return [log_row_length, log_row_count, row_length, row_count]


def derive_challenges(transcript: bytes, extended_row_length: int) -> list[int]:
    """
    Derives the Fiat-Shamir challenge column indices from transcript: repeatedly hashes
    transcript with an incrementing counter and takes the result mod extended_row_length,
    skipping duplicates, until NUM_CHALLENGES distinct indices are collected (or every column
    has been picked).
    """
    challenges = []
    counter = 0
    while len(challenges) < min(NUM_CHALLENGES, extended_row_length):
        candidate = int.from_bytes(hash(transcript + counter.to_bytes(4, 'little')), 'little') % extended_row_length
        if candidate not in challenges:
            challenges.append(candidate)
        counter += 1
    return challenges


def simple_binius_proof(evaluations: list[BinaryFieldElement], point: list[BinaryFieldElement]) -> dict[str, Any]:
    """
    Commits to evaluations (a multilinear polynomial's values on the boolean hypercube) and
    proves its evaluation at point, using a row/column commitment scheme.
    """
    log_evaluation_count = log_2(len(evaluations))
    if log_evaluation_count is None:
        raise ValueError("Must be a positive integer of power of two")
    if len(point) != log_evaluation_count:
        raise ValueError(f"point has {len(point)} coordinates. {log_evaluation_count} expected.")
    log_row_length, log_row_count, row_length, row_count = choose_row_length_and_count(log_evaluation_count)

    # 1. Rearrange evaluations into a row_count x row_length matrix.
    rows = []
    for i in range(0, len(evaluations), row_length):
        rows.append(evaluations[i : i + row_length])

    # 2. Extend each row using the Reed-Solomon code from utils.py.
    extended_rows = [extend(row, EXPANSION_FACTOR) for row in rows]
    extended_row_length = row_length * EXPANSION_FACTOR

    # 3. Fold the rows together into t_prime, using the row half of the evaluation point.
    row_combo = evaluation_tensor_product(point[log_row_length:])
    t_prime = [sum([rows[i][j] * row_combo[i] for i in range(row_count)], BinaryFieldElement(0)) for j in range(row_length)]

    # 4. Pack the extended columns and commit to them with a Merkle tree.
    columns = [[row[j] for row in extended_rows] for j in range(extended_row_length)]
    packed_columns = [b''.join(element.to_bytes(BYTES_PER_ELEMENT, 'little') for element in column) for column in columns]
    tree = build_merkle(packed_columns)
    root = get_root(tree)

    # 5. Bind everything into the Fiat-Shamir transcript and derive the challenge columns.
    claimed_eval = multilinear_polynomial_evaluation(evaluations, point)
    transcript = root
    for coord in point:
        transcript += coord.to_bytes(BYTES_PER_ELEMENT, 'little')
    transcript += claimed_eval.to_bytes(BYTES_PER_ELEMENT, 'little')
    for t in t_prime:
        transcript += t.to_bytes(BYTES_PER_ELEMENT, 'little')

    challenges = derive_challenges(transcript, extended_row_length)

    return {
        "root": root,
        "evaluation_point": point,
        "eval": claimed_eval,
        "t_prime": t_prime,
        "columns": [columns[c] for c in challenges],
        "branches": [get_branch(tree, c) for c in challenges]
    }


def verify(proof: dict[str, Any]) -> bool:
    """
    Re-derives the challenge columns from the proof's transcript, checks each opened column's
    Merkle branch against the committed root, re-extends t_prime and cross-checks it against
    the opened columns, then folds t_prime down to a single value and checks it against the
    claimed evaluation.
    """
    # 1. Recompute the grid dimensions from the evaluation point. The verifier never sees the
    #    original evaluations, only the point, and len(point) is already the log2 count, so no
    #    extra call is needed.
    log_row_length, log_row_count, row_length, row_count = choose_row_length_and_count(len(proof["evaluation_point"]))
    extended_row_length = row_length * EXPANSION_FACTOR

    # 2. Rebuild the transcript exactly as the prover did, and re-derive the challenge columns.
    transcript = proof["root"]
    for coord in proof["evaluation_point"]:
        transcript += coord.to_bytes(BYTES_PER_ELEMENT, 'little')
    transcript += proof["eval"].to_bytes(BYTES_PER_ELEMENT, 'little')
    for t in proof["t_prime"]:
        transcript += t.to_bytes(BYTES_PER_ELEMENT, 'little')

    challenges = derive_challenges(transcript, extended_row_length)

    if len(proof["columns"]) != len(challenges) or len(proof["branches"]) != len(challenges):
        return False

    expected_branch_length = log_2(extended_row_length)
    for branch in proof["branches"]:
        if len(branch) != expected_branch_length:
            return False

    # 3. Verify each opened column's Merkle branch against the committed root. Pairs up the
    #    three lists by position: the i-th challenge to the i-th branch and column.
    for challenge, branch, column in zip(challenges, proof["branches"], proof["columns"]):
        packed_column = b''.join(x.to_bytes(BYTES_PER_ELEMENT, 'little') for x in column)

        if not verify_branch(proof["root"], challenge, packed_column, branch):
            return False

    # 4. Re-extend t_prime and cross-check it against the opened columns: the prover folded the
    #    original rows into t_prime, so folding the opened column values the same way should
    #    reproduce the same extension. This relies on the code's linearity, and on row_combo
    #    being computed with the same formula the prover used, from an evaluation point the
    #    verifier also has.
    extended_t_prime = extend(proof["t_prime"], EXPANSION_FACTOR)
    row_combo = evaluation_tensor_product(proof["evaluation_point"][log_row_length:])
    for column, challenge in zip(proof["columns"], challenges):
        expected = sum([column[i] * row_combo[i] for i in range(row_count)], BinaryFieldElement(0))

        if expected != extended_t_prime[challenge]:
            return False

    # 5. Fold t_prime down to a single value using the other half of the evaluation point. Since
    #    t_prime already collapsed the row dimension during proving, this collapses the column
    #    dimension, leaving a single value that should match the prover's claimed eval.
    col_combo = evaluation_tensor_product(proof["evaluation_point"][:log_row_length])
    computed_evaluation = sum([proof["t_prime"][i] * col_combo[i] for i in range(row_length)], BinaryFieldElement(0))
    if computed_evaluation != proof["eval"]:
        return False
    return True
