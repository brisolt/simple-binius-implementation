# Author: Brian Soltani
# Date: 08/06/2026
# Description: its way too late to make this rn just trust me on this.

from utils import *
from merkle import *

EXPANSION_FACTOR = 8
NUM_CHALLENGES = 32
BYTES_PER_ELEMENT = 2

# Straightforward
def choose_row_length_and_count(log_evaluation_count):
    log_row_length = log_evaluation_count // 2
    log_row_count = (log_evaluation_count + 1) // 2
    row_length = 2 ** log_row_length
    row_count = 2 ** log_row_count
    return [log_row_length, log_row_count, row_length, row_count]

def derive_challenges(transcript, extended_rl):
    challenges = []
    counter = 0
    while len(challenges) < min(NUM_CHALLENGES, extended_rl):
        c = int.from_bytes(hash(transcript + counter.to_bytes(4, 'little')), 'little') % extended_rl
        if c not in challenges:
            challenges.append(c)
        counter += 1
    return challenges


def simple_binius_proof(evaluations, point):
    # lrl, lrc, rl, rc = choose_row_length_and_count(log_2(len(evaluations)))

    log_evaluation_count = log_2(len(evaluations))
    if log_evaluation_count is None:
        raise ValueError("Must be a positive integer of power of two")
    if len(point) != log_evaluation_count:
        raise ValueError(f"point has {len(point)} coordinates. {log_evaluation_count} expected.")
    lrl, lrc, rl, rc = choose_row_length_and_count(log_evaluation_count)

    rows = []
    for i in range(0, len(evaluations), rl): # Rearranges evaluations into a rl x rc matrix.
        rows.append(evaluations[i : i + rl])
    extended_rows = [extend(row, EXPANSION_FACTOR) for row in rows] # Extend rows using Reed-Soloman code from utilities
    extended_rl = rl * EXPANSION_FACTOR
    row_combo = evaluation_tensor_product(point[lrl:]) # Linear combo
    t_prime = [sum([rows[i][j] * row_combo[i] for i in range(rc)], BinaryFieldElement(0)) for j in range(rl)] # Fold rows to t_prime
    columns = [[row[j] for row in extended_rows] for j in range(extended_rl)]
    packed_columns = [b''.join(element.to_bytes(BYTES_PER_ELEMENT, 'little') for element in column) for column in columns] #
    tree = build_merkle(packed_columns) # Pack extended columns and comit to Merkle tree
    root = get_root(tree)

    # T/S: Fixing Fiat-Shamir binding
    claimed_eval = multilinear_polynomial_evaluation(evaluations, point)
    transcript = root
    for coord in point:
        transcript += coord.to_bytes(BYTES_PER_ELEMENT, 'little')
    transcript += claimed_eval.to_bytes(BYTES_PER_ELEMENT, 'little')
    for t in t_prime:
        transcript += t.to_bytes(BYTES_PER_ELEMENT, 'little')

    challenges = derive_challenges(transcript, extended_rl)

    return  { # finally.
        "root": root,
        "evaluation_point": point,
        "eval": claimed_eval,
        "t_prime": t_prime,
        "columns": [columns[c] for c in challenges],
        "branches": [get_branch(tree, c) for c in challenges]
    }
    

def verify(proof):
    lrl, lrc, rl, rc = choose_row_length_and_count(len(proof["evaluation_point"])) # grid dimensions. V doenst have evaluations, only the point. also len(evalpt) is already log2 count so no call needed
    extended_rl = rl * EXPANSION_FACTOR

    # T/S: more fiat fixes
    transcript = proof["root"]
    for coord in proof["evaluation_point"]:
        transcript += coord.to_bytes(BYTES_PER_ELEMENT, 'little')
    transcript += proof["eval"].to_bytes(BYTES_PER_ELEMENT, 'little')
    for t in proof["t_prime"]:
        transcript += t.to_bytes(BYTES_PER_ELEMENT, 'little')

    challenges = derive_challenges(transcript, extended_rl) # challenges (identical formula to what the prover computed)


    if len(proof["columns"]) != len(challenges) or len(proof["branches"]) != len(challenges):
        return False

    expected_branch_length = log_2(extended_rl)
    for branch in proof["branches"]:
        if len(branch) != expected_branch_length:
            return False
    

    for challenge, branch, column in zip(challenges, proof["branches"], proof["columns"]): # verify each opened column's merkle branch + pairs up three lists by position, ith challenge to ith branch and column. (matching order)
        packed_column = b''.join(x.to_bytes(BYTES_PER_ELEMENT, 'little') for x in column) 

        if not verify_branch(proof["root"], challenge, packed_column, branch):
            return False

    extended_t_prime = extend(proof["t_prime"], EXPANSION_FACTOR)
    row_combo = evaluation_tensor_product(proof["evaluation_point"][lrl:])
    for column, challenge in zip(proof["columns"], challenges): #re-extend t-prime and cross check against the opened columns
        expected = sum([column[i] * row_combo[i] for i in range(rc)], BinaryFieldElement(0))

        if expected != extended_t_prime[challenge]:
            return False

        # prover computed t' by folding OG rows together, now folding opened column valutes together same way. linearity and R-S being so is also important
        # + row combo is computed with same formula P used, depending on eval_pt, which V has.
    col_combo = evaluation_tensor_product(proof["evaluation_point"][:lrl])
    computed_evaluation = sum([proof["t_prime"][i] * col_combo[i] for i in range(rl)], BinaryFieldElement(0)) # fold t_prime down to final value and...
    if computed_evaluation != proof["eval"]:
        return False
    return True # check if it matches
    # folded across cols using back half of eval_pt. since t' already collapsed row dim during proving, this collapses the column dim, leaving
    # V with a single value, which should look like the claim "eval" that P claimed.