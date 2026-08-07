# Author: Brian Soltani
# Description:
#   A short guided tour of this simple Binius implementation

import copy
import random
import sys
import time

from binary_fields import BinaryFieldElement as B, binary_multiplication
from merkle import build_merkle, get_root, get_branch, verify_branch
from utils import extend, multilinear_polynomial_evaluation
from simple_binius import (simple_binius_proof, verify, choose_row_length_and_count,
                           EXPANSION_FACTOR, NUM_CHALLENGES, BYTES_PER_ELEMENT)

random.seed(2024)

WIDTH = 70
RESULTS = []

BOLD, DIM, GREEN, RED, RESET = "\033[1m", "\033[2m", "\033[32m", "\033[31m", "\033[0m"


def title(text):
    print(f"\n{BOLD}{text}{RESET}")
    print("=" * WIDTH)


def section(number, text):
    print(f"\n{BOLD}{number}. {text}{RESET}")
    print("-" * WIDTH)


def show(label, value):
    print(f"   {label:<34} {value}")


def note(text):
    print(f"   {DIM}{text}{RESET}")


def check(label, condition):
    RESULTS.append(condition)
    mark = f"{GREEN}PASS{RESET}" if condition else f"{RED}FAIL{RESET}"
    print(f"   [{mark}]  {label}")


# ------------------------------------------------------------------ 1. field
def demo_field():
    section(1, "BINARY TOWER FIELD")
    note("Elements are integers. Addition is XOR, so every element is its own")
    note("negative. Multiplication recurses down the tower with Karatsuba.")
    print()

    a, b = B(13), B(7)
    show("a, b", f"{a}, {b}")
    show("a + b   (XOR)", a + b)
    show("a - b   (same as a + b)", a - b)
    show("a * b   (tower multiply)", a * b)
    show("a + a   (char. 2)", a + a)
    print()

    show("a^-1", a.multiplicative_inverse())
    show("a * a^-1", a * a.multiplicative_inverse())
    check("every nonzero element of GF(2^4) inverts",
          all(B(x) * B(x).multiplicative_inverse() == 1 for x in range(1, 16)))
    print()

    note("The tower is built so each level satisfies  X_k^2 = X_k * X_{k-1} + 1,")
    note("and each level is closed: GF(2^4) products stay inside GF(2^4).")
    X, prev = B(1 << 2), B(1 << 1)
    show("X^2", X * X)
    show("X * X_prev + 1", X * prev + B(1))
    check("defining relation holds", X * X == X * prev + B(1))
    check("GF(2^4) is closed under multiply",
          all(binary_multiplication(x, y) < 16 for x in range(16) for y in range(16)))


# ----------------------------------------------------------------- 2. merkle
def demo_merkle():
    section(2, "MERKLE COMMITMENT")
    note("Leaves are hashed with a 0x00 prefix and internal nodes with 0x01, so")
    note("no internal node can ever be replayed as a leaf.")
    print()

    leaves = [bytes([i]) * 8 for i in range(8)]
    tree = build_merkle(leaves)
    root = get_root(tree)
    branch = get_branch(tree, 3)

    show("leaves committed", len(leaves))
    show("root", root.hex()[:32] + "...")
    show("branch length for leaf 3", f"{len(branch)}  (= log2 of {len(leaves)})")
    print()

    check("leaf 3 opens against the root", verify_branch(root, 3, leaves[3], branch))
    check("a forged leaf value is rejected", not verify_branch(root, 3, b"forgery!", branch))
    check("the right leaf at the wrong slot is rejected",
          not verify_branch(root, 4, leaves[3], branch))
    check("flipping one bit changes the root",
          get_root(build_merkle(leaves[:7] + [b"\xff" * 8])) != root)


# ------------------------------------------------------- 3. reed-solomon code
def demo_extend():
    section(3, "REED-SOLOMON EXTENSION")
    note("Rows are stretched to 8x their length by interpolating a polynomial")
    note("through them and evaluating it further out. Redundancy is what makes")
    note("spot-checking a few columns meaningful.")
    print()

    values = [B(3), B(9), B(12), B(5)]
    extended = extend(values, EXPANSION_FACTOR)

    show("original symbols", values)
    show("extended length", f"{len(values)} -> {len(extended)}  ({EXPANSION_FACTOR}x)")
    show("first 8 extended", extended[:8])
    print()

    check("original symbols are preserved", extended[:len(values)] == values)

    # Linearity is the property the whole protocol leans on.
    other = [B(1), B(4), B(15), B(2)]
    scalar = B(6)
    combined = extend([values[i] + scalar * other[i] for i in range(4)], EXPANSION_FACTOR)
    separate = [extend(values, EXPANSION_FACTOR)[i] + scalar * extend(other, EXPANSION_FACTOR)[i]
                for i in range(len(extended))]
    check("the code is linear:  ext(a + c*b) == ext(a) + c*ext(b)", combined == separate)


# --------------------------------------------------------------- 4. protocol
def demo_protocol():
    section(4, "COMMIT, PROVE, VERIFY")
    note("The prover commits to a multilinear polynomial, then proves what it")
    note("evaluates to at one point without revealing the whole polynomial.")
    print()

    num_vars = 8
    size = 2 ** num_vars
    evaluations = [B(random.randrange(65536)) for _ in range(size)]
    point = [B(random.randrange(65536)) for _ in range(num_vars)]
    lrl, lrc, rl, rc = choose_row_length_and_count(num_vars)

    show("polynomial", f"multilinear, {num_vars} variables")
    show("evaluations on the hypercube", size)
    show("arranged as a grid", f"{rc} rows x {rl} columns")
    show("after Reed-Solomon expansion", f"{rc} rows x {rl * EXPANSION_FACTOR} columns")
    print()

    start = time.time()
    proof = simple_binius_proof(evaluations, point)
    prove_time = time.time() - start

    start = time.time()
    accepted = verify(proof)
    verify_time = time.time() - start

    branch_bytes = sum(len(b) * 32 for b in proof["branches"])
    proof_bytes = (32
                   + len(proof["t_prime"]) * BYTES_PER_ELEMENT
                   + sum(len(c) for c in proof["columns"]) * BYTES_PER_ELEMENT
                   + branch_bytes)

    show("claimed evaluation", proof["eval"])
    show("columns opened", f"{len(proof['columns'])} of {rl * EXPANSION_FACTOR}")
    show("proof size", f"{proof_bytes:,} bytes")
    show("prove / verify time", f"{prove_time:.2f}s / {verify_time:.3f}s")
    print()

    check("the honest proof verifies", accepted is True)
    check("the claimed value is the true evaluation",
          proof["eval"] == multilinear_polynomial_evaluation(evaluations, point))
    return proof


# -------------------------------------------------------------- 5. soundness
def demo_soundness(proof):
    section(5, "TAMPERING IS CAUGHT")
    note("Every check in verify() returns False rather than asserting, so the")
    note("verifier still rejects when Python runs with -O and strips asserts.")
    print()

    def rejects(label, mutate):
        forged = copy.deepcopy(proof)
        mutate(forged)
        try:
            accepted = verify(forged) is True
        except Exception:
            accepted = False
        check(label, not accepted)

    rejects("a lie about the evaluation", lambda p: p.update(eval=p["eval"] + B(1)))
    rejects("a doctored commitment root", lambda p: p.update(root=bytes(32)))
    rejects("a single altered field element",
            lambda p: p["columns"][0].__setitem__(0, p["columns"][0][0] + B(1)))
    rejects("a tampered Merkle branch", lambda p: p["branches"][0].__setitem__(0, bytes(32)))
    rejects("withholding the opened columns",
            lambda p: (p.update(columns=[]), p.update(branches=[])))
    rejects("reordering the openings",
            lambda p: (p.update(columns=p["columns"][::-1]),
                       p.update(branches=p["branches"][::-1])))


def main():
    title("SIMPLE BINIUS  -  a polynomial commitment scheme over binary towers")
    print(f"{DIM}   Brian Soltani{RESET}")

    demo_field()
    demo_merkle()
    demo_extend()
    proof = demo_protocol()
    demo_soundness(proof)

    passed, total = sum(RESULTS), len(RESULTS)
    print("\n" + "=" * WIDTH)
    if passed == total:
        print(f"{BOLD}{GREEN}   {passed}/{total} checks passed{RESET}")
    else:
        print(f"{BOLD}{RED}   {passed}/{total} checks passed{RESET}")
    print("=" * WIDTH)
    print(f"{DIM}   Full stress-test suite: python3 stress_tests.py{RESET}\n")
    return 0 if passed == total else 1


sys.exit(main())
