# Author: Brian Soltani
# Date: 08/06/2026
# Description:
#   Markle Tree utilities which builds a tree from a list of leaf values, gets root as a commitment, and verifies branch proofs
#   that a leaf was part of the committed data.

from hashlib import sha256
from utils import log_2

# Implements a SHA-256 hashing.
def hash(x):
    return sha256(x).digest()

# First n slots are internal nodes filledin and last n slots are hashed leaves in original order.
def build_merkle(values):
    n = len(values)
    assert n > 0 and (n & (n - 1)) == 0
    # tree = [None] * n + [hash(v) for v in values]
    tree = [None] * n + [hash(b'\x00' + v) for v in values]
    for i in range(n - 1, 0, -1):
        # tree[i] = hash(tree[2 * i] + tree[2 * i + 1])
        tree[i] = hash(b'\x01' + tree[2 * i] + tree[2 * i + 1])
    return tree

# Pretty trivial
def get_root(tree):
    return tree[1]

# Collects sibling hashes neede to prove that a specific leaf belongs to a tree, walking up and grabbing one sibling per level.
def get_branch(tree, position):
    n = len(tree) // 2 # Number of leaves
    pos = position + n # Convert leaf index to actual array
    level_count = log_2(n)
    branch = []
    for level in range(0, level_count):
        curr = pos >> level
        sibling = curr ^ 1
        branch.append(tree[sibling])
    return branch

# Given a claimed leaf value and its branch, recomputes the hash patch up to the root and checks if it matches.
#   -> ONLY RETURNS TRUE IF VALUE WAS PART OF DATA COMMITTED BY ROOT!!!
def verify_branch(root, position, value, branch):
    current = hash(b'\x00' + value)
    for sibling in branch:
        if position % 2 != 0: # We are the right child
            current = hash(b'\x01' + sibling + current)
        else: # Meaning position is even so left child
            current = hash(b'\x01' + current + sibling)
        position = position // 2
    return current == root