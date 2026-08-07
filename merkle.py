# Author: Brian Soltani
# Date: 08/07/2026
# Description:
#   Merkle tree utilities which build a tree from a list of leaf values, get the root as a
#   commitment, and verify branch proofs that a leaf was part of the committed data.

from hashlib import sha256

from utils import log_2


def hash(x: bytes) -> bytes:
    """SHA-256 hash of x."""
    return sha256(x).digest()


def build_merkle(values: list[bytes]) -> list[bytes | None]:
    """
    Builds a Merkle tree over values as a 1-indexed array: index 1 holds the root, and the
    last len(values) slots (in original order) hold the hashed leaves. Leaves are hashed with
    a 0x00 prefix and internal nodes with a 0x01 prefix, so an internal node can never be
    replayed as a leaf.
    """
    leaf_count = len(values)
    assert leaf_count > 0 and (leaf_count & (leaf_count - 1)) == 0
    tree = [None] * leaf_count + [hash(b'\x00' + v) for v in values]
    for i in range(leaf_count - 1, 0, -1):
        tree[i] = hash(b'\x01' + tree[2 * i] + tree[2 * i + 1])
    return tree


def get_root(tree: list[bytes | None]) -> bytes:
    """Returns the Merkle root: the hash stored at index 1."""
    return tree[1]


def get_branch(tree: list[bytes | None], position: int) -> list[bytes]:
    """Collects the sibling hashes needed to prove that the leaf at position belongs to tree, walking up from the leaf and taking one sibling per level."""
    leaf_count = len(tree) // 2
    node_index = position + leaf_count
    level_count = log_2(leaf_count)
    branch = []
    for level in range(level_count):
        current_index = node_index >> level
        sibling = current_index ^ 1
        branch.append(tree[sibling])
    return branch


def verify_branch(root: bytes, position: int, value: bytes, branch: list[bytes]) -> bool:
    """
    Recomputes the hash path from a claimed leaf value and its branch up to the root, and
    checks it against root. Only returns True if value was part of the data committed by root.
    """
    current = hash(b'\x00' + value)
    for sibling in branch:
        if position % 2 != 0:  # right child
            current = hash(b'\x01' + sibling + current)
        else:  # left child
            current = hash(b'\x01' + current + sibling)
        position = position // 2
    return current == root
