# Author: Brian Soltani
# Date: 08/06/2026
# Description:
#   binary_fields class for my implementation of a simple Binius.

from __future__ import annotations

from typing import Optional


def binary_multiplication(v1: int, v2: int, length: Optional[int] = None) -> int:
    """
    Multiplies two binary tower field elements, represented as integers, using the tower's
    Karatsuba-style folding rule: split each value into a low half L and a high half R
    (value = L + R*X), then combine as L1*L2 + (L1*R2 + R1*L2)*X + R1*R2*X^2, where X^2
    reduces via the tower relation X_k^2 = X_k*X_(k-1) + 1.

    Base case is plain integer multiplication for elements of GF(2), i.e. 0 or 1.
    """
    if v1 < 2 or v2 < 2:
        return v1 * v2
    if length is None:
        length = 1 << (max(v1, v2).bit_length() - 1).bit_length()
    half = length // 2
    cut = (1 << half) - 1
    l1 = v1 & cut
    r1 = v1 >> half
    l2 = v2 & cut
    r2 = v2 >> half

    quarter = length // 4
    x_kminus1 = 1 << quarter

    l1l2 = binary_multiplication(l1, l2, half)
    r1r2_reg = binary_multiplication(r1, r2, half)
    r1r2_reduced = binary_multiplication(x_kminus1, r1r2_reg, half)

    cross_combined = binary_multiplication(l1 ^ r1, l2 ^ r2, half)
    cross_term = cross_combined ^ l1l2 ^ r1r2_reg

    low_half = l1l2 ^ r1r2_reg
    high_half = cross_term ^ r1r2_reduced

    return low_half ^ (high_half << half)


class BinaryFieldElement:
    """
    An element of a binary tower field, stored as a non-negative int. Addition is XOR (so
    every element is its own additive inverse); multiplication follows the tower rule
    implemented in binary_multiplication.
    """

    def __init__(self, value: int) -> None:
        """Wraps a raw integer as a field element."""
        self.value = value

    def __add__(self, other: BinaryFieldElement | int) -> BinaryFieldElement:
        """Adds via XOR. A plain int on the right is treated as its raw value."""
        if isinstance(other, BinaryFieldElement):
            other_value = other.value
        else:
            other_value = other
        return BinaryFieldElement(self.value ^ other_value)

    __sub__ = __add__  # subtraction is the same as addition in characteristic two
    __radd__ = __add__
    __rsub__ = __add__

    def __repr__(self) -> str:
        """Renders as the underlying integer value."""
        return f"{self.value}"

    def __mul__(self, other: BinaryFieldElement | int) -> BinaryFieldElement:
        """Multiplies via the tower field multiplication rule."""
        if isinstance(other, BinaryFieldElement):
            other_value = other.value
        else:
            other_value = other
        return BinaryFieldElement(binary_multiplication(self.value, other_value))

    def __pow__(self, other: int) -> BinaryFieldElement:
        """Exponentiation by repeated squaring, recursing down to powers of 0, 1, and 2."""
        if other < 0:
            raise ValueError("Negative exponents are unsupported. Use multiplicative_inverse()")

        if other == 0:
            return BinaryFieldElement(1)
        elif other == 1:
            return self
        elif other == 2:
            return self * self
        else:
            return self.__pow__(other % 2) * self.__pow__(other // 2) ** 2
            # pow(3, 5) = pow(3, 1) * pow(3, 2) = 3 * 3 = 2 in the binary field.

    def __eq__(self, other: object) -> bool:
        """Compares by underlying value, against another BinaryFieldElement or a raw int."""
        if isinstance(other, BinaryFieldElement):
            other_value = other.value
        else:
            other_value = other
        return self.value == other_value

    def __hash__(self) -> int:
        """Hashes by underlying value, so a BinaryFieldElement and the equal raw int hash the same."""
        return hash(self.value)

    def __truediv__(self, other: BinaryFieldElement | int) -> BinaryFieldElement:
        """Divides by multiplying with the multiplicative inverse of other."""
        if isinstance(other, BinaryFieldElement):
            divisor = other
        else:
            divisor = BinaryFieldElement(other)
        return self * divisor.multiplicative_inverse()

    def multiplicative_inverse(self) -> BinaryFieldElement:
        """
        Computes the multiplicative inverse via Fermat's little theorem, x^(N-2), where N is
        the size of the smallest sub-field of the tower that self.value lives in. The bit
        length of self.value (rounded up to the nearest tower tier) gives that sub-field's
        level, and N = 2^(2^level).
        """
        if self.value == 0:
            raise ZeroDivisionError("No multiplicative inverse")

        value_bit_length = self.value.bit_length()
        tower_level = (value_bit_length - 1).bit_length()
        subfield_bit_width = 2 ** tower_level
        field_size = 2 ** subfield_bit_width
        return self ** (field_size - 2)

    def to_bytes(self, length: int, byte_order: str) -> bytes:
        """Serializes the element's value to raw bytes."""
        return self.value.to_bytes(length, byte_order)

    @classmethod
    def from_bytes(cls, b: bytes, byte_order: str) -> BinaryFieldElement:
        """Deserializes raw bytes back into a BinaryFieldElement."""
        return cls(int.from_bytes(b, byte_order))
