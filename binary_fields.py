# Author: Brian Soltani
# Date: 08/06/2026
# Description:
#   binary_fields class for my implementation of a simple Binius.

# Base case: regular multiplication under GF(2). Checks if v1 or v2, which are the values of the two inputs, are 0 or 1.
# Uses folding and Karatsuba trick to make multiplicaion more efficient.
# Multiplying (L1 + R1 * X)(L2 + R2 * X) expands to L1L2 + (L1R2 + R1L2)X + R1R2X^2.
#     -> Tower rule: X^2 reduces to X * X_(k-1) + 1
#     -> L1L2 is low by low, R1R2_reg is high by high, R1R2_reduced applies tower reduction to R1R2 since it is high
# 
def binary_multiplication(v1, v2, length=None):
    if (v1 < 2 or v2 < 2):
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
    # return low_half + high_half[:half_length - len(low_half)]


class BinaryFieldElement:

    # Constructor class.
    def __init__(self, value):
        self.value = value

    # Adding class forcing exclusive OR (XOR) with the ^ bitwise operator. Adding two BinaryFieldElements will apply this XOR
    # but two non-BFEs will have regular addition in this class.
    def __add__(self, other):

        if isinstance(other, BinaryFieldElement):
            ov = other.value
        else:
            ov = other
        return BinaryFieldElement(self.value ^ ov)
    __sub__ = __add__ # Makes subtraction and addition the same structure.
    __radd__ = __add__
    __rsub__ = __add__

    # Printing reformating to make it clear to the reader.
    def __repr__(self):
        return f"{self.value}"

    # Wraps as BinaryFieldElement and calls a recursive multiplication class as a wrapper.
    def __mul__(self, other): 
        if isinstance(other, BinaryFieldElement):
            ov = other.value
        else:
            ov = other
        return BinaryFieldElement(binary_multiplication(self.value, ov))
        
    # Turns regular "**" to the binary-ok version where it recursively brings it down to powers of 0, 1, and 2.
    def __pow__(self, other):

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
            # pow(3, 5) = pow(3, 1) * pow(3, 2) = 3 * 3 = 2 in Binary Field.

    # Checks if two BinaryFieldElements are equal. Preeeetty straightforward.
    def __eq__(self, other):
        if isinstance(other, BinaryFieldElement):
            value = other.value
        else:
            value = other
        return self.value == value

    def __hash__(self):
        return hash(self.value)

    def __truediv__(self, other):
        if isinstance(other, BinaryFieldElement):
            wrapped = other
        else:
            wrapped = BinaryFieldElement(other)
        return self * wrapped.multiplicative_inverse()

    # Fermat's little theorem generalizes the inverse as x^(N-2) where N is the size of the smallest sb-field that self.value
    # lives in. Taking the highest set bit condition (ignoring leading zeros) then rounding to the nearest "tier" of the Binary
    # tower finds a value 'level' with is the bit-width of the smallest level containing the value. The bit size would be 2^level
    def multiplicative_inverse(self):

        if self.value == 0:
            raise ZeroDivisionError("No multiplicative inverse")

        first = self.value.bit_length()
        second = (first - 1).bit_length()
        level = 2 ** second
        N = 2 ** level
        return self ** (N - 2)

    # BinaryFieldElement -> raw bytes so it can be hashed
    # Careful about length though-?
    def to_bytes(self, length, byte_order):
        return self.value.to_bytes(length, byte_order)

    # Reverse of to_bytes- raw bytes -> BinaryFieldElement
    @classmethod
    def from_bytes(cls, b, byte_order):
        return cls(int.from_bytes(b, byte_order))

# Test driver for each method:
# def main():

#     a = BinaryFieldElement(3)
#     b = BinaryFieldElement(5)
#     c = BinaryFieldElement(13)

#     print(binary_multiplication(3, 3))         # Expected result is 2

#     print(a + b)                                # Expected result is 6  (3 XOR 5)

#     print(a * a)                                # Expected result is 2  (wraps binary_multiplication)

#     print(a ** 5)                               # Expected result is 2

#     print(c.multiplicative_inverse())           # Expected result is 8
#     print(c * c.multiplicative_inverse())       # Expected result is 1  (sanity check: x * mult_inv(x) == 1)

#     print(BinaryFieldElement(5).to_bytes(2, 'little'))              # Expected result is b'\x05\x00'
#     print(BinaryFieldElement.from_bytes(b'\x05\x00', 'little'))     # Expected result is 5

# main()