# Author: Brian Soltani
# Date: 08/07/2026
# Description:
#   Utilities class for my implementation of a simple Binius, consisting of helper methods.

import math

from binary_fields import BinaryFieldElement


def log_2(x: int) -> int | None:
    """Returns log2(x) if x is a power of two, else None."""
    if x > 0 and (x & (x - 1)) == 0:
        return int(math.log2(x))
    return None


def recursive_evaluations(evals: list[BinaryFieldElement], point: list[BinaryFieldElement]) -> BinaryFieldElement:
    """
    Evaluates a multilinear polynomial, given as its values on the boolean hypercube, at point.
    Recursively splits evals in half on the last coordinate of point and linearly interpolates
    between the two halves.
    """
    if len(point) == 0:
        return evals[0]
    half = len(evals) // 2
    top = recursive_evaluations(evals[:half], point[:-1])
    bottom = recursive_evaluations(evals[half:], point[:-1])

    return top + (bottom - top) * point[-1]


def multilinear_polynomial_evaluation(evals: list[BinaryFieldElement], point: list[BinaryFieldElement]) -> BinaryFieldElement:
    """Entry point for recursive_evaluations: checks evals has the expected length, 2^len(point), before recursing."""
    assert len(evals) == 2 ** len(point)
    return recursive_evaluations(evals, point)


def evaluation_tensor_product(point: list[BinaryFieldElement]) -> list[BinaryFieldElement]:
    """
    Builds the tensor-product basis expansion of point: the length-2^len(point) vector whose
    entries are, for every subset of coordinates, the product of that coordinate (if included)
    or its complement (if not). Used to fold a multilinear polynomial's evaluations down to a
    single value via a dot product.
    """
    result = [1]
    for coordinate in point:
        low = [(1 - coordinate) * value for value in result]
        high = [coordinate * value for value in result]
        result = low + high
    return result


def evaluate_polynomial_at_point(coefficients: list[BinaryFieldElement], point: BinaryFieldElement) -> BinaryFieldElement:
    """Evaluates a univariate polynomial, given as coefficients in ascending degree order, at point, using Horner's rule."""
    result = BinaryFieldElement(0)
    for coefficient in reversed(coefficients):
        result = result * point + coefficient
    return result


def multiply_polynomials(poly1: list[BinaryFieldElement], poly2: list[BinaryFieldElement]) -> list[BinaryFieldElement]:
    """Multiplies two polynomials, given as coefficient lists in ascending degree order, via convolution."""
    length = len(poly1) + len(poly2) - 1
    result = [BinaryFieldElement(0) for _ in range(length)]
    for i, coefficient1 in enumerate(poly1):
        for j, coefficient2 in enumerate(poly2):
            result[i + j] += coefficient1 * coefficient2
    return result


def divide_polynomials(dividend: list[BinaryFieldElement], divisor: list[BinaryFieldElement]) -> list[BinaryFieldElement]:
    """Polynomial long division: divides dividend by divisor and returns the quotient's coefficients, ascending degree order."""
    assert len(divisor) > 0 and len(dividend) >= len(divisor)
    remainder = [coeff if isinstance(coeff, BinaryFieldElement) else BinaryFieldElement(coeff) for coeff in dividend]
    quotient = [BinaryFieldElement(0) for _ in range(len(dividend) - len(divisor) + 1)]
    for i in range(len(dividend) - len(divisor), -1, -1):
        factor = remainder[i + len(divisor) - 1] / divisor[-1]
        quotient[i] = factor
        for j, coefficient in enumerate(divisor):
            remainder[i + j] -= coefficient * factor
    return quotient


def compute_lagrange_polynomial(size: int, point: BinaryFieldElement | int) -> list[BinaryFieldElement]:
    """
    Builds the coefficients (ascending degree) of the unique degree-(size-1) polynomial that
    equals 1 at x = point and 0 at every other x in range(size).
    """
    # accept a raw int here so callers don't have to wrap it themselves
    if not isinstance(point, BinaryFieldElement):
        point = BinaryFieldElement(point)

    polynomial = [BinaryFieldElement(1)]
    denominator = BinaryFieldElement(1)

    for i in range(size):
        if i == point:
            continue
        polynomial = multiply_polynomials(polynomial, [i, 1])
        denominator = denominator * (point - i)

    return [coefficient / denominator for coefficient in polynomial]


def extend(values: list[BinaryFieldElement], expansion_factor: int) -> list[BinaryFieldElement]:
    """
    Reed-Solomon-style extension: interpolates the unique polynomial passing through values
    (treated as evaluations at x = 0..n-1), evaluates it at the next (expansion_factor - 1) * n
    points, and appends those to the original values.
    """
    n = len(values)
    lagranges = [compute_lagrange_polynomial(n, BinaryFieldElement(i)) for i in range(n)]

    output = values[:]
    for x in range(n, expansion_factor * n):
        total = 0
        for value, lagrange_poly in zip(values, lagranges):
            total += value * evaluate_polynomial_at_point(lagrange_poly, BinaryFieldElement(x))
        output.append(total)
    return output
