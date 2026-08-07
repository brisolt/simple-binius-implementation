# Author: Brian Soltani
# Date: 08/06/2026
# Description:
#   Utilities class for my implementation of a simple Binius, consisting of helper methods.

import math
from binary_fields import BinaryFieldElement

# Only works for powers of two.
def log_2(x):
    if x > 0 and (x & (x - 1)) == 0: # Checks if x is a power of two.
        return int(math.log2(x))
    return None

def recursive_evaluations(evals, point):
    if len(point) == 0:
        return evals[0]
    half = len(evals) // 2
    top = recursive_evaluations(evals[:half], point[:-1])
    bottom = recursive_evaluations(evals[half:], point[:-1])

    return top + (bottom - top) * point[-1]

def multilinear_polynomial_evaluation(evals, point):
    assert len(evals) == 2 ** len(point)
    return recursive_evaluations(evals, point)


def evaluation_tensor_product(point):
    result = [1]
    for coordinate in point:
        low = [(1 - coordinate) * v for v in result]
        high = [coordinate * v for v in result]
        result = low + high
    return result

# OD: Evaluates a UNIVERIATE polynomial with coefficients listed in polynomials[] at a single point, point.
# Relies on BinaryFieldElement.__radd__ to handle the first addition.
# The above comment could be optimized by  "extending" from binary_fields.py

# Curr: evaluates a univariate polynomial whose coeff are listed in ascending degree order, at a single point, using Horner's rule
def evaluate_polynomial_at_point(polynomials, point):
    # evaluation = 0
    # for i in range(len(polynomials)):
    #     evaluation += polynomials[i] * (point ** i)
    # return evaluation
    result = BinaryFieldElement(0)
    for coeff in reversed(polynomials):
        result = result * point + coeff
    return result # homer's rule

# Pretty straightforward: corresponding terms get multiplied and added to a final list.
def multiply_polynomials(p1, p2):
    length = len(p1) + len(p2) - 1
    result = [BinaryFieldElement(0) for _ in range(length)]
    for i, p1_coefficients in enumerate(p1):
        for j, p2_coefficients in enumerate(p2):
            result[i + j] += p1_coefficients * p2_coefficients
    return result

def divide_polynomials(p1, p2):
    assert len(p2) > 0 and len(p1) >= len(p2)
    r = [c if isinstance(c, BinaryFieldElement) else BinaryFieldElement(c) for c in p1]
    q = [BinaryFieldElement(0) for _ in range(len(p1) - len(p2) + 1)]
    for i in range(len(p1) - len(p2), -1, -1):
        factor = r[i + len(p2) - 1] / p2[-1]
        q[i] = factor
        for j, coefficient in enumerate(p2):
            r[i + j] -= coefficient * factor
    return q

# Builds a polynomial that equals 1 at point and zero at all other points in the range(size).
def compute_lagrange_polynomial(size, point):
    # T/S: kills at source rather than relying on future caller to wrap input
    if not isinstance(point, BinaryFieldElement):
        point = BinaryFieldElement(point)

    polynomial = [BinaryFieldElement(1)]
    denominator = BinaryFieldElement(1)

    for i in range(0, size):
        if i == point:
            continue
        polynomial = multiply_polynomials(polynomial, [i, 1])
        denominator = denominator * (point - i)

    return [coefficient / denominator for coefficient in polynomial]

def extend(values, expansion_factor):
    n = len(values)
    lagranges = [compute_lagrange_polynomial(n, BinaryFieldElement(i)) for i in range(n)] ### T/S: Forces point to always be a real field element.

    output = values[:]
    for x in range(n, expansion_factor * n):
        total = 0
        for v, L in zip(values, lagranges):
            total += v * evaluate_polynomial_at_point(L, BinaryFieldElement(x))
        output.append(total)
    return output


# def main():

#     poly = compute_lagrange_polynomial(4, BinaryFieldElement(2))
#     print(poly)  # Expected: [0, 3, 2, 1]

# main()