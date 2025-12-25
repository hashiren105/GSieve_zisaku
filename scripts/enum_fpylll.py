#!/usr/bin/env python3.10
"""Use fpylll enumeration to target a vector near the Gaussian heuristic."""

from __future__ import annotations

import argparse
import math
from typing import Iterable, List, Sequence

from fpylll import GSO, IntegerMatrix, LLL
from fpylll.fplll.enumeration import Enumeration, EnumerationError


def parse_basis(path: str) -> List[List[int]]:
    rows: List[List[int]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if not stripped.startswith("["):
                break
            content = stripped.lstrip("[").rstrip("]")
            if not content:
                break
            rows.append([int(entry) for entry in content.split()])
    if not rows:
        raise ValueError(f"no basis rows found in {path}")
    length = len(rows[0])
    for idx, row in enumerate(rows):
        if len(row) != length:
            raise ValueError(
                f"row {idx} has {len(row)} entries but row 0 has {length}"
            )
    return rows


def load_integer_matrix(rows: Sequence[Sequence[int]]) -> IntegerMatrix:
    matrix = IntegerMatrix(len(rows), len(rows[0]))
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            matrix[i, j] = value
    return matrix


def copy_basis(matrix: IntegerMatrix) -> List[List[int]]:
    rows = []
    for i in range(matrix.nrows):
        rows.append([int(matrix[i, j]) for j in range(matrix.ncols)])
    return rows


def compute_gaussian_heuristic(matrix: GSO.MatGSO) -> float:
    dimension = matrix.d
    log_det = matrix.get_log_det(0, dimension)
    log_volume = log_det / 2
    return math.exp(log_volume / dimension) * math.sqrt(
        dimension / (2 * math.pi * math.e)
    )


def combine_vector(
    coeffs: Sequence[float], basis: Sequence[Sequence[int]]
) -> List[int]:
    dimension = min(len(coeffs), len(basis))
    vector = [0] * len(basis[0])
    for idx in range(dimension):
        coef = int(round(coeffs[idx]))
        if coef == 0:
            continue
        row = basis[idx]
        for jdx, base_val in enumerate(row):
            vector[jdx] += coef * base_val
    return vector


def search_with_factor(
    enum: Enumeration,
    basis: Sequence[Sequence[int]],
    gh_length: float,
    factor: float,
    report_limit: int,
) -> tuple[int, List[tuple[float, List[int], float]]]:
    radius = (factor * gh_length) ** 2
    mantissa, expo = math.frexp(radius)
    if mantissa == 0.0:
        return 0, []
    try:
        solutions = enum.enumerate(0, len(basis), mantissa, expo)
    except EnumerationError:
        return 0, []
    results: List[tuple[float, List[int], float]] = []
    for length, coeffs in solutions:
        vector = combine_vector(coeffs, basis)
        if all(v == 0 for v in vector):
            continue
        norm_sq = sum(v * v for v in vector)
        norm = math.sqrt(norm_sq)
        results.append((norm, vector, length))
    results.sort(key=lambda entry: entry[0])
    return len(results), results[:report_limit]


def build_factor_sequence(start: float, stop: float, step: float) -> List[float]:
    factors = []
    current = start
    while current <= stop + step * 1e-9:
        factors.append(round(current, 9))
        current += step
    return factors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run fpylll enumeration on a basis and compare to the GH length."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="input/ideallatticeindex23seed0.txt",
        help="path to the basis file",
    )
    parser.add_argument(
        "--factor",
        type=float,
        default=1.05,
        help="starting multiple of GH to enumerate (default: 1.05)",
    )
    parser.add_argument(
        "--max-factor",
        type=float,
        help="maximum multiple to try (defaults to --factor)",
    )
    parser.add_argument(
        "--step",
        type=float,
        default=0.01,
        help="step size between factors when searching (default: 0.01)",
    )
    parser.add_argument(
        "--show",
        type=int,
        default=2,
        help="how many solution vectors to print per factor (default: 2)",
    )
    args = parser.parse_args()

    max_factor = args.max_factor if args.max_factor is not None else args.factor
    if max_factor < args.factor:
        parser.error("--max-factor must be >= --factor")
    if args.step <= 0:
        parser.error("--step must be positive")

    raw_basis = parse_basis(args.input)
    matrix = load_integer_matrix(raw_basis)
    LLL.reduction(matrix)
    reduced_basis = copy_basis(matrix)

    gso_matrix = GSO.Mat(matrix)
    gso_matrix.update_gso()
    gh_length = compute_gaussian_heuristic(gso_matrix)

    print(f"input: {args.input}")
    print(f"dimension: {matrix.nrows}")
    print(f"GH estimate: {gh_length:.6g}")

    enum = Enumeration(gso_matrix)
    searched = False
    for factor in build_factor_sequence(args.factor, max_factor, args.step):
        total, solutions = search_with_factor(
            enum, reduced_basis, gh_length, factor, args.show
        )
        if not total:
            print(f"factor {factor:.4f}: no short vectors found")
            continue
        searched = True
        print(
            f"factor {factor:.4f} (~{factor * gh_length:.6g} target length) => {total} solution(s)"
        )
        for idx, (norm, vector, raw_length) in enumerate(solutions, start=1):
            ratio = norm / gh_length
            print(
                f"  solution {idx}: norm={norm:.6g} (~{ratio:.4f}×GH, raw length={raw_length})"
            )
            print(f"    vector: {vector}")
        break
    if not searched:
        print(f"no vectors found up to factor {max_factor:.3f}")


if __name__ == "__main__":
    main()
