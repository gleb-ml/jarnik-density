"""Reproduce the computational tables for the Jarnik polygon density paper."""

import argparse
from collections import Counter
from functools import lru_cache
from math import inf, isqrt, prod


DEFAULT_LIMIT = 1_000_000
DEFAULT_MAX_K = 4
DEFAULT_SEED_LIMIT = 30_000
PAPER_CUTOFFS = (30_000, 100_000, 300_000, 1_000_000)


def allowed(limit, max_k):
    """List admissible odd radii and their split prime powers in order."""
    prime = bytearray(b"\1") * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            prime[p * p::p] = b"\0" * ((limit - p * p) // p + 1)
    primes = [p for p in range(5, limit + 1, 4) if prime[p]]

    def generate(m, factors, start):
        yield m, factors
        if len(factors) == max_k:
            return
        for j in range(start, len(primes)):
            p = primes[j]
            if m * p > limit:
                break
            q = p
            while m * q <= limit:
                yield from generate(m * q, factors + ((p, q),), j + 1)
                q *= p

    return sorted(generate(1, (), 0))


@lru_cache(None)
def gaussian_power(p, q):
    """Return coordinates of a Gaussian prime power of norm q = p**a."""
    nonresidue = 2
    while pow(nonresidue, (p - 1) // 2, p) == 1:
        nonresidue += 1
    a, b = p, pow(nonresidue, (p - 1) // 4, p)
    while b * b > p:
        a, b = b, a % b
    c = isqrt(p - b * b)
    x, y = 1, 0
    while q > 1:
        x, y = b * x - c * y, c * x + b * y
        q //= p
    return x, y


def representatives(factors):
    """Return one primitive vector from each opposite pair on a circle."""
    vectors = [(1, 0)]
    for p, q in factors:
        a, b = gaussian_power(p, q)
        vectors = [
            (a * x - c * y, c * x + a * y)
            for x, y in vectors
            for c in (b, -b)
        ]
    return [w for x, y in vectors for w in ((x, y), (-y, x))]


@lru_cache(96)
def small_sums(factors):
    """Map each attainable signed sum to its minimum number of vectors."""
    results = {(0, 0): 0}
    for a, b in representatives(factors):
        previous = results.copy()
        for (x, y), length in previous.items():
            for dest in ((x + a, y + b), (x - a, y - b)):
                results[dest] = min(results.get(dest, inf), length + 1)
    return results


def mixed_relation(factors, p, q, best):
    """Improve best by matching nonzero signed sums across the split."""
    length = 2 ** (len(factors) + 1)
    if q > length**2 * prod(r for _, r in factors):
        return best
    a, b = gaussian_power(p, q)
    if length <= 8:
        sums = small_sums(factors)
        for (x, y), s in sums.items():
            if s == 0 or s + 1 >= best:
                continue
            u, rem = divmod(a * x - b * y, q)
            if rem == 0:
                v = (b * x + a * y) // q
                best = min(
                    best, s + sums.get((b * v - a * u, -b * u - a * v), inf)
                )
                if best == 6:
                    break
        return best

    # Track reachable sums separately for each number of selected vectors.
    vectors = representatives(factors)
    bound = sum(abs(x) for x, _ in vectors)
    stride = 2 * bound + 1
    origin = bound * (stride + 1)
    cap = length if best == inf else min(length, best - 2)
    states = [1 << origin] + [0] * cap
    for index, (x, y) in enumerate(vectors, 1):
        shift = abs(x + stride * y)
        for s in range(min(index, cap), 0, -1):
            states[s] |= (states[s - 1] << shift) | (states[s - 1] >> shift)
    rows = [state.to_bytes(1 + 2 * origin // 8, "little") for state in states]

    def size(x, y, maximum):
        if max(abs(x), abs(y)) > bound:
            return inf
        cell, bit = divmod(origin + x + stride * y, 8)
        for s in range(1, min(cap, maximum) + 1):
            if rows[s][cell] & (1 << bit):
                return s
        return inf

    radius = (abs(a) + abs(b)) * bound // q
    for u in range(radius + 1):
        for v in range(1 if u == 0 else -radius, radius + 1):
            s = size(a * u + b * v, a * v - b * u, cap)
            if s < inf:
                r = size(
                    a * u - b * v,
                    a * v + b * u,
                    cap if best == inf else best - s - 1,
                )
                best = min(best, s + r)
                if best == 6:
                    return best
    return best


def run(
    limit=DEFAULT_LIMIT, max_k=DEFAULT_MAX_K, seed_limit=DEFAULT_SEED_LIMIT
):
    """Compute relation sizes, inherited sizes, and minimal radii."""
    data, cores = {}, []
    for m, factors in allowed(limit, max_k):
        parents = [data[m // q] for _, q in factors]
        inherited = min((v[1] for v in parents), default=inf)
        strict = min((v[2] for v in parents), default=inf)
        s, seed = inherited, strict
        k = len(factors)
        t = 2 ** (k + 1)
        if k >= 2 and s != 6:
            p, q = max(factors, key=lambda f: f[1])
            s = mixed_relation(
                tuple(f for f in factors if f != (p, q)), p, q, s
            )
        if s < inf and inherited == inf:
            cores.append((m, factors, t, s))
            if m <= seed_limit:
                seed = s
        data[m] = (t, s, seed, strict)
    return data, cores


def summarize(data, cores, limit):
    """Return weight sums, maximum even n, and both counts at a cutoff."""
    weights = Counter()
    for m, factors, t, _ in cores:
        if m <= limit:
            weights[len(factors)] += t / m
    pairs = marked = strict = 0
    for m, (t, s, seed, inherited) in data.items():
        odd, even = int(m <= limit), int(2 * m <= limit)
        pairs += (odd + even) * t
        marked += (odd + even) * max(0, t - s + 1)
        strict += odd * max(0, t - inherited + 1) + even * max(0, t - seed + 1)
    return dict(sorted(weights.items())), 2 * pairs, marked, strict


def format_factors(factors):
    """Format the stored (prime, prime power) pairs for display."""
    terms = []
    for prime, power in factors:
        exponent = 0
        while power > 1:
            power //= prime
            exponent += 1
        terms.append(str(prime) if exponent == 1 else f"{prime}^{exponent}")
    return " * ".join(terms)


def print_table(headers, rows, left_aligned=()):
    """Print a table with column widths fitted to its headers and rows."""
    widths = [max(map(len, column)) for column in zip(headers, *rows)]
    for row in [headers, ["-" * width for width in widths], *rows]:
        line = "  ".join(
            value.ljust(width) if index in left_aligned else value.rjust(width)
            for index, (value, width) in enumerate(zip(row, widths))
        )
        print(line.rstrip())


def main():
    """Parse options, run the enumeration, and print the paper tables."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--limit", type=int, default=DEFAULT_LIMIT,
        help="largest squared radius to enumerate",
    )
    parser.add_argument(
        "--max-k", type=int, default=DEFAULT_MAX_K,
        help="maximum number of distinct prime factors congruent to 1 mod 4",
    )
    parser.add_argument(
        "--seed-limit", type=int, default=DEFAULT_SEED_LIMIT,
        help="largest minimal radius used to certify strict descendants",
    )
    parser.add_argument(
        "--show-cores", type=int, default=32,
        help="number of minimal radii to display",
    )
    args = parser.parse_args()
    if args.limit < 1:
        parser.error("--limit must be at least 1")
    for name in ("max_k", "seed_limit", "show_cores"):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be nonnegative")

    data, cores = run(args.limit, args.max_k, args.seed_limit)
    if args.show_cores:
        print("First minimal radii")
        rows = [
            [f"{m:,}", format_factors(factors), str(t), str(s)]
            for m, factors, t, s in cores[:args.show_cores]
        ]
        print_table(["m", "Factorization", "t", "s"], rows, left_aligned=(1,))
        print()
    print(f"Total minimal radii: {len(cores):,}")

    cutoffs = {cutoff for cutoff in PAPER_CUTOFFS if cutoff <= args.limit}
    cutoffs.add(args.limit)
    summaries = [
        (cutoff, summarize(data, cores, cutoff)) for cutoff in sorted(cutoffs)
    ]
    weight_indices = sorted({
        k for _, (weights, _, _, _) in summaries for k in weights
    })
    headers = [
        "Cutoff m", "Max even n", "Even n >= 4", "Marked", "Marked %",
        "Strict inh.", "Strict inh. %",
    ] + [f"S_{k}" for k in weight_indices]
    rows = []
    for cutoff, (weights, max_n, marked, strict) in summaries:
        even_count = max_n // 2 - 1
        rows.append([
            f"{cutoff:,}",
            f"{max_n:,}",
            f"{even_count:,}",
            f"{marked:,}",
            f"{marked / even_count:.4%}",
            f"{strict:,}",
            f"{strict / even_count:.4%}",
        ] + [f"{weights.get(k, 0):.6f}" for k in weight_indices])
    print()
    print("Summary by squared radius")
    print_table(headers, rows)


if __name__ == "__main__":
    main()
