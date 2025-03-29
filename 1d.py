import itertools
from tqdm import tqdm
import math

# Parameters for the SHC: n cells, l levels, k bits.
# Here we assume n = 1, l = 32, k = 5.
n = 1
l = 32
k = 5

# We only allow up to this many thresholds per bit.
MAX_THRESHOLDS_NUM = 7


def decision(v, R):
    """Compute the decision bit for cell level v given threshold set R (a frozenset)."""
    b = 1
    for r in sorted(R):
        if r <= v:
            b ^= 1
        else:
            break
    return b


def is_balanced_threshold(R_tuple):
    """
    Check that for n = 1 (R_tuple = (R,)), as v runs from 0 to l-1,
    the decision vector has exactly l/2 ones and l/2 zeros.
    """
    R = R_tuple[0]
    counts = [0, 0]
    for v in range(l):
        counts[decision(v, R)] += 1
    return counts[0] == (l // 2) and counts[1] == (l // 2)


def is_minimal_threshold_set(R_tuple):
    """
    For n = 1, a threshold set R is minimal if removing any one threshold
    causes the decision vector to lose balance.
    """
    R = R_tuple[0]
    for r in R:
        reduced = frozenset(set(R) - {r})
        if is_balanced_threshold((reduced,)):
            return False
    return True


def generate_minimal_candidate_thresholds():
    """
    For n = 1 and l = 32, generate only those threshold sets that are:
      - Of odd size (since only an odd number of thresholds can yield balance).
      - Balanced: the decision vector (over 0..31) has 16 zeros and 16 ones.
      - Minimal: no threshold can be removed without breaking balance.
    """
    candidates = []
    # Only odd sizes can yield balance.
    for r in range(1, MAX_THRESHOLDS_NUM + 1, 2):
        # Iterate over combinations of r thresholds from {1,...,l-1}.
        for combo in itertools.combinations(range(1, l), r):
            R = frozenset(combo)
            R_tuple = (R,)
            if is_balanced_threshold(R_tuple) and is_minimal_threshold_set(R_tuple):
                candidates.append(R)
    return candidates


def generate_candidate_D(n):
    """
    Generate non-constant D functions for n cells, with symmetry reduction.
    For n = 1, there are only 2^2 = 4 possibilities; we exclude the constant ones.
    """
    size = 2 ** n
    seen = set()
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if all(b == bits[0] for b in bits):
            continue
        min_form = min(bits, tuple(1 - b for b in bits))
        if min_form not in seen:
            seen.add(min_form)
            candidates.append(bits)
    return candidates


def canonicalize_pair(R, D):
    """
    Produce a canonical representation for a candidate pair.
    R is a tuple (one per cell); we convert each frozenset into a sorted tuple.
    """
    canonical_R = tuple(tuple(sorted(r)) for r in R)
    return (canonical_R, tuple(D))


def decode_bit(candidate_pair, v):
    """Decode a bit for a given candidate pair and cell level v (for n cells)."""
    R, D = candidate_pair
    bits = [decision(v[j], R[j]) for j in range(n)]
    index = sum(bits[j] * (2 ** j) for j in range(n))
    return D[index]


def mapping(candidate_pairs_list):
    """
    Create a mapping from cell levels (as tuples) to a k-tuple of bits,
    using candidate pairs (one for each bit).
    """
    return {v: tuple(decode_bit(candidate_pairs_list[i], v) for i in range(k))
            for v in itertools.product(range(l), repeat=n)}


def is_surjective(mapping_dict):
    """Check that the mapping covers all 2^k possible k-bit outputs."""
    return len(set(mapping_dict.values())) == 2 ** k


def check_prop1(candidate_r, candidate_d):
    """
    For n = 1 and l = 32, check that applying decode_bit produces exactly
    16 zeros and 16 ones.
    """
    counts = [0, 0]
    for v in range(l):
        counts[decode_bit((candidate_r, candidate_d), (v,))] += 1
    return counts[0] == (l // 2) and counts[1] == (l // 2)


# For n = 1, generate candidate D functions.
candidate_Ds = generate_candidate_D(n)
# For n = 1 and l = 32, generate minimal candidate thresholds.
minimal_R_candidates = generate_minimal_candidate_thresholds()


def generate_candidate_pairs():
    """
    Generate candidate pairs (R, D) by pairing each minimal balanced threshold set R
    (wrapped as a one-tuple) with each candidate D, but only keeping those that satisfy
    the necessary condition (check_prop1) and are canonical (to remove duplicates).
    """
    seen_canonical = set()
    candidate_pairs = []
    # For n = 1, each candidate R is a frozenset from minimal_R_candidates.
    for R in minimal_R_candidates:
        R_tuple = (R,)
        for d in candidate_Ds:
            if not check_prop1(R_tuple, d):
                continue
            canon = canonicalize_pair(R_tuple, d)
            if canon not in seen_canonical:
                seen_canonical.add(canon)
                candidate_pairs.append((R_tuple, d))
    return candidate_pairs


def average_reads(candidate_solution):
    """
    For a candidate solution (a tuple of candidate pairs, one per bit), return the maximum
    number of thresholds (reads) used among the pairs.
    """
    return max(len(r) for R, _ in candidate_solution for r in R)


def main():
    print(f"Searching for SHC({n}, {l}, {k}) with aggressive preprocessing")
    candidate_pairs = generate_candidate_pairs()
    print(f"Candidate pairs after filtering: {len(candidate_pairs)}")

    total_combinations = math.comb(len(candidate_pairs), k)
    print(f"Total combinations to search: {total_combinations}")

    min_reads = float('inf')
    # Exhaustively search over all k-tuples of candidate pairs.
    for candidate_solution in tqdm(itertools.combinations(candidate_pairs, k), total=total_combinations,
                                   desc="Searching"):
        avg_reads = average_reads(candidate_solution)
        if avg_reads > min_reads:
            continue
        mapping_dict = mapping(candidate_solution)
        if is_surjective(mapping_dict):
            print("\nValid Solution Found:")
            for i, (R, D) in enumerate(candidate_solution, 1):
                # For readability, convert each frozenset to a sorted list.
                readable_R = [sorted(list(r)) for r in R]
                print(f"Bit {i}: R={readable_R}, D={D}")
            print(f"Average reads per bit: {avg_reads}")
            min_reads = avg_reads

    print(f"\nBest average reads per bit: {min_reads}")


if __name__ == '__main__':
    main()
