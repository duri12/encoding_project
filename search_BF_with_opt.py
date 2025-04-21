import itertools
import random
from tqdm import tqdm
import math
import argparse

MAX_THRESHOLDS_NUM = 18
MIN_THRESHOLDS_NUM = 1
BATCH_SIZE = 10000


def generate_candidate_thresholds(l):
    s = list(range(1, l))
    for r in range(MIN_THRESHOLDS_NUM, MAX_THRESHOLDS_NUM + 1):
        thresholds = []
        for batch_start in range(0, math.comb(len(s), r), BATCH_SIZE):
            # Generate a batch of combinations of r thresholds from s
            thresholds += [frozenset(subset) for subset in
                           itertools.islice(itertools.combinations(s, r), batch_start, batch_start + BATCH_SIZE)]
        random.shuffle(thresholds)
        yield thresholds


def generate_candidate_D(n):
    size = 2 ** n
    seen = set()
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        # Skip constant vectors
        if all(b == bits[0] for b in bits):
            continue
        # Use symmetry to avoid duplicates (bitwise negation)
        min_form = min(bits, tuple(1 - b for b in bits))
        if min_form not in seen:
            seen.add(min_form)
            candidates.append(bits)
    return candidates


# New optimization: identify bits that D doesn't depend on
def d_not_depend_on(candidate_d, n):
    dependent_bits = set()
    for bit in range(n):
        for i in range(len(candidate_d)):
            j = i ^ (1 << bit)  # Flip the `bit`-th position
            if j > i and candidate_d[i] != candidate_d[j]:
                dependent_bits.add(bit)
    not_dependent_bits = set(range(n)) - dependent_bits
    return not_dependent_bits


def decision(v, R):
    b = 1
    for r in sorted(R):
        if r <= v:
            b ^= 1
        else:
            break
    return b


def decode_bit(candidate_pair, v):
    R, D = candidate_pair
    # Compute each bit using the threshold sets
    bits = [decision(v[j], R[j]) for j in range(len(v))]
    # Convert bit vector to an index for D
    index = sum(bits[j] * (2 ** j) for j in range(len(v)))
    return D[index]


def mapping(candidate_pairs_list, n, l):
    # Create a mapping from input vectors to output k-bit vectors
    return {v: tuple(decode_bit(candidate_pairs_list[i], v) for i in range(len(candidate_pairs_list))) for v in
            itertools.product(range(l), repeat=n)}


def is_surjective(mapping_dict, k):
    # Check if all 2^k possible outputs are covered
    return len(set(mapping_dict.values())) == 2 ** k


def check_prop1(candidate_r, candidate_d, n, l, k):
    counters = [0, 0]
    THRESHOLD = 2 ** (k - 1)
    for v in itertools.product(range(l), repeat=n):
        out = decode_bit((candidate_r, candidate_d), v)
        counters[out] += 1
        # Check if both outputs appear at least 2^(k-1) times
        if counters[0] >= THRESHOLD and counters[1] >= THRESHOLD:
            return True
    return False


def generate_candidate_pairs(n, l, k, candidate_Ds, candidate_thresholds_batches):
    seen_canonical = set()
    candidate_pairs = []

    # Pre-calculate all possible Rs
    all_Rs = list(itertools.product(list(itertools.chain(*candidate_thresholds_batches)), repeat=n))

    random.shuffle(candidate_Ds)
    for d in candidate_Ds:
        # Optimization: identify bits that D doesn't depend on
        not_dependent_bits = d_not_depend_on(d, n)
        for r in all_Rs:
            # Optimization: remove unnecessary thresholds for bits D doesn't depend on
            new_r = list(r)
            for bit in not_dependent_bits:
                new_r[bit] = frozenset()
            new_r = tuple(new_r)
            if check_prop1(new_r, d, n, l, k):
                canonical_form = (new_r, d)
                if canonical_form not in seen_canonical:
                    seen_canonical.add(canonical_form)
                    candidate_pairs.append((new_r, d))

    random.shuffle(candidate_pairs)
    return candidate_pairs


def average_reads(candidate_solution):
    # Maximum number of thresholds used for any individual read operation
    return max(len(r) for R, _ in candidate_solution for r in R)


def main():
    parser = argparse.ArgumentParser(description="Search for SHC(n, l, k) codes.")
    parser.add_argument("--n", type=int, required=True, help="Number of dimensions")
    parser.add_argument("--l", type=int, required=True, help="Alphabet size per dimension")
    parser.add_argument("--k", type=int, required=True, help="Number of output bits")
    args = parser.parse_args()

    n = args.n
    l = args.l
    k = args.k

    if n < 1 or l < 2 or k < 1:
        raise ValueError("n must be >= 1, l >= 2, k >= 1")

    print(f"Searching for SHC({n}, {l}, {k})")

    candidate_Ds = generate_candidate_D(n)
    print(f"Generated {len(candidate_Ds)} candidate D functions")

    # Collect all threshold batches
    all_threshold_batches = list(generate_candidate_thresholds(l))

    candidate_pairs = generate_candidate_pairs(n, l, k, candidate_Ds, all_threshold_batches)
    print(f"After applying optimization: {len(candidate_pairs)} candidate pairs")
    print(f"Before optimization (estimated): {len(candidate_Ds) * (l ** n)}")

    min_reads = float('inf')
    total_combinations = math.comb(len(candidate_pairs), k)

    for candidate_solution in tqdm(itertools.combinations(candidate_pairs, k), total=total_combinations,
                                   desc="Searching"):
        avg_reads = average_reads(candidate_solution)
        if avg_reads > min_reads:
            continue
        mapping_dict = mapping(candidate_solution, n, l)
        if is_surjective(mapping_dict, k):
            print("\nValid Solution Found:")
            for i, (R, D) in enumerate(candidate_solution, 1):
                print(f"Bit {i}: R={[set(r) for r in R]}, D={D}")
            print(f"Average reads per bit: {avg_reads}")
            min_reads = avg_reads

    print(f"\nBest average reads per bit: {min_reads}")


if __name__ == '__main__':
    main()