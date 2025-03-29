import itertools
import random
from tqdm import tqdm
import math

# Parameters for the SHC: n cells, l levels, k bits.
n = 1  # number of cells per codeword
l = 32  # levels
k = 5  # number of bits to store

MAX_THRESHOLDS_NUM = 7  # Max number of thresholds per bit
MIN_THRESHOLDS_NUM = 6  # Minimum number of thresholds per bit (new parameter)
BATCH_SIZE = 10000  # Process threshold combinations in batches


def generate_candidate_thresholds():
    """Generate threshold sets in batches to avoid MemoryError."""
    s = list(range(1, l))
    for r in range(MIN_THRESHOLDS_NUM, MAX_THRESHOLDS_NUM + 1):  # Start from MIN_THRESHOLDS_NUM
        thresholds = []
        for batch_start in range(0, math.comb(len(s), r), BATCH_SIZE):
            thresholds += [frozenset(subset) for subset in
                           itertools.islice(itertools.combinations(s, r), batch_start, batch_start + BATCH_SIZE)]
        random.shuffle(thresholds)  # Shuffle the thresholds to introduce randomness
        yield thresholds


def generate_candidate_D(n):
    """Generate non-constant D functions, considering symmetries."""
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
    bits = [decision(v[j], R[j]) for j in range(n)]
    index = sum(bits[j] * (2 ** j) for j in range(n))
    return D[index]


def mapping(candidate_pairs_list):
    return {v: tuple(decode_bit(candidate_pairs_list[i], v) for i in range(k))
            for v in itertools.product(range(l), repeat=n)}


def is_surjective(mapping_dict):
    return len(set(mapping_dict.values())) == 2 ** k


def check_prop1(candidate_r, candidate_d):
    counters = [0, 0]
    THRESHOLD = pow(2, k - 1)
    for v in itertools.product(range(l), repeat=n):
        out = decode_bit((candidate_r, candidate_d), v)
        counters[out] += 1
        if counters[0] >= THRESHOLD and counters[1] >= THRESHOLD:
            return True
    return False


candidate_Ds = generate_candidate_D(n)


def generate_candidate_pairs(n, candidate_Ds):
    seen_canonical = set()
    candidate_pairs = []

    for candidate_thresholds_batch in generate_candidate_thresholds():
        random.shuffle(candidate_Ds)  # Shuffle the D functions to add randomness
        for d in candidate_Ds:
            for r in itertools.product(candidate_thresholds_batch, repeat=n):
                if check_prop1(r, d):
                    canonical_form = (tuple(r), tuple(d))
                    if canonical_form not in seen_canonical:
                        seen_canonical.add(canonical_form)
                        candidate_pairs.append((r, d))
    random.shuffle(candidate_pairs)  # Shuffle candidate pairs to randomize the selection order
    return candidate_pairs


def average_reads(candidate_solution):
    return max(len(r) for R, _ in candidate_solution for r in R)


def main():
    print(f"Searching for SHC({n}, {l}, {k})")
    candidate_pairs = generate_candidate_pairs(n, candidate_Ds)
    print(f"Candidate pairs: {len(candidate_pairs)}")

    min_reads = float('inf')
    total_combinations = min(90000000, math.comb(len(candidate_pairs), k))  # Limit search space

    for candidate_solution in tqdm(itertools.combinations(candidate_pairs, k), total=total_combinations,
                                   desc="Searching"):
        avg_reads = average_reads(candidate_solution)
        if avg_reads > min_reads:
            continue
        mapping_dict = mapping(candidate_solution)
        if is_surjective(mapping_dict):
            print("\nValid Solution Found:")
            for i, (R, D) in enumerate(candidate_solution, 1):
                print(f"Bit {i}: R={[set(r) for r in R]}, D={D}")
            print(f"Average reads per bit: {avg_reads}")
            min_reads = avg_reads

    print(f"\nBest average reads per bit: {min_reads}")


if __name__ == '__main__':
    main()
