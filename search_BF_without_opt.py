import itertools
import math
import argparse
from multiprocessing import Pool, cpu_count, Manager
from tqdm import tqdm
from functools import lru_cache

MAX_THRESHOLDS_NUM = 8

# Global variables for multiprocessing
global_candidate_pairs = []
global_cell_values = []
global_n = 0
global_l = 0
global_k = 0


@lru_cache(maxsize=4096)
def decision(v, R):
    b = 1
    for r in sorted(R):
        if r > v:
            break
        b ^= 1
    return b


@lru_cache(maxsize=1024)
def is_balanced(R):
    mask = 0
    for v in global_cell_values:
        mask |= decision(v, R) << v
    expected_count = len(global_cell_values) // 2
    return bin(mask).count('1') == expected_count


@lru_cache(maxsize=1024)
def is_minimal(R):
    for r in R:
        reduced = tuple(sorted(set(R) - {r}))
        if reduced and is_balanced(reduced):
            return False
    return True


def generate_candidates(l):
    candidates = []
    for size in range(1, MAX_THRESHOLDS_NUM + 1, 2):
        for combo in itertools.combinations(range(1, l), size):
            if combo[0] > 2 or combo[-1] < (l - 2):
                continue
            if is_balanced(combo) and is_minimal(combo):
                candidates.append(combo)
    return candidates


def d_not_depend_on(d_mapping, n):
    return set()


def decode_bit(candidate_pair, v):
    R, D = candidate_pair
    result = decision(v, R)
    return D[result]


def test_combination(combination_index):
    combination = [global_candidate_pairs[i] for i in combination_index]
    mapping_dict = {}
    for v in itertools.product(range(global_l), repeat=global_n):
        output = tuple(decode_bit(combination[i], v[0]) for i in range(global_k))
        mapping_dict[v] = output
    if len(set(mapping_dict.values())) == 2 ** global_k:
        avg_reads = sum(len(R) for R, _ in combination) / global_k
        return True, combination_index, avg_reads
    return False, combination_index, float('inf')


def process_chunk(chunk):
    results = []
    for i in chunk:
        R, D = global_candidate_pairs[i]
        vec = sum((decision(v, R) ^ D[0]) << v for v in global_cell_values)
        unique_bits = bin(vec).count('1')
        if unique_bits >= (len(global_cell_values) // 3):
            results.append((i, vec))
    return len(chunk), results


def parallel_search_candidates(candidate_pairs, l):
    total = len(candidate_pairs)
    if total == 0:
        print("No candidate pairs to process.")
        return []

    chunk_size = max(1000, total // (cpu_count() * 10))
    chunks = [range(i, min(i + chunk_size, total)) for i in range(0, total, chunk_size)]
    valid = []

    with Manager() as manager:
        progress = manager.Value('i', 0)
        with tqdm(total=total, desc="Screening candidates") as pbar:
            with Pool(initializer=init_pool,
                      initargs=(progress, candidate_pairs, list(range(l)), 0, l, 0)) as pool:
                for processed, result in pool.imap_unordered(process_chunk, chunks):
                    pbar.update(processed)
                    valid.extend(result)

    valid.sort(key=lambda x: bin(x[1]).count('1'), reverse=True)
    return [idx for idx, _ in valid]


def process_combinations(combinations):
    results = []
    for combo in combinations:
        valid, indices, avg_reads = test_combination(combo)
        if valid:
            results.append((valid, indices, avg_reads))
    return len(combinations), results


def parallel_search_combinations(candidate_indices, candidate_pairs, n, l, k):
    k_combinations = list(itertools.combinations(candidate_indices, k))
    total = len(k_combinations)

    if total == 0:
        print("No combinations to process.")
        return None

    print(f"Testing {total} combinations of {k} candidate pairs")
    chunk_size = max(100, total // (cpu_count() * 10))
    chunks = [k_combinations[i:i + chunk_size] for i in range(0, total, chunk_size)]
    valid_solutions = []

    with Manager() as manager:
        progress = manager.Value('i', 0)
        with tqdm(total=total, desc="Searching for valid SHC") as pbar:
            with Pool(initializer=init_pool,
                      initargs=(progress, candidate_pairs, list(range(l)), n, l, k)) as pool:
                for processed, results in pool.imap_unordered(process_combinations, chunks):
                    pbar.update(processed)
                    valid_solutions.extend(results)

    if not valid_solutions:
        return None

    valid_solutions.sort(key=lambda x: x[2])
    best_solution = valid_solutions[0]
    return [candidate_pairs[i] for i in best_solution[1]]


def init_pool(progress_ref, candidate_pairs, cell_values, n, l, k):
    global progress, global_candidate_pairs, global_cell_values, global_n, global_l, global_k
    progress = progress_ref
    global_candidate_pairs = candidate_pairs
    global_cell_values = cell_values
    global_n = n
    global_l = l
    global_k = k


def optimize_candidate_pairs(R_candidates, D_mappings, n):
    optimized_pairs = []
    for r in R_candidates:
        for d in D_mappings:
            not_dependent_bits = d_not_depend_on(d, n)
            if not_dependent_bits:
                pass
            optimized_pairs.append((r, d))
    print(f"Optimized candidate pairs: {len(optimized_pairs)}")
    return optimized_pairs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--l", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    args = parser.parse_args()

    if args.n < 1 or args.l < 1 or args.k < 1:
        raise ValueError("n, l, and k must all be positive integers.")

    n = args.n
    l = args.l
    k = args.k

    global_cell_values = list(range(l))
    global_n = n
    global_l = l
    global_k = k

    R_candidates = generate_candidates(l)
    print(f"Generated {len(R_candidates)} R candidates")

    D_mappings = [(0, 1), (1, 0)]
    candidate_pairs = optimize_candidate_pairs(R_candidates, D_mappings, n)

    valid_indices = parallel_search_candidates(candidate_pairs, l)
    print(f"Found {len(valid_indices)} promising candidate pairs")

    solution = parallel_search_combinations(valid_indices, candidate_pairs, n, l, k)

    if solution:
        print("\nValid SHC Solution Found:")
        for i, (R, D) in enumerate(solution, 1):
            print(f"Bit {i}: R={set(R)}, D={D}")
        avg_reads = sum(len(R) for R, _ in solution) / k
        print(f"Average reads per bit: {avg_reads}")
        print("Optimization successful!")
    else:
        print("No valid SHC configuration found")
