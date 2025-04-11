import itertools
import math
import argparse
from multiprocessing import Pool, cpu_count, Manager
from tqdm import tqdm
from functools import lru_cache

MAX_THRESHOLDS_NUM = 8

@lru_cache(maxsize=4096)
def decision(v, R):
    b = 1
    for r in sorted(R):
        if r > v:
            break
        b ^= 1
    return b

def is_balanced(R):
    mask = 0
    for v in CELL_VALUES:
        mask |= decision(v, R) << v
    return bin(mask).count('1') == 16

@lru_cache(maxsize=1024)
def is_minimal(R):
    for r in R:
        reduced = tuple(sorted(set(R) - {r}))
        if reduced and is_balanced(reduced):
            return False
    return True

def generate_candidates():
    candidates = []
    for size in range(1, MAX_THRESHOLDS_NUM + 1, 2):
        for combo in itertools.combinations(range(1, 31), size):
            if combo[0] > 2 or combo[-1] < 30:
                continue
            if is_balanced(combo) and is_minimal(combo):
                candidates.append(combo)
    return candidates

def process_chunk(chunk):
    results = []
    for i in chunk:
        R, D = candidate_pairs[i]
        vec = sum((decision(v, R) ^ D[0]) << v for v in CELL_VALUES)
        unique_bits = bin(vec).count('1')
        if unique_bits >= 10:
            results.append((i, vec))
    return len(chunk), results

def parallel_search():
    total = len(candidate_pairs)
    chunk_size = max(1000, total // (cpu_count() * 10))
    chunks = [range(i, min(i + chunk_size, total)) for i in range(0, total, chunk_size)]
    valid = []
    with Manager() as manager:
        progress = manager.Value('i', 0)
        with tqdm(total=total, desc="Screening") as pbar:
            with Pool(initializer=init_pool, initargs=(progress,)) as pool:
                for result in pool.imap_unordered(process_chunk, chunks):
                    pbar.update(result[0])
                    valid.extend(result[1])
        valid.sort(key=lambda x: bin(x[1]).count('1'), reverse=True)

def init_pool(progress_ref):
    global progress
    progress = progress_ref

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

    CELL_VALUES = list(range(l))
    R_candidates = generate_candidates()
    D_mappings = [(0, 1), (1, 0)]
    candidate_pairs = [(r, d) for r in R_candidates for d in D_mappings]

    if parallel_search():
        print("Optimization successful!")
    else:
        print("No valid configuration found")
