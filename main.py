import itertools
from tqdm import tqdm

# Parameters for SHC: n cells, l levels, k bits
n = 2  # Number of cells per codeword
l = 4  # Levels (MLC)
k = 4  # Number of bits to store

# --- Step 1: Generate candidate thresholds (R) ---
# Each cell's thresholds are subsets of {1, 2, ..., l-1}
candidate_thresholds = [frozenset()]
for r in range(1, l):
    candidate_thresholds.append(frozenset({r}))
candidate_Rs = list(itertools.product(candidate_thresholds, repeat=n))

# --- Step 2: Generate candidate decoding rules (D) ---
# All non-constant functions {0,1}^n -> {0,1}
def generate_candidate_D(n):
    size = 2 ** n
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if not all(b == bits[0] for b in bits):  # Skip constant functions
            candidates.append(bits)
    return candidates

candidate_Ds = generate_candidate_D(n)

# --- Step 3: Filter candidate pairs (R, D) using Proposition 1 ---
def decision(v, R):
    """Algorithm 1: Compute decision bit for cell level `v` and thresholds `R`."""
    b = 1
    for r in sorted(R):
        if r <= v:
            b ^= 1
        else:
            break
    return b

def satisfies_proposition1(candidate_pair, n, l, k):
    """Check if a (R, D) pair meets the necessary condition."""
    R, D = candidate_pair
    required_min = 2 ** (k - 1)  # Minimum occurrences for 0/1
    outputs = []
    # Iterate over all possible cell level vectors (v1, v2, ..., vn)
    for v in itertools.product(range(l), repeat=n):
        bits = [decision(v[j], R[j]) for j in range(n)]  # Decision bits
        index = sum(bit * (2 ** j) for j, bit in enumerate(bits))  # Index into D
        outputs.append(D[index])
    # Check if 0/1 counts meet the requirement
    count0 = outputs.count(0)
    count1 = outputs.count(1)
    return count0 >= required_min and count1 >= required_min

# Filter candidate pairs
print("Filtering candidate pairs using Proposition 1...")
filtered_candidate_pairs = [
    pair for pair in tqdm(
        itertools.product(candidate_Rs, candidate_Ds),
        total=len(candidate_Rs) * len(candidate_Ds),
        desc="Filtering"
    ) if satisfies_proposition1(pair, n, l, k)
]

# --- Step 4: Exhaustive search over filtered candidates ---
def mapping(candidate_solution):
    """Compute the mapping for a set of k pairs."""
    mapping_dict = {}
    for v in itertools.product(range(l), repeat=n):
        output = tuple(
            D[sum(decision(v[j], R[j]) * (2 ** j) for j in range(n))]
            for (R, D) in candidate_solution
        )
        mapping_dict[v] = output
    return mapping_dict

def average_reads(candidate_solution):
    """Calculate average number of reads per bit."""
    total = 0
    for (R, _) in candidate_solution:
        total += max(len(r) for r in R)
    return total / k

def main():
    print(f"Filtered candidate pairs: {len(filtered_candidate_pairs)}")
    valid_solutions = []
    total_combinations = len(filtered_candidate_pairs) ** k

    # Iterate over all k-length combinations of filtered pairs
    for candidate_solution in tqdm(
        itertools.product(filtered_candidate_pairs, repeat=k),
        total=total_combinations,
        desc="Searching"
    ):
        mapping_dict = mapping(candidate_solution)
        outputs = set(mapping_dict.values())
        if len(outputs) == 2 ** k:  # Check surjectivity
            avg_reads = average_reads(candidate_solution)
            valid_solutions.append((candidate_solution, avg_reads))

    # Output results
    if valid_solutions:
        best_avg = min(sol[1] for sol in valid_solutions)
        best_solutions = [sol for sol in valid_solutions if sol[1] == best_avg]
        print(f"Found {len(valid_solutions)} valid SHCs. Best avg reads: {best_avg:.2f}")
        for sol in best_solutions[:1]:  # Print first best solution
            print("\nExample SHC configuration:")
            for i, (R, D) in enumerate(sol[0], 1):
                print(f"Bit {i}: R={[set(r) for r in R]}, D={D}")
    else:
        print("No valid SHC found.")

if __name__ == '__main__':
    main()