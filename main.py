import itertools
from tqdm import tqdm
import math

# Parameters for the SHC: n cells, l levels, k bits.
n = 2  # number of cells per codeword
l = 4  # levels
k = 4  # number of bits to store

# Generate candidate thresholds dynamically based on l (1 to l-1 as singletons + empty)
candidate_thresholds = [frozenset()]
for r in range(1, l):
    candidate_thresholds.append(frozenset({r}))

# All possible assignments of thresholds for n cells.
candidate_Rs = list(itertools.product(candidate_thresholds, repeat=n))

# Generate all possible decoding functions D: {0,1}^n -> {0,1}.
# Exclude constant functions.
def generate_candidate_D(n):
    size = 2 ** n
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if not (all(b == bits[0] for b in bits)):# check if the function is not a constant
            candidates.append(bits)
    return candidates


candidate_Ds = generate_candidate_D(n)

# Candidate pairs are combinations of R and D
candidate_pairs = list(itertools.product(candidate_Rs, candidate_Ds)) #TODO: check redundant groups and switch them to empty


# Decision function (Algorithm 1)
def decision(v, R):
    b = 1
    for r in sorted(R):
        if r <= v:
            b ^= 1
        else:
            break
    return b


# Compute the output bit for a candidate pair and cell vector v
def decode_bit(candidate_pair, v):
    R, D = candidate_pair
    bits = [decision(v[j], R[j]) for j in range(n)]
    index = sum(bits[j] * (2 ** j) for j in range(n))
    return D[index]


# Compute the overall mapping for a set S of k pairs
def mapping(candidate_pairs_list):
    mapping_dict = {}
    for v in itertools.product(range(l), repeat=n):
        output = tuple(decode_bit(candidate_pairs_list[i], v) for i in range(k))
        mapping_dict[v] = output
    return mapping_dict


# Check if the mapping covers all possible outputs
def is_surjective(mapping_dict):
    outputs = set(mapping_dict.values())
    return len(outputs) == 2 ** k


# Calculate average number of reads per bit
def average_reads(candidate_solution):
    total_reads = 0
    for R, _ in candidate_solution:
        total_reads += max(len(r) for r in R)
    return total_reads / k


def main():
    print(f"Searching for SHC({n}, {l}, {k}) with thresholds: {[set(t) for t in candidate_thresholds]}")
    print(f"Possible R per cell: {len(candidate_thresholds)}, Total R combinations: {len(candidate_Rs)}")
    print(f"Possible D functions: {len(candidate_Ds)}, Candidate pairs per bit: {len(candidate_pairs)}")

    #TODO: make the check of all valid pairs and remove is_surjective. use memoization
    valid_solutions = []
    total_combinations = math.comb(len(candidate_pairs), k)

    min_reads = l
    # Iterate over all possible combinations of k pairs with a progress bar
    for candidate_solution in tqdm(itertools.combinations(candidate_pairs,k),total=total_combinations, desc="Searching"):
        avg_reads = average_reads(candidate_solution)
        if avg_reads > min_reads:
            continue
        mapping_dict = mapping(candidate_solution)
        if is_surjective(mapping_dict):#check that capacity is 1 and change the check
            valid_solutions.append((candidate_solution, mapping_dict, avg_reads))
            min_reads = avg_reads

    if valid_solutions:
        min_avg = min(sol[2] for sol in valid_solutions)
        best_solutions = [sol for sol in valid_solutions if sol[2] == min_avg]
        print(f"Found {len(valid_solutions)} valid SHCs. Best avg reads: {min_avg}")
        for sol in best_solutions:
            print("\nSolution:")
            for i, (R, D) in enumerate(sol[0], 1):
                print(f"Bit {i}: R={[set(r) for r in R]}, D={D}")
    else:
        print("No valid SHC found.")


if __name__ == '__main__':
    main()