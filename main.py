import itertools

# Parameters for the SHC: n cells, l levels, k bits.
n = 2  # number of cells per codeword
l = 4  # levels: 0, 1, 2
k = 4  # number of bits to store

# For each cell, allowed read threshold sets.
# We restrict each R(i)_j to either be empty (no threshold) or a singleton from {1, 2}.
candidate_thresholds = [frozenset(), frozenset({1}), frozenset({2})]

# All possible assignments of thresholds for n cells.
# Each candidate R is a tuple of n sets.
candidate_Rs = list(itertools.product(candidate_thresholds, repeat=n))


# Generate all possible decoding functions D: {0,1}^n -> {0,1}.
# Represent D as a tuple (lookup table) of length 2^n.
# Exclude constant functions.
def generate_candidate_D(n):
    size = 2 ** n
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if all(b == bits[0] for b in bits):
            continue  # skip constant functions
        candidates.append(bits)
    return candidates


candidate_Ds = generate_candidate_D(n)

# A candidate pair for one bit is (R, D) where R is a tuple of length n (one per cell)
candidate_pairs = list(itertools.product(candidate_Rs, candidate_Ds))


# Decision function as defined in Algorithm 1.
# For a given cell level v and a set of thresholds R, return the decision bit.
def decision(v, R):
    b = 1
    for r in sorted(R):
        if r <= v:
            b ^= 1  # XOR with 1
        else:
            break
    return b


# Given a candidate pair (R, D) and an input cell vector v,
# compute the output bit.
def decode_bit(candidate_pair, v):
    R, D = candidate_pair  # R is a tuple of n threshold sets, D is the lookup table.
    bits = [decision(v[j], R[j]) for j in range(n)]
    # Convert binary vector to an integer (with cell 0 as LSB)
    index = sum(bits[j] * (2 ** j) for j in range(n))
    return D[index]


# Compute the overall mapping f_S for a set S of k candidate pairs.
def mapping(candidate_pairs_list):
    mapping_dict = {}
    for v in itertools.product(range(l), repeat=n):
        output = tuple(decode_bit(candidate_pairs_list[i], v) for i in range(k))
        mapping_dict[v] = output
    return mapping_dict


# Check if the mapping is surjective: does its image cover all 2^k possible output vectors?
def is_surjective(mapping_dict):
    outputs = set(mapping_dict.values())
    return len(outputs) == 2 ** k


# Compute the "read cost" for one candidate pair.
# This is defined as the maximum number of thresholds used among the n cells.
def reads_for_bit(R):
    return max(len(r) for r in R)


# Compute the average number of reads for a candidate solution (list of k candidate pairs)
def average_reads(candidate_solution):
    total_reads = 0
    for candidate_pair in candidate_solution:
        R, _ = candidate_pair
        total_reads += reads_for_bit(R)
    return total_reads / k


def main():
    print(f"Searching for an SHC({n}, {l}, {k}) using a simple brute-force approach.")
    print(f"Number of candidate pairs per bit: {len(candidate_pairs)}")

    valid_solutions = []

    # Iterate over every combination of candidate pairs for the k bits.
    for candidate_solution in itertools.product(candidate_pairs, repeat=k):
        mapping_dict = mapping(candidate_solution)
        if is_surjective(mapping_dict):
            avg_reads = average_reads(candidate_solution)
            valid_solutions.append((candidate_solution, mapping_dict, avg_reads))

    if valid_solutions:
        # Determine the minimum average reads among all valid solutions.
        min_avg = min(sol[2] for sol in valid_solutions)
        best_solutions = [sol for sol in valid_solutions if sol[2] == min_avg]
        print(f"\nFound {len(valid_solutions)} valid SHCs in total.")
        print(f"The best solutions (with average reads = {min_avg}) are:\n")
        for idx, (solution, mapping_dict, avg_reads) in enumerate(best_solutions, start=1):
            print(f"Solution #{idx}: (Average reads: {avg_reads})")
            for i, (R, D) in enumerate(solution, start=1):
                # Convert each threshold set to a sorted list for readability.
                R_str = [sorted(list(r)) for r in R]
                print(f"  Bit {i}:")
                print(f"    Read thresholds: {R_str}")
                print(f"    Decoding rule (lookup table): {D}")
            print("  Mapping from cell levels v to output bits:")
            for v in sorted(mapping_dict):
                print(f"    v = {v}  ->  output = {mapping_dict[v]}")
            print("-" * 40)
    else:
        print("No valid SHC was found in the search space.")


if __name__ == '__main__':
    main()
