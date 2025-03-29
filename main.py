import itertools
from tqdm import tqdm
import math

# Parameters for the SHC: n cells, l levels, k bits.
n = 1  # number of cells per codeword
l = 16  # levels
k = 4  # number of bits to store

MAX_THRESHOLDS_NUM = 3  # NOTE: set the maximal amount of threshold for a bit in a cell
# possible inner-groups for R (single cell): only singletons and empty set.
# candidate_thresholds = [frozenset()]
# for r in range(1, l):
#     candidate_thresholds.append(frozenset({r}))

s = list(range(1, l))

candidate_thresholds = [frozenset(subset) for r in range(MAX_THRESHOLDS_NUM + 1) for subset in
                        itertools.combinations(s, r)]
# get all values for R.
candidate_Rs = list(itertools.product(candidate_thresholds, repeat=n))


# Generate all possible decoding functions D: {0,1}^n -> {0,1}.
# Exclude constant functions.
def generate_candidate_D(n):
    size = 2 ** n
    # Generate all possible non-constant functions
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if not all(b == bits[0] for b in bits):  # is constant function?
            candidates.append(bits)
    return candidates

# on which bits does d not depend on
def d_not_depend_on(candidate_d):
    num_of_bits = n
    f = candidate_d
    dependent_bits = set()
    for bit in range(num_of_bits):
        for i in range(len(f)):  # all indexes
            j = i ^ (1 << bit)  # Flip the `bit`-th position
            if j < len(f) and f[i] != f[j]:  # Check j is within bounds
                dependent_bits.add(bit)
    not_dependent_bits = set(range(num_of_bits)) - dependent_bits
    return not_dependent_bits

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


def check_prop1(candidate_r, candidate_d):
    counters = [0, 0]
    THRESHOLD = pow(2, k - 1)
    # print(candidate_d, candidate_r)
    for v in itertools.product(range(l), repeat=n):
        out = decode_bit((candidate_r, candidate_d), v)
        # print(v,out)
        counters[out] += 1
        if counters[0] >= THRESHOLD and counters[1] >= THRESHOLD:
            return True
    return False

candidate_Ds = generate_candidate_D(n)

# Candidate pairs are combinations of R and D
candidate_pairs_original = list(
    itertools.product(candidate_Rs, candidate_Ds))  # TODO: check redundant groups and switch them to empty


def get_canonical_form(R, D, n):
    """For n=2, return the lex smallest permutation of R and its corresponding D."""
    if n != 2:
        return (R, D)  # Only implemented for n=2

    # Generate swapped R and adjusted D
    R_swapped = (R[1], R[0])
    # For swapped cells, D's truth table is [D[0], D[2], D[1], D[3]]
    D_swapped = (D[0], D[2], D[1], D[3])

    # Compare original and swapped forms lexicographically
    original = (R, D)
    swapped = (R_swapped, D_swapped)
    return min(original, swapped)


# Modify the candidate_pairs generation loop:
seen_canonical = set()
candidate_pairs = []
for d in candidate_Ds:
    not_dependent_bits = d_not_depend_on(d)
    for r in candidate_Rs:
        new_r = list(r)
        # Apply optimizations
        for bit in not_dependent_bits:
            new_r[bit] = frozenset()
        if not check_prop1(new_r, d):
            continue

        # Convert to tuple for hashability
        new_r_tuple = tuple(new_r)

        # Get canonical form (for n=2)
        canonical_r, canonical_d = get_canonical_form(new_r_tuple, d, n)

        # Check if canonical form already seen
        if (canonical_r, canonical_d) not in seen_canonical:
            seen_canonical.add((canonical_r, canonical_d))
            candidate_pairs.append((new_r_tuple, d))
print("Before acceleration:", len(candidate_pairs_original))
print("After applying acceleration and canonicalization:", len(candidate_pairs))


# Decision function (Algorithm 1)
# Calculate average number of reads per bit
def average_reads(candidate_solution):
    total_reads = 0
    for R, _ in candidate_solution:
        total_reads = max(total_reads, max(len(r) for r in R))
    return total_reads


def main():
    print(f"Searching for SHC({n}, {l}, {k}) with thresholds: {[set(t) for t in candidate_thresholds]}")
    print(f"Possible R per cell: {len(candidate_thresholds)}, Total R combinations: {len(candidate_Rs)}")
    print(f"Possible D functions: {len(candidate_Ds)}, Candidate pairs per bit: {len(candidate_pairs)}")

    # TODO: make the check of all valid pairs and remove is_surjective. use memoization
    valid_solutions = []
    total_combinations = math.comb(len(candidate_pairs), k)

    min_reads = 100000000  # l
    # Iterate over all possible combinations of k pairs with a progress bar
    for candidate_solution in tqdm(itertools.combinations(candidate_pairs, k), total=total_combinations,
                                   desc="Searching"):
        avg_reads = average_reads(candidate_solution)
        if avg_reads > min_reads:
            continue
        mapping_dict = mapping(candidate_solution)
        if is_surjective(mapping_dict):  # check that capacity is 1 and change the check
            if min_reads > avg_reads:
                min_reads = avg_reads
                valid_solutions = [(candidate_solution, mapping_dict, avg_reads)]
            elif min_reads == avg_reads:
                valid_solutions.append((candidate_solution, mapping_dict, avg_reads))

    if valid_solutions:
        min_avg = min(sol[2] for sol in valid_solutions)
        best_solutions = [sol for sol in valid_solutions if sol[2] == min_avg]
        for sol in best_solutions:
            print("\nSolution:")
            for i, (R, D) in enumerate(sol[0], 1):
                print(f"Bit {i}: R={[set(r) for r in R]}, D={D}")
            print(f"Average reads per bit: {sol[2]}")
        print(f"Found {len(valid_solutions)} valid SHCs. Best avg reads: {min_avg}")

    else:
        print("No valid SHC found.")


if __name__ == '__main__':
    main()