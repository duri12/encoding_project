import itertools
from tqdm import tqdm
import math

# Parameters for the SHC: n cells, l levels, k bits.
n = 2  # number of cells per codeword
l = 4 # levels
k = 4  # number of bits to store

MAX_THRESHOLDS_NUM = l #NOTE: set the maximal amount of threshold for a bit in a cell.
# possible inner-groups for R (single cell): only singletons and empty set.
# candidate_thresholds = [frozenset()]
# for r in range(1, l):
#     candidate_thresholds.append(frozenset({r}))

s = list(range(1, l))

candidate_thresholds = [frozenset(subset) for r in range(MAX_THRESHOLDS_NUM+1) for subset in itertools.combinations(s, r)]
print(candidate_thresholds)

# get all values for R.
candidate_Rs = list(itertools.product(candidate_thresholds, repeat=n))

# Generate all possible decoding functions D: {0,1}^n -> {0,1}.
# Exclude constant functions.
def generate_candidate_D(n):
    size = 2 ** n
    candidates = []
    for bits in itertools.product([0, 1], repeat=size):
        if not (all(b == bits[0] for b in bits)):# is constant function? i.e. f(i) == f(0) for every i
            candidates.append(bits)
    return candidates

#on which bits does d not depend on
def d_not_depend_on(candidate_d):
    num_of_bits = n
    f = candidate_d
    dependent_bits = set()
    for bit in range(num_of_bits):
        for i in range(len(f)): # all indexes
            j = i ^ (1 << bit)  # Flip the `bit`-th position
            if j > i and f[i] != f[j]:
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


def check_prop1(candidate_r,candidate_d):
    counters =[0,0]
    THRESHOLD = pow(2,k-1)
    # print(candidate_d, candidate_r)
    for v in itertools.product(range(l), repeat=n):
        out = decode_bit((candidate_r,candidate_d), v)
        # print(v,out)
        counters[out] += 1
        if counters[0] >= THRESHOLD and counters[1] >= THRESHOLD:
            return True
    return False

candidate_Ds = generate_candidate_D(n)


# Candidate pairs are combinations of R and D
candidate_pairs_original = list(itertools.product(candidate_Rs, candidate_Ds)) #TODO: check redundant groups and switch them to empty
candidate_pairs = []
for d in candidate_Ds:
    not_dependent_bits = d_not_depend_on(d)
    for r in candidate_Rs:
        new_r = list(r)
        #optimization NEW 1: remove unnecessary R
        for bit in not_dependent_bits:
            new_r[bit] = frozenset()
        
        #optimization PAPER 1: check surjective of R & D.
        if not check_prop1(new_r, d):
            continue
        if (new_r,d) not in candidate_pairs:
            candidate_pairs.append((new_r, d))

print("After applying acceleration ",len(candidate_pairs))
print("Before applying acceleration ",len(candidate_pairs_original))

# Decision function (Algorithm 1)
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