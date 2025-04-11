import itertools

def decision(v, thresholds):
    b = 1
    for r in sorted(thresholds):
        if r <= v:
            b ^= 1  # Flip bit if threshold is passed
        else:
            break
    return b

def decode_bits(cell_levels, shc):
    bits = []
    for R, D in shc:
        # Apply decision rule for each cell level using its threshold set
        decision_bits = [decision(cell_levels[i], R[i]) for i in range(len(R))]
        # Convert binary vector to index
        index = sum(bit << i for i, bit in enumerate(decision_bits))
        bits.append(D[index])
    return tuple(bits)

def is_valid_shc(shc, n, l, k):
    all_outputs = set()
    # Check output for all possible cell states
    for cell_levels in itertools.product(range(l), repeat=n):
        output = decode_bits(cell_levels, shc)
        all_outputs.add(output)
    # Must produce all 2^k possible outputs
    return len(all_outputs) == 2 ** k

def average_reads_per_bit(shc):
    total = 0
    for R, _ in shc:
        # Each bit reads the cell with the largest number of thresholds
        max_thresholds = max(len(cell_thresholds) for cell_thresholds in R)
        total += max_thresholds
    return total / len(shc)

if __name__ == "__main__":
    example_config = [
        ((frozenset({1, 3, 6, 10, 14, 18, 23, 27}),), (0, 1)),
        ((frozenset({2, 5, 9, 13, 17, 21, 25, 29}),), (0, 1)),
        ((frozenset({4, 8, 12, 16, 20, 24, 28}),), (0, 1)),
        ((frozenset({7, 15, 23}),), (0, 1)),
        ((frozenset({16}),), (0, 1))
    ]
    n = 1
    l = 32
    k = 5

    valid = is_valid_shc(example_config, n, l, k)
    print(f"SHC is valid: {valid}")

    avg_reads = average_reads_per_bit(example_config)
    print(f"Average reads per bit: {avg_reads:.2f}")
