import itertools

def decision(v, thresholds):
    """
    Compute the decision bit for a cell level 'v' using the given thresholds (Algorithm 1).
    thresholds: a set of integers representing the read thresholds.
    Returns: 0 or 1.
    """
    b = 1
    for r in sorted(thresholds):
        if r <= v:
            b ^= 1
        else:
            break
    return b

def decode_bits(cell_levels, shc):
    """
    Compute the k output bits for a given cell_levels vector using the SHC.
    cell_levels: tuple of cell levels (v1, v2, ..., vn).
    shc: list of tuples [(R1, D1), (R2, D2), ..., (Rk, Dk)].
    Returns: tuple of bits (u1, u2, ..., uk).
    """
    bits = []
    for R, D in shc:
        # Compute decision bits for each cell in R
        decision_bits = [decision(cell_levels[i], R[i]) for i in range(len(R))]
        # Convert decision_bits to an index (LSB first)
        index = sum(bit << i for i, bit in enumerate(decision_bits))
        # Get the decoded bit from D
        bits.append(D[index])
    return tuple(bits)

def is_valid_shc(shc, n, l, k):
    """
    Check if the SHC is valid (surjective).
    Returns: True if valid, False otherwise.
    """
    all_outputs = set()
    # Generate all possible cell level combinations (v1, v2, ..., vn)
    for cell_levels in itertools.product(range(l), repeat=n):
        output = decode_bits(cell_levels, shc)
        all_outputs.add(output)
    # Check if all 2^k possible outputs are covered
    return len(all_outputs) == 2 ** k

def average_reads_per_bit(shc):
    """
    Calculate the average number of reads per bit.
    """
    total = 0
    for R, _ in shc:
        # Max thresholds per cell for this bit
        max_thresholds = max(len(cell_thresholds) for cell_thresholds in R)
        total += max_thresholds
    return total / len(shc)

# Example usage:
if __name__ == "__main__":
    # Example SHC(2,3,3) from the paper (simplified)
    # Each entry is (R, D), where R is thresholds per cell, D is decoding table
    shc_example = [
        ((frozenset({0, 5, 10, 15, 20, 25, 30}),), (0, 1)),
        ((frozenset({1, 6, 11, 16, 21, 26, 29}),), (0, 1)),
        ((frozenset({2, 7, 12, 17, 22, 27, 31}),), (0, 1)),
        ((frozenset({3, 8, 13, 18, 23, 28, 30}),), (0, 1)),
        ((frozenset({4, 9, 14, 19, 24, 29, 31}),), (0, 1))
    ]

    n = 1
    l = 32
    k = 5

    # Verify the SHC
    valid = is_valid_shc(shc_example, n, l, k)
    print(f"SHC is valid: {valid}")

    # Calculate average reads
    avg_reads = average_reads_per_bit(shc_example)
    print(f"Average reads per bit: {avg_reads:.2f}")