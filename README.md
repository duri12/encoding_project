# SHC Search and Verification

This project implements and searches for **SHC (SHC(n, l, k)) codes**, a family of codes that map multi-dimensional inputs (or _cell levels_) into a binary output space. In an SHC code, the key parameters are:

- **n:** The number of dimensions or cells.
- **l:** The alphabet size for each dimension (i.e. the number of distinct values each cell can take).
- **k:** The number of output bits.

An SHC code must be **surjective**—that is, when applied to every combination of inputs (cell levels), it produces all possible values in the binary output space (of size 2^k). The search and verification processes in this repository help discover valid SHC configurations while also minimizing the average read cost per bit.

## Project Overview

At a high level, the project performs the following:

1. **Candidate Generation:**  
   Generate candidate threshold sets and bit mapping (denoted as D) functions. The thresholds decide when a cell’s contribution toggles the bit value (using a decision function that iterates over sorted thresholds).

2. **Mapping and Validation:**  
   Each candidate configuration (a pairing of threshold sets and a D mapping) is used to construct a mapping over the input space. A configuration is considered valid if the mapping is surjective (covers every possible output). Additionally, an average “read cost” (how many thresholds each bit uses) is computed and minimized.

3. **Search Optimization:**  
   Two different search approaches are provided:
   - **With optimizations (`search_BF_with_opt.py`):** Uses batching, shuffling, and a candidate evaluation process to reduce computation.
   - **Without optimizations (`search_BF_without_opt.py`):** Explores candidate configurations using a multiprocessing approach and applies criteria such as balance and minimality of threshold sets.

4. **Verification:**  
   A verification script (`verify_valid_code.py`) takes an SHC configuration and tests whether it truly is surjective (i.e. every binary output is produced) and reports the average read cost per bit.

## File Breakdown

### 1. `search_BF_with_opt.py`

- **Purpose:**  
  Searches for valid SHC configurations using an optimized brute-force approach.
  
- **Key Functions:**  
  - `generate_candidate_thresholds(l)`: Generates candidate threshold sets for a given alphabet size.
  - `generate_candidate_D(n)`: Constructs candidate mappings (D arrays) for n dimensions.
  - `decision(v, R)`: Computes the decision bit based on the input value and thresholds.
  - `mapping()`: Builds the full mapping from cell level inputs to output bits.
  - The main loop iterates over combinations of candidate pairs (each pair is a threshold set and a D mapping) and prints the solution with the lowest average read cost.

- **Usage:**  
  Run from the command line with arguments:
  ```bash
  python search_BF_with_opt.py --n <number_of_dimensions> --l <alphabet_size> --k <num_output_bits>
  ```

### 2. `search_BF_without_opt.py`

- **Purpose:**  
  Provides an alternative search method using multiprocessing. It leverages candidate generation that enforces balanced and minimal threshold sets.
  
- **Key Functions:**  
  - `generate_candidates()`: Iterates over possible candidate sets, verifying they are balanced and minimal.
  - `process_chunk(chunk)`: Processes a chunk of candidate pairs and screens them based on the uniqueness of their bit mapping.
  - `parallel_search()`: Runs the search in parallel using multiple processors.
  
- **Usage:**  
  Execute the script from the command line:
  ```bash
  python search_BF_without_opt.py --n <number_of_dimensions> --l <alphabet_size> --k <num_output_bits>
  ```
  This script uses Python's multiprocessing module to improve search efficiency.

### 3. `verify_valid_code.py`

- **Purpose:**  
  Verifies whether a given SHC configuration is valid. It checks:
  - **Surjectivity:** Ensures that the mapping covers all 2^k possible outputs.
  - **Read Cost:** Computes the average number of thresholds (reads) per bit.
  
- **Key Functions:**  
  - `is_valid_shc(shc, n, l, k)`: Iterates over all input combinations and confirms that the outputs are surjective.
  - `average_reads_per_bit(shc)`: Calculates the average cost in terms of the maximum thresholds needed per bit.
  
- **Usage:**  
  You can run the script as follows to test an example configuration:
  ```bash
  python verify_valid_code.py
  ```
  The script prints whether the configuration is valid and the average reads per bit. You can replace the `example_config` variable with your own candidate configuration to test different SHC codes.

## How SHC Works

In more detail, an SHC configuration consists of several candidate pairs, each made up of:
- A set of **thresholds** for each cell (R), which determine when a value flip occurs (using the `decision` function).
- A bit mapping **D** that assigns a binary digit based on the decision results.

By combining multiple candidate pairs, an overall mapping is produced using a multi-dimensional product over the input space. The combined output must cover every possible binary string of length *k* exactly once (hence being surjective).

The search scripts explore various candidates by:
- Shuffling and combining candidate thresholds.
- Evaluating the effect on the complete mapping from inputs to binary outputs.
- Keeping track of and minimizing the “read” cost—the number of thresholds that need to be processed for each bit decision.

## Dependencies

- **Python 3.6 or higher**
- **tqdm:** For progress bars during searches.
- **argparse:** For command-line argument parsing.

Install required packages (if not already installed) using pip:
```bash
pip install tqdm
```

