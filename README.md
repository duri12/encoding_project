# SHC Code Search

This project searches for valid **SHC(n, l, k)** codes—mappings from multi-dimensional inputs to binary outputs. A code is valid if it's **surjective** (covers all 2^k outputs) and uses a **low average read cost** (how many thresholds per bit).

---

## Key Scripts

### `search_BF_with_opt.py` — Optimized

- Batches threshold generation (faster, less memory)
- Shuffles candidates to hit good ones early
- Skips duplicate D vectors via symmetry
- Prunes early if current read cost isn't better than best


---

### `search_BF_without_opt.py` — Unoptimized

- Generates only **balanced + minimal** threshold sets
- Uses two simple D mappings
- Processes in **parallel chunks** (multiprocessing)
- No pruning — checks everything


---

### `verify_valid_code.py`
(inner tool)
Validates any SHC config:
- Surjectivity check
- Average read cost
---


## How to Use

To run the exhaustive search script for finding SHC codes, execute `search_BF_with_opt.py` from the command line with the following required arguments:

```bash
python search_BF_with_opt.py --n <number_of_cells> --l <alphabet_size> --k <number_of_output_bits>
```

where:
- `--n` = number of cells (dimensions),
- `--l` = number of levels per cell (alphabet size),
- `--k` = number of bits to encode.

### Example

```bash
python search_BF_with_opt.py --n 1 --l 8 --k 3
```

This command will search for an SHC(1,8,3) code using the optimized search algorithm.


## Requirements

- **Python 3.6+**
- **tqdm** for progress bars  
  Install with:
  ```bash
  pip install tqdm
  ```
