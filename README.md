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

Validates any SHC config:
- Surjectivity check
- Average read cost
---

## Requirements

- **Python 3.6+**
- **tqdm** for progress bars  
  Install with:
  ```bash
  pip install tqdm
  ```