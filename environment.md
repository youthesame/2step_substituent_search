# Execution Environment

## 1. PC / VM Resources
| resource | value | how-to-reproduce |
|----------|-------|------------------|
| CPU      | Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz (8 cores) | `lscpu \| grep -E 'Model name\|^CPU'` |
| RAM      | 31Gi total, 27Gi available | `free -h` |
| GPU (CUDA) | NVIDIA RTX 4000 Ada Generation (20475MiB VRAM, CUDA 12.4) | `nvidia-smi` |
| OS / Kernel | Linux 6.6.87.2-microsoft-standard-WSL2 (Ubuntu 24.04) | `uname -a` |

## 2. Python Environment
| item | value | how-to-reproduce |
|------|-------|------------------|
| Python | Python 3.11.13 | `uv run python --version` |
| Package manager | uv (see `uv.lock`) |  |
| Virtual-env path | `.venv/` |  |
| Locked dependencies | see [`uv.lock`](./uv.lock) |  |

> Tip – if you ever switch to conda, replace this table with `conda env export --no-builds`.

## 3. Directory Layout (top level)
```text
workspace/
├── CLAUDE.md
├── README.md
├── _docs/         # documentation logs
│   └── templates/
├── environment.md # (this file)
├── out/           # computation output files
│   └── DNTT_fullopt.out
├── pyproject.toml # project metadata / build settings
├── src/           # source code (empty)
└── uv.lock        # python deps lock-file
```

_Generate with:_ `tree -L 2 -I ".git|.venv"`