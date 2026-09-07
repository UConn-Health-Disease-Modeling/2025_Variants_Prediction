#!/usr/bin/env python3
"""Run the binary model comparison for the growth response only."""

from pathlib import Path

from modeling_comparison import main


if __name__ == "__main__":
    main(
        fixed_outcome="growth",
        default_output_dir=Path("code/modeling_results/growth"),
    )
