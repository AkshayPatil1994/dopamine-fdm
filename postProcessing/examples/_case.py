"""Shared by the examples: the case directory is argv[1] when given, otherwise the shipped example case."""
import sys
from pathlib import Path


def case_dir(default):
    return Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "examples" / default
