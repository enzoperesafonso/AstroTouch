#!/usr/bin/env python3
"""Run the AstroTouch command line from a clone (no install needed). Installed: `astrotouch`."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from astrotouch.core import main  # noqa: E402

if __name__ == "__main__":
    main()
