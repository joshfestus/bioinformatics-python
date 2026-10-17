"""Check that this computer is ready for the rest of the book."""

import shutil
import sys
from importlib.metadata import PackageNotFoundError, version

PACKAGES = ["numpy", "pandas", "matplotlib", "biopython"]

# Each program, with the chapter that first needs it
PROGRAMS = {"blastn": 11, "mafft": 12, "fastqc": 14, "minimap2": 15}
PROGRAMS |= {"samtools": 15, "bcftools": 16}

found = sys.version_info[:2]
status = "ok" if found >= (3, 14) else "too old, install 3.14"
print(f"Python {sys.version.split()[0]:12} {status}")

print("\nPython packages")
for name in PACKAGES:
    try:
        print(f"  {name:12} {version(name):12} ok")
    except PackageNotFoundError:
        print(f"  {name:12} {'missing':12} run uv sync")

print("\nCommand line programs")
for name, chapter in PROGRAMS.items():
    if shutil.which(name):
        print(f"  {name:12} {'found':12} ok")
    else:
        print(
            f"  {name:12} {'not yet':12} installed in Chapter {chapter}"
        )
