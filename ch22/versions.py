# Record exactly which versions produced a result, for the methods section
# of a paper and for anyone who reruns the analysis.
import platform
import subprocess
from importlib.metadata import version

PACKAGES = ["biopython", "pysam", "pandas", "numpy", "snakemake"]
TOOLS = {
    "minimap2": ["--version"],
    "samtools": ["--version"],
    "bcftools": ["--version"],
}

lines = [f"python {platform.python_version()}"]
lines += [f"{name} {version(name)}" for name in PACKAGES]
for tool, arguments in TOOLS.items():
    output = subprocess.run(
        [tool, *arguments],
        capture_output=True,
        text=True,
        errors="replace",  # one tool prints a byte that is not UTF-8
        check=True,
    ).stdout
    first = output.splitlines()[0].split()
    lines.append(f"{tool} {first[-1]}")

with open("versions.txt", "w") as handle:
    handle.write("\n".join(lines) + "\n")
print("\n".join(lines))
