# Counts per million of four genes dexamethasone is known to switch on.
from expression import load

MARKERS = {2289: "FKBP5", 7704: "ZBTB16", 28999: "KLF15", 5187: "PER1"}

expression, samples = load()
cpm = 2 ** expression[list(MARKERS)] - 1  # back from log2(cpm + 1)
cpm.columns = list(MARKERS.values())
table = samples[["cell_line", "label"]].join(cpm.round(0).astype(int))
print(table.sort_values(["cell_line", "label"]).to_string())
