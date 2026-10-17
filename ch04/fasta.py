"""A FASTA reader, written from scratch."""


def read_fasta(path: str) -> dict[str, str]:
    """Map each record's header line to its sequence, in capitals."""
    sequences = {}
    header = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                sequences[header] = []
            elif header is None:
                raise ValueError(f"{path}: sequence before any header")
            else:
                sequences[header].append(line.upper())
    return {header: "".join(lines) for header, lines in sequences.items()}


def gene_name(header: str) -> str:
    """The value of [gene=...] in an NCBI header, or its first word."""
    if "[gene=" not in header:
        return header.split()[0]
    return header.split("[gene=")[1].split("]")[0]
