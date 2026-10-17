"""Name the biology behind a list of changed genes.

A list of 945 gene symbols is not an answer. Enrichment asks, for each
set of genes that share a known role, whether more of that set changed
than chance would give. The sets here are Gene Ontology biological
process terms, and the test is the hypergeometric one: drawing 492 genes
from 14,053, how surprising is it to draw 7 of the 34 that respond to
glucocorticoids?

uv run enrich.py            writes enrichment.tsv
"""

import csv
import gzip

import pandas as pd
from scipy.stats import hypergeom

PROCESS = "biological_process"
SMALLEST, LARGEST = 10, 500  # terms outside this range are not testable
MIN_FOUND = 3  # a term needs at least this many changed genes to report


def ontology(path="go-basic.obo"):
    """{term: (name, namespace)} and {term: parents}, from the OBO."""
    names, parents, current = {}, {}, None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                current = {"parents": set()}
            elif line == "" and current is not None:
                if current.get("id") and not current.get("obsolete"):
                    names[current["id"]] = (current.get("name", ""),
                                            current.get("namespace", ""))
                    parents[current["id"]] = current["parents"]
                current = None
            elif current is not None and ": " in line:
                key, value = line.split(": ", 1)
                if key in ("id", "name", "namespace"):
                    current[key] = value
                elif key == "is_a":
                    current["parents"].add(value.split(" ! ")[0])
                elif key == "relationship":
                    kind, parent = value.split()[:2]
                    if kind == "part_of":
                        current["parents"].add(parent)
                elif key == "is_obsolete":
                    current["obsolete"] = True
    return names, parents


def annotations(path="goa_human.gaf.gz"):
    """{gene symbol: terms it is directly annotated to}."""
    direct = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            column = line.split("\t")
            if "NOT" in column[3]:  # an annotation saying it does not
                continue
            direct.setdefault(column[2], set()).add(column[4])
    return direct


def with_ancestors(direct, names, parents):
    """Add every parent term, so a specific role implies the general one.

    A gene annotated to "cellular response to glucocorticoid stimulus" is
    also responding to a hormone and to a chemical, and a test that does
    not know this finds nothing at the levels people think in.
    """
    reachable, full = {}, {}

    def above(term):
        if term not in reachable:
            found, stack = set(), [term]
            while stack:
                for parent in parents.get(stack.pop(), ()):
                    if parent not in found:
                        found.add(parent)
                        stack.append(parent)
            reachable[term] = found
        return reachable[term]

    for symbol, terms in direct.items():
        every = set()
        for term in terms:
            if term in names:
                every |= {term} | above(term)
        full[symbol] = {term for term in every
                        if names[term][1] == PROCESS}
    return full


def benjamini_hochberg(pvalues):
    """The step-up correction of Chapter 19, written out."""
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted, smallest = [1.0] * len(pvalues), 1.0
    for rank, index in enumerate(reversed(order), 1):
        position = len(pvalues) - rank + 1
        smallest = min(smallest, pvalues[index] * len(pvalues) / position)
        adjusted[index] = smallest
    return adjusted


def test(sets, genes, background, names, direction):
    rows = []
    for term, members in sets.items():
        found = members & genes
        if len(found) < MIN_FOUND:
            continue
        rows.append({
            "term": term,
            "name": names[term][0],
            "direction": direction,
            "found": len(found),
            "in_set": len(members),
            "expected": round(
                len(members) * len(genes) / len(background), 1),
            "p": hypergeom.sf(len(found) - 1, len(background),
                              len(members), len(genes)),
            "genes": " ".join(sorted(found)[:8]),
        })
    for row, q in zip(rows, benjamini_hochberg([r["p"] for r in rows])):
        row["q"] = q
    return rows


def main():
    names, parents = ontology()
    full = with_ancestors(annotations(), names, parents)
    results = pd.read_csv("results.tsv", sep="\t").dropna(
        subset=["padj", "symbol"])

    background = {s for s in results["symbol"] if s in full}
    changed = results[(results["padj"] < 0.05)
                      & (results["log2FoldChange"].abs() >= 1)]
    sets = {}
    for symbol in background:
        for term in full[symbol]:
            sets.setdefault(term, set()).add(symbol)
    sets = {term: members for term, members in sets.items()
            if SMALLEST <= len(members) <= LARGEST}

    rows = []
    for direction, wanted in (("up", changed["log2FoldChange"] > 0),
                              ("down", changed["log2FoldChange"] < 0)):
        genes = {s for s in changed[wanted]["symbol"] if s in full}
        rows += test(sets, genes, background, names, direction)
        print(f"{len(genes)} genes {direction}")
    rows.sort(key=lambda row: (row["direction"], row["p"]))
    with open("enrichment.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(background):,} genes carry an annotation, and "
          f"{len(sets):,} terms are testable")
    for direction in ("up", "down"):
        found = [row for row in rows if row["direction"] == direction
                 and row["q"] < 0.05]
        print(f"\n{len(found)} terms enriched among the genes that went "
              f"{direction}:")
        for row in found[:8]:
            print(f"  q={row['q']:.0e} {row['found']:3} of "
                  f"{row['in_set']:3}  {row['name'][:42]}")


if __name__ == "__main__":
    main()
