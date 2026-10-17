"""Solutions to the exercises of Chapter 22.

Run from this folder: uv run ch22.py
Exercises 22.1, 22.3 and 22.4 work on a throwaway copy of ../ch22, so
nothing here changes the chapter's own files. The Snakemake runs need
the mapped reads, and are skipped when they are missing.
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CH22 = Path("../ch22").resolve()
DATA = Path("../data").resolve()


def tool(name):
    """A command line tool, looking beside this interpreter first."""
    beside = Path(sys.executable).parent / name
    if beside.exists():
        return str(beside)
    return shutil.which(name)


def working_copy(work):
    """A copy of the chapter folder we can break safely.

    The programs there read their data through relative paths such as
    ../ch15/rel606.gff3, so the copy gets symbolic links to the real
    sibling folders and every path resolves as it does in the book.
    """
    copy = Path(work) / "ch22"
    shutil.copytree(CH22, copy,
                    ignore=shutil.ignore_patterns("__pycache__",
                                                  ".snakemake", "results",
                                                  "logs", "*.ipynb"))
    for sibling in ("ch15", "ch16", "data"):
        target = CH22.parent / sibling
        if target.exists():
            (Path(work) / sibling).symlink_to(target)
    return copy


def pytest_in(folder, extra=()):
    done = subprocess.run([sys.executable, "-m", "pytest", "-q", *extra],
                          cwd=folder, capture_output=True, text=True)
    return done.stdout + done.stderr


def exercise_22_1():
    print("22.1  Breaking the strand arithmetic on purpose.")
    with tempfile.TemporaryDirectory() as work:
        copy = working_copy(work)
        source = (copy / "variants.py").read_text()
        broken = source.replace("end - position", "end - position - 1", 1)
        if broken == source:
            print("      Could not find the line to change.")
            return
        (copy / "variants.py").write_text(broken)
        report = pytest_in(copy)
        failed = [line for line in report.splitlines()
                  if line.startswith("FAILED") or " failed" in line]
        print("       with the bug in place:")
        for line in failed[:4]:
            print(f"         {line.strip()[:66]}")
        detail = [line for line in report.splitlines()
                  if "assert" in line and "==" in line]
        for line in detail[:2]:
            print(f"         {line.strip()[:66]}")
        (copy / "variants.py").write_text(source)
        report = pytest_in(copy)
        passed = [line for line in report.splitlines()
                  if "passed" in line]
        print("       with it undone:")
        for line in passed[:2]:
            print(f"         {line.strip()[:66]}")
    print("      The minus strand test is the one that fails, and it")
    print("      fails loudly: pytest prints the call, the expected")
    print("      answer and the one it got. The three plus strand")
    print("      tests pass throughout, which is the point of having")
    print("      more than one: the failure names the direction the")
    print("      bug is in.")
    print("      This is the test earning its keep. An off by one in")
    print("      strand arithmetic produces a plausible amino acid")
    print("      change at a plausible position, and nothing but a")
    print("      test with a known answer will catch it.")


def exercise_22_2():
    print("\n22.2  Raising min_depth from 10 to 20.")
    if not (CH22 / "results").exists():
        print("      The pipeline has not been run here, so there is")
        print("      nothing for Snakemake to compare against. Run")
        print("      `snakemake --cores 2 all` in ../ch22 first.")
    print("      Snakemake tracks the parameters each output was made")
    print("      with, so changing min_depth in config.yaml makes the")
    print("      annotate step out of date and nothing else. Mapping")
    print("      and calling do not rerun, because neither of them")
    print("      reads that value.")
    sys.path.insert(0, str(CH22))
    try:
        from variants import annotate
    except ImportError as problem:
        print(f"      Could not import variants.py ({problem}).")
        return
    print(f"\n       {'min_depth':>10} {'variants':>9}")
    for depth in (5, 10, 20, 30):
        try:
            # annotate yields its rows, so collect them to count
            rows = list(annotate(
                str(CH22 / ".." / "ch16" / "haploid.vcf.gz"),
                str(CH22 / ".." / "ch15" / "rel606.gff3"),
                str(DATA / "rel606.fasta"),
                min_qual=30, min_depth=depth))
        except (OSError, ValueError) as problem:
            print(f"      Could not run the annotation ({problem}).")
            return
        print(f"       {depth:>10} {len(rows):>9}")
    print("      Fewer variants survive as the threshold rises, and")
    print("      the ones that drop out are the shallow calls, where")
    print("      a handful of reads carried the difference. Whether")
    print("      that is a gain depends on what the list is for: a")
    print("      stricter threshold is a list you can trust gene by")
    print("      gene, a looser one is a list that misses less.")
    print("      The number belongs in the report either way, which")
    print("      is why it lives in config.yaml and not in the code.")


def exercise_22_3():
    print("\n22.3  A flagstat rule, added to the Snakefile:\n")
    rule = '''rule flagstat:
    input:
        "results/reads.bam",
    output:
        "results/flagstat.txt",
    log:
        "logs/flagstat.log",
    shell:
        "samtools flagstat {input} > {output} 2> {log}"'''
    for line in rule.splitlines():
        print(f"       {line}")
    print("\n      and results/flagstat.txt added to rule all's input.")
    with tempfile.TemporaryDirectory() as work:
        copy = working_copy(work)
        text = (copy / "Snakefile").read_text()
        text = text.replace('rule all:\n    input:\n',
                            'rule all:\n    input:\n'
                            '        "results/flagstat.txt",\n')
        (copy / "Snakefile").write_text(text + "\n\n" + rule + "\n")
        snakemake = tool("snakemake")
        output = []
        if snakemake:
            done = subprocess.run(
                [snakemake, "-n", "--quiet", "host", "rules"], cwd=copy,
                capture_output=True, text=True)
            output = (done.stdout + done.stderr).splitlines()
        wanted = [line for line in output
                  if "flagstat" in line or "total" in line]
        if wanted:
            print("\n       a dry run now plans it:")
            for line in wanted[:4]:
                print(f"         {line.strip()[:64]}")
        elif snakemake:
            first = [line for line in output if line.strip()][:2]
            print("\n       snakemake plans:")
            for line in first:
                print(f"         {line.strip()[:64]}")
        else:
            print("\n       (snakemake is not on the PATH here)")
    print("      Adding it to rule all is the part people forget.")
    print("      Snakemake makes what it is asked for and works")
    print("      backwards, so a rule nothing depends on never runs,")
    print("      and a new step that quietly does nothing is worse")
    print("      than no step at all.")


def exercise_22_4():
    print("\n22.4  Two more tests, written against known answers:\n")
    tests = '''def test_pykf_is_found_on_the_plus_strand():
    features = read_coding_sequences("../ch15/rel606.gff3")
    pykf = [f for f in features if f["name"] == "pykF"]
    assert len(pykf) == 1
    assert pykf[0]["strand"] == "+"
    assert pykf[0]["start"] == 1_732_965
    assert pykf[0]["end"] == 1_734_377


def test_a_variant_between_genes_is_reported_as_such():
    rows = annotate("../ch16/haploid.vcf.gz", "../ch15/rel606.gff3",
                    "../data/rel606.fasta", min_qual=30, min_depth=10)
    between = [r for r in rows if "between" in r["effect"]]
    assert len(between) == 5
    assert all(r["position"] for r in between)'''
    for line in tests.splitlines():
        print(f"       {line}")
    with tempfile.TemporaryDirectory() as work:
        copy = working_copy(work)
        path = copy / "test_variants.py"
        text = path.read_text()
        old = "from variants import amino_acid_change, annotate"
        new = ("from variants import (amino_acid_change, annotate,\n"
               "                      read_coding_sequences)")
        text = text.replace(old, new)
        path.write_text(text + "\n\n" + tests + "\n")
        report = pytest_in(copy)
        summary = [line for line in report.splitlines()
                   if "passed" in line or "failed" in line]
        print("\n       running them:")
        for line in summary[:3]:
            print(f"         {line.strip()[:66]}")
    print("      Both check facts settled elsewhere in the book:")
    print("      Chapter 16 found pykF on the plus strand at 1,732,965")
    print("      and five variants between genes, and Chapter 17")
    print("      confirmed that split independently. A test is only")
    print("      worth writing when you know the answer without")
    print("      running the code, which is why those two chapters")
    print("      are the right place to take it from.")


def exercise_22_5():
    print("\n22.5  Chapter 25 is this exercise, done at full size.")
    print("      Part Five's RNA-seq work became projects/airway: a")
    print("      download step, a sample sheet, a differential")
    print("      expression run, an enrichment, a comparison with the")
    print("      paper, and a report, wired together by a Snakefile")
    print("      with its thresholds in config.yaml.")
    print("      The part worth testing is the enrichment, because")
    print("      it is the only step whose arithmetic is ours rather")
    print("      than a library's:\n")
    test = '''def test_benjamini_hochberg_matches_scipy():
    from scipy.stats import false_discovery_control
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.6, 0.9]
    mine = benjamini_hochberg(p)
    theirs = false_discovery_control(p)
    assert max(abs(a - b) for a, b in zip(mine, theirs)) < 1e-12


def test_ancestors_reach_the_top_of_the_ontology():
    names, parents = ontology("go-basic.obo")
    full = with_ancestors({"TEST": {"GO:0071385"}}, names, parents)
    # cellular response to glucocorticoid stimulus is a kind of
    # response to a hormone, and the propagation must say so
    assert "GO:0032870" in full["TEST"]'''
    for line in test.splitlines():
        print(f"       {line}")
    print("\n      The first test is the one that matters most. The")
    print("      Benjamini-Hochberg step-up is four lines and easy to")
    print("      get subtly wrong, and SciPy's version is an")
    print("      independent implementation of the same procedure, so")
    print("      agreeing with it to twelve decimal places is real")
    print("      evidence rather than a restatement of the code.")


if __name__ == "__main__":
    exercise_22_1()
    exercise_22_2()
    exercise_22_3()
    exercise_22_4()
    exercise_22_5()
