"""Global alignment by the Needleman-Wunsch method, from scratch."""

MATCH, MISMATCH, GAP = 1, -1, -2


def score_pair(x: str, y: str) -> int:
    return MATCH if x == y else MISMATCH


def fill(a: str, b: str) -> list[list[int]]:
    """The scoring grid: grid[i][j] is the best score for a[:i], b[:j]."""
    grid = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    for i in range(1, len(a) + 1):
        grid[i][0] = i * GAP
    for j in range(1, len(b) + 1):
        grid[0][j] = j * GAP
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            grid[i][j] = max(
                grid[i - 1][j - 1] + score_pair(a[i - 1], b[j - 1]),
                grid[i - 1][j] + GAP,  # a letter of a against a gap
                grid[i][j - 1] + GAP,  # a letter of b against a gap
            )
    return grid


def traceback(a: str, b: str, grid: list[list[int]]) -> tuple[str, str]:
    """Walk back from the bottom right corner to recover the alignment."""
    top, bottom = [], []
    i, j = len(a), len(b)
    while i > 0 or j > 0:
        step = score_pair(a[i - 1], b[j - 1]) if i and j else None
        if step is not None and grid[i][j] == grid[i - 1][j - 1] + step:
            top.append(a[i - 1])
            bottom.append(b[j - 1])
            i, j = i - 1, j - 1
        elif i > 0 and grid[i][j] == grid[i - 1][j] + GAP:
            top.append(a[i - 1])
            bottom.append("-")
            i -= 1
        else:
            top.append("-")
            bottom.append(b[j - 1])
            j -= 1
    return "".join(reversed(top)), "".join(reversed(bottom))


def mark(x: str, y: str) -> str:
    """| for a match, . for a mismatch, a space for a gap."""
    if "-" in (x, y):
        return " "
    return "|" if x == y else "."


def show(top: str, bottom: str) -> None:
    marks = "".join(mark(x, y) for x, y in zip(top, bottom))
    print(top, marks, bottom, sep="\n")


if __name__ == "__main__":
    a, b = "ATGAGCA", "AGACCA"
    grid = fill(a, b)
    for row in grid:
        print(" ".join(f"{value:3}" for value in row))
    print("Score:", grid[-1][-1])
    show(*traceback(a, b, grid))
