# Group the samples by how alike their expression is, with no labels used.
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from expression import load

expression, samples = load()
variable = expression.var().sort_values(ascending=False).index[:500]

distances = pdist(expression[variable], metric="correlation")
tree = linkage(distances, method="average")
samples["cluster"] = fcluster(tree, t=2, criterion="maxclust")
print(
    samples.sort_values(["cluster", "cell_line"])[
        ["cell_line", "label", "cluster"]
    ].to_string()
)
