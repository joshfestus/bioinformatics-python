# Principal components of the sixteen airway samples.
from sklearn.decomposition import PCA

from expression import load

expression, samples = load()
variable = expression.var().sort_values(ascending=False).index[:500]
centered = expression[variable] - expression[variable].mean()

pca = PCA(n_components=3)
scores = pca.fit_transform(centered)
shares = ", ".join(f"{s:.0%}" for s in pca.explained_variance_ratio_)
print(f"share of variation on PC1, PC2, PC3: {shares}")

samples[["PC1", "PC2", "PC3"]] = scores.round(1)
print(samples[["cell_line", "label", "PC1", "PC2", "PC3"]].to_string())
samples.to_csv("pca_scores.tsv", sep="\t")
