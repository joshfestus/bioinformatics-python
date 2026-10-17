# How good does a model look when the labels are random? Twice: once with
# gene selection inside cross validation, once with it done beforehand.
import numpy as np
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from expression import load

expression, samples = load()
groups = samples["cell_line"]
labels = samples["treatment"].str.contains("Dexamethasone").to_numpy()
random = np.random.default_rng(1)
cv = LeaveOneGroupOut()

honest, leaky = [], []
for _ in range(20):
    shuffled = random.permutation(labels)  # no real signal left

    inside = make_pipeline(
        StandardScaler(),
        SelectKBest(f_classif, k=50),
        LogisticRegression(max_iter=1000),
    )
    honest.append(
        cross_val_score(
            inside, expression, shuffled, groups=groups, cv=cv
        ).mean()
    )

    # The mistake: choose genes using every sample, test samples included
    chosen = SelectKBest(f_classif, k=50).fit(expression, shuffled)
    outside = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=1000)
    )
    leaky.append(
        cross_val_score(
            outside,
            chosen.transform(expression),
            shuffled,
            groups=groups,
            cv=cv,
        ).mean()
    )

print(f"random labels, genes chosen inside:  {np.mean(honest):.0%} right")
print(f"random labels, genes chosen outside: {np.mean(leaky):.0%} right")
