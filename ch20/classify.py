# Learn to recognize dexamethasone treated samples from expression, then
# ask the model about samples it has never seen.
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from expression import load

expression, samples = load()
treated = samples["treatment"].str.contains("Dexamethasone").astype(int)
known = samples["label"].isin(["none", "dex"])

# Every step that learns from the data sits inside the pipeline, so each
# round of cross validation repeats it using its training samples only.
model = make_pipeline(
    StandardScaler(),
    SelectKBest(f_classif, k=50),
    LogisticRegression(max_iter=1000),
)

predicted = cross_val_predict(
    model,
    expression[known],
    treated[known],
    groups=samples["cell_line"][known],
    cv=LeaveOneGroupOut(),
)
correct = (predicted == treated[known]).sum()
print(
    f"leaving out one donor at a time: {correct} of {known.sum()} right"
)

model.fit(expression[known], treated[known])
unseen = samples[~known][["cell_line", "label"]].copy()
unseen["chance_of_dex"] = model.predict_proba(expression[~known])[:, 1]
print(unseen.round(3).to_string())
