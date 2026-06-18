import os, re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr

glmnet, pROC, base = importr("glmnet"), importr("pROC"), importr("base")
r_predict, r_coef, p2r = ro.r["predict"], ro.r["coef"], numpy2ri.py2rpy

SEED = 1
base.set_seed(SEED)

filename = "HUMAN_ES.txt"
target_file = "HIGH_AUC_PPV_FILTERED_FUNCTIONS/GO0005730_nucleolus.tsv"

data = pd.read_csv(filename, sep="\t").fillna(0)  # set missing data to 0. R's regression models don't always work well with missing data.
data = data.set_index("idr_name", drop=False)

gene2row, gene2idx = [], {}
for i in range(len(data)):
    gene = str(data["idr_name"].iloc[i]).split("_")[0]
    if gene not in gene2idx:
        gene2idx[gene] = len(gene2row)
        gene2row.append([])
    gene2row[gene2idx[gene]].append(i)

X = data.iloc[:, 2:].to_numpy(dtype=float)

initw = np.zeros(len(data))  # will store initial weights
multi = []
for g in range(len(gene2row)):
    initw[gene2row[g][0]] = 1
    if len(gene2row[g]) > 1:  # if there is more than one idr for this gene ...
        multi.append(g + 1)
        initw[gene2row[g]] = 1 / len(gene2row[g])  # downweight them all equally

targets = pd.read_csv(target_file, sep="\t").set_index("idr_name", drop=False).loc[data.index]
y = targets.iloc[:, 1].to_numpy(dtype=float)
category_name = re.sub(r"\.(txt|tsv)$", "", Path(target_file).name)


def fit_FAIDR(train_ind, initw, lam):
    w, z, pr, w1 = initw.copy(), None, None, None
    mod = glmnet.glmnet(p2r(X[train_ind, :]), p2r(y[train_ind]), family="gaussian", weights=p2r(w[train_ind]))  # intial model
    for i in range(20):
        for j in range(5):  # M-step: fit the weighted logistic regression by iterative least squares
            pr = 1 / (1 + np.exp(-np.asarray(r_predict(mod, p2r(X), s=lam)).ravel()))
            w1 = pr * (1 - pr)
            w1[w1 == 0] = np.exp(-20)  # avoid 0s
            z = np.asarray(r_predict(mod, p2r(X), s=lam)).ravel() + (y - pr) / w1
            mod = glmnet.glmnet(p2r(X[train_ind, :]), p2r(z[train_ind]), family="gaussian", weights=p2r((w1 * w)[train_ind]))
        for g in multi:  # E-step: recalculate weights given latest model
            rows = gene2row[g - 1]
            w[rows] = y[rows] * pr[rows] + (1 - y[rows]) * (1 - pr[rows])
            w[rows] = w[rows] / w[rows].sum()
    return {"mod": mod, "pr": pr, "w": w, "w1": w1, "z": z}


lam, N_splits = 0.2, 6
N_reps = int(os.environ.get("N_REPS", "100"))
Balance, Specificity = 3, 0.99
Tstats = np.zeros((X.shape[1] + 1, 1))
PostIDR, Pred, geney = np.zeros(len(data)), np.zeros(len(gene2row)), np.zeros(len(gene2row))

IDR_level_pr = pd.DataFrame({"idr_name": targets["idr_name"].values, "pr": 0.0})
IDR_level_pred_counts = pd.DataFrame({"idr_name": targets["idr_name"].values, "pred": 0.0})
classi_info = np.zeros((1, 7), dtype=object)
classi_info[0, 0] = category_name
for g in range(len(gene2row)):
    geney[g] = y[gene2row[g][0]]

auc_cv_reps = np.zeros(N_reps)
roc_cv_list = [None] * N_reps
for rep in range(N_reps):
    gene_split_assign = np.array(base.cut(base.sample(ro.IntVector(list(range(1, len(gene2row) + 1)))), breaks=N_splits, labels=ro.r("FALSE")), dtype=int)
    unbiased_pr = np.full(len(data), np.nan)
    for this_split in range(1, N_splits + 1):
        in_split = np.where(gene_split_assign == this_split)[0] + 1
        not_split = np.where(gene_split_assign != this_split)[0] + 1
        test_pos = np.intersect1d(np.where(y == 1)[0], np.concatenate([gene2row[g - 1] for g in in_split]))
        neg_genes = np.array(base.sample(ro.IntVector(np.intersect1d(np.where(geney == 0)[0] + 1, in_split)), size=int(geney[in_split - 1].sum() * Balance)), dtype=int)
        test_neg = np.concatenate([gene2row[g - 1] for g in neg_genes])
        train_pos = np.intersect1d(np.where(y == 1)[0], np.concatenate([gene2row[g - 1] for g in not_split]))
        neg_genes = np.array(base.sample(ro.IntVector(np.intersect1d(np.where(geney == 0)[0] + 1, not_split)), size=int(geney[not_split - 1].sum() * Balance)), dtype=int)
        train_neg = np.concatenate([gene2row[g - 1] for g in neg_genes])
        train = np.concatenate([train_pos, train_neg])
        split_fit = fit_FAIDR(train, initw, lam)
        split_pr = 1 / (1 + np.exp(-np.asarray(r_predict(split_fit["mod"], p2r(X), s=lam)).ravel()))
        unbiased_pr[test_pos], unbiased_pr[test_neg] = split_pr[test_pos], split_pr[test_neg]
    mask = ~np.isnan(unbiased_pr)
    roc_cv_list[rep] = pROC.roc(p2r(y[mask].astype(float)), p2r(unbiased_pr[mask].astype(float)))
    auc_cv_reps[rep] = float(np.asarray(pROC.auc(roc_cv_list[rep])).ravel()[0])

classi_info[0, 4] = auc_cv_reps.mean()
if os.environ.get("SKIP_PLOT", "") != "1":
    fig, ax = plt.subplots()
    for rep in range(N_reps):
        roc_cv = roc_cv_list[rep]
        ax.plot(1 - np.asarray(roc_cv.rx2("specificities"), dtype=float), np.asarray(roc_cv.rx2("sensitivities"), dtype=float), color=(0.23, 0.45, 0.71, 0.2), linewidth=0.8)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"CV ROC ({N_reps} reps), mean AUC = {auc_cv_reps.mean():.3f}")
    fig.tight_layout()
    plt.show()

for rep in range(N_reps):
    train_pos = np.where(y == 1)[0]  # use all the positives
    gene_neg = np.array(base.sample(ro.IntVector(np.where(geney == 0)[0] + 1), size=int(geney.sum() * Balance)), dtype=int) - 1
    train = np.concatenate([train_pos, np.concatenate([gene2row[g] for g in gene_neg])])
    fit = fit_FAIDR(train, initw, lam)
    mod, pr, w, w1, z = fit["mod"], fit["pr"], fit["w"], fit["w1"], fit["z"]
    for g in multi:  # final weights using a formula that doesn't depend on y
        rows = gene2row[g - 1]
        w[rows] = (pr[rows] + 0.05) / (pr[rows] + 0.05).sum()
    PostIDR = w
    for g in range(len(gene2row)):
        Pred[g] = np.sum(pr[gene2row[g]] * w[gene2row[g]])
    cols = np.abs(np.asarray(r_coef(mod, s=lam)).ravel()[1:]) > 0
    if cols.sum() > 0:
        subX = X[:, cols]
        ro.r.assign("z_glm", z); ro.r.assign("subX_glm", subX); ro.r.assign("w_glm", w * w1)
        tstats = np.asarray(ro.r("summary(glm(z_glm ~ subX_glm, weights=w_glm))$coefficients[,3]"), dtype=float).ravel()
        Tstats[np.concatenate([[0], np.where(cols)[0] + 1]), 0] = tstats
    roc_train = pROC.roc(p2r(y[train].astype(float)), p2r(pr[train].astype(float)))
    specificities, thresholds = np.asarray(roc_train.rx2("specificities"), dtype=float), np.asarray(roc_train.rx2("thresholds"), dtype=float)
    thresh = thresholds[np.where(specificities > Specificity)[0].min()]
    classi_info[0, 2] = float(classi_info[0, 2]) + len(np.intersect1d(np.where(y == 1)[0], np.where(pr > thresh)[0])) / N_reps
    classi_info[0, 3] = float(classi_info[0, 3]) + len(np.intersect1d(np.where(y == 0)[0], np.where(pr > thresh)[0])) / N_reps
    classi_info[0, 1] = float(classi_info[0, 1]) + thresh / N_reps
    prot_idx = np.concatenate([np.where(geney == 1)[0], gene_neg])
    classi_info[0, 5] = float(classi_info[0, 5]) + float(np.asarray(pROC.auc(pROC.roc(p2r(geney[prot_idx].astype(float)), p2r(Pred[prot_idx].astype(float))))).ravel()[0]) / N_reps
    IDR_level_pr["pr"] = IDR_level_pr["pr"] + pr / N_reps
    IDR_level_pred_counts["pred"] = IDR_level_pred_counts["pred"] + (pr > thresh) / N_reps

classi_info[0, 6] = y.sum()
IDR_level_pr.columns = ["idr_name", category_name]
IDR_level_pred_counts.columns = ["idr_name", category_name]
classi_info_df = pd.DataFrame([classi_info[0]], columns=["category", "thresh", "above_thresh_in_annotated", "above_thresh_not_in_annotated", "IDR_AUC_CV", "Prot_AUC_CV", "N_annotated_IDRs"])
IDR_level_pred_counts.to_csv("IDR_above_cut_bal_control_no_cv.txt", sep="\t", index=False)
classi_info_df.to_csv("IDR_classification_bal_control_no_cv.txt", sep="\t", index=False)
