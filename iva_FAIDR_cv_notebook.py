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

FEATURE_DIR = Path(__file__).resolve().parent / "feature_matrices"

FEATURE_DEFAULTS = {
    "evolutionary": FEATURE_DIR / "HUMAN_ES.txt",
    "signature": FEATURE_DIR / "FS_UP000005640_9606_SPOTD_MIN_30AA.txt",
    "rohit": FEATURE_DIR / "Rohit's_features.tsv",
}


def uniprot_from_idr(idr_name):
    s = str(idr_name)
    m = re.search(r"\|([^|]+)\|", s)
    return m.group(1) if m else s.split("_")[0]


def _load_features(path):
    data = pd.read_csv(path, sep="\t").fillna(0)  # set missing data to 0. R's regression models don't always work well with missing data.
    if "idr_name" not in data.columns:
        raise ValueError(f"{path}: expected idr_name column")
    data = data[data["idr_name"].notna() & data["idr_name"].astype(str).str.len().gt(0)]
    return data.set_index("idr_name", drop=False)


def run_faidr_cv(
    target_file,
    feature_type="evolutionary",
    features_path=None,
    specificity=0.99,
    n_reps=100,
    consistency_cutoff=0.5,
    show_progress=False,
):
    target_file = Path(target_file)
    filename = Path(features_path) if features_path is not None else Path(FEATURE_DEFAULTS[feature_type])

    N_reps = int(n_reps if os.environ.get("N_REPS", "") == "" else os.environ.get("N_REPS", n_reps))
    Balance, Specificity = 3, float(specificity)
    consistency_cutoff = float(consistency_cutoff)

    base.set_seed(1)

    data = _load_features(filename)

    gene2row, gene2idx = [], {}
    for i in range(len(data)):
        gene = str(data["idr_name"].iloc[i]).split("_")[0]
        if gene not in gene2idx:
            gene2idx[gene] = len(gene2row)
            gene2row.append([])
        gene2row[gene2idx[gene]].append(i)

    idr_index = data.index
    X = data.iloc[:, 1:].to_numpy(dtype=float)
    del data

    targets = pd.read_csv(target_file, sep="\t").set_index("idr_name", drop=False).loc[idr_index]
    y = targets.iloc[:, 1].to_numpy(dtype=float)
    category_name = re.sub(r"\.(txt|tsv)$", "", target_file.name)

    ro.r.assign("X_faidr", ro.r.matrix(ro.FloatVector(X.ravel()), nrow=X.shape[0], ncol=X.shape[1], byrow=True))
    ro.r.assign("y_faidr", ro.FloatVector(y))
    X_r = ro.r["X_faidr"]

    initw = np.zeros(len(X))  # will store initial weights
    multi = []
    for g in range(len(gene2row)):
        initw[gene2row[g][0]] = 1
        if len(gene2row[g]) > 1:  # if there is more than one idr for this gene ...
            multi.append(g + 1)
            initw[gene2row[g]] = 1 / len(gene2row[g])  # downweight them all equally

    def _r_train_rows(train_ind):
        ro.r.assign("i_tr", ro.IntVector((np.asarray(train_ind, dtype=int) + 1).tolist()))
        return ro.r("X_faidr[i_tr,,drop=FALSE]"), ro.r("y_faidr[i_tr]")

    def fit_FAIDR(train_ind, initw, lam):
        X_train, y_train = _r_train_rows(train_ind)
        w, z, pr, w1 = initw.copy(), None, None, None
        mod = glmnet.glmnet(X_train, y_train, family="gaussian", weights=p2r(w[train_ind]))  # intial model
        for i in range(20):
            for j in range(5):  # M-step: fit the weighted logistic regression by iterative least squares
                pred_lin = np.asarray(r_predict(mod, X_r, s=lam)).ravel()
                pr = 1 / (1 + np.exp(-pred_lin))
                w1 = pr * (1 - pr)
                w1[w1 == 0] = np.exp(-20)  # avoid 0s
                z = pred_lin + (y - pr) / w1
                mod = glmnet.glmnet(X_train, p2r(z[train_ind]), family="gaussian", weights=p2r((w1 * w)[train_ind]))
            for g in multi:  # E-step: recalculate weights given latest model
                rows = gene2row[g - 1]
                w[rows] = y[rows] * pr[rows] + (1 - y[rows]) * (1 - pr[rows])
                w[rows] = w[rows] / w[rows].sum()
        return {"mod": mod, "pr": pr, "w": w, "w1": w1, "z": z}

    lam, N_splits = 0.2, 6
    Tstats = np.zeros((X.shape[1] + 1, 1))
    PostIDR, Pred, geney = np.zeros(len(X)), np.zeros(len(gene2row)), np.zeros(len(gene2row))

    IDR_level_pr = pd.DataFrame({"idr_name": targets["idr_name"].values, "pr": 0.0})
    IDR_level_pred_counts = pd.DataFrame({"idr_name": targets["idr_name"].values, "pred": 0.0})
    classi_info = np.zeros((1, 7), dtype=object)
    classi_info[0, 0] = category_name
    for g in range(len(gene2row)):
        geney[g] = y[gene2row[g][0]]

    auc_cv_reps = np.zeros(N_reps)
    roc_cv_curves = []
    cv_reps = range(N_reps)
    if show_progress:
        from tqdm.auto import tqdm
        cv_reps = tqdm(cv_reps, desc="Step 1/2 — Cross-validation (ROC / AUC)", unit="run")
    for rep in cv_reps:
        gene_split_assign = np.array(base.cut(base.sample(ro.IntVector(list(range(1, len(gene2row) + 1)))), breaks=N_splits, labels=ro.r("FALSE")), dtype=int)
        unbiased_pr = np.full(len(X), np.nan)
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
            split_pr = 1 / (1 + np.exp(-np.asarray(r_predict(split_fit["mod"], X_r, s=lam)).ravel()))
            unbiased_pr[test_pos], unbiased_pr[test_neg] = split_pr[test_pos], split_pr[test_neg]
        mask = ~np.isnan(unbiased_pr)
        roc_cv = pROC.roc(p2r(y[mask].astype(float)), p2r(unbiased_pr[mask].astype(float)), quiet=True)
        auc_cv_reps[rep] = float(np.asarray(pROC.auc(roc_cv)).ravel()[0])
        fpr = 1 - np.asarray(roc_cv.rx2("specificities"), dtype=float)
        tpr = np.asarray(roc_cv.rx2("sensitivities"), dtype=float)
        order = np.argsort(fpr)
        roc_cv_curves.append((fpr[order], tpr[order]))

    fig, ax = plt.subplots()
    for fpr, tpr in roc_cv_curves:
        ax.plot(fpr, tpr, color=(0.23, 0.45, 0.71, 0.2), linewidth=0.8)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"CV ROC ({N_reps} reps), mean AUC = {auc_cv_reps.mean():.3f}")
    fig.tight_layout()
    plt.show()

    final_reps = range(N_reps)
    if show_progress:
        from tqdm.auto import tqdm
        final_reps = tqdm(final_reps, desc="Step 2/2 — Final training (IDR predictions)", unit="run")
    for rep in final_reps:
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
        roc_train = pROC.roc(p2r(y[train].astype(float)), p2r(pr[train].astype(float)), quiet=True)
        classi_info[0, 4] = float(classi_info[0, 4]) + float(np.asarray(pROC.auc(roc_train)).ravel()[0]) / N_reps
        specificities, thresholds = np.asarray(roc_train.rx2("specificities"), dtype=float), np.asarray(roc_train.rx2("thresholds"), dtype=float)
        thresh = thresholds[np.where(specificities > Specificity)[0].min()]
        classi_info[0, 2] = float(classi_info[0, 2]) + len(np.intersect1d(np.where(y == 1)[0], np.where(pr > thresh)[0])) / N_reps
        classi_info[0, 3] = float(classi_info[0, 3]) + len(np.intersect1d(np.where(y == 0)[0], np.where(pr > thresh)[0])) / N_reps
        classi_info[0, 1] = float(classi_info[0, 1]) + thresh / N_reps
        prot_idx = np.concatenate([np.where(geney == 1)[0], gene_neg])
        classi_info[0, 5] = float(classi_info[0, 5]) + float(np.asarray(pROC.auc(pROC.roc(p2r(geney[prot_idx].astype(float)), p2r(Pred[prot_idx].astype(float)), quiet=True))).ravel()[0]) / N_reps
        IDR_level_pr["pr"] = IDR_level_pr["pr"] + pr / N_reps
        IDR_level_pred_counts["pred"] = IDR_level_pred_counts["pred"] + (pr > thresh) / N_reps

    classi_info[0, 6] = y.sum()
    classi_info_df = pd.DataFrame([classi_info[0]], columns=["category", "thresh", "above_thresh_in_annotated", "above_thresh_not_in_annotated", "IDR_AUC_train", "Prot_AUC_train", "N_annotated_IDRs"])

    idr_results = pd.DataFrame({
        "uniprot_accession": [uniprot_from_idr(x) for x in IDR_level_pr["idr_name"]],
        "idr_name": IDR_level_pr["idr_name"],
        "mean_probability": IDR_level_pr["pr"],
        "consistency_fraction": IDR_level_pred_counts["pred"],
    })

    out_dir = Path("FAIDR_output") / category_name
    out_dir.mkdir(parents=True, exist_ok=True)
    classi_info_df.to_csv(out_dir / "run_summary.tsv", sep="\t", index=False)
    idr_results.to_csv(out_dir / "idr_predictions.tsv", sep="\t", index=False)
    idr_results_filtered = idr_results[idr_results["consistency_fraction"] >= consistency_cutoff]
    filtered_path = out_dir / f"idr_predictions_consistency_{consistency_cutoff:g}.tsv"
    idr_results_filtered.to_csv(filtered_path, sep="\t", index=False)

    return {
        "classi_info": classi_info_df,
        "idr_results": idr_results,
        "idr_results_filtered": idr_results_filtered,
        "out_dir": out_dir,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--target-file", default="annotation_targets/GO0005730_nucleolus.tsv")
    parser.add_argument("--feature-type", choices=sorted(FEATURE_DEFAULTS), default="evolutionary")
    parser.add_argument("--features-path", default=None)
    parser.add_argument("--specificity", type=float, default=0.99)
    parser.add_argument("--n-reps", type=int, default=100)
    parser.add_argument("--consistency-cutoff", type=float, default=0.5)
    args = parser.parse_args()
    run_faidr_cv(
        target_file=args.target_file,
        feature_type=args.feature_type,
        features_path=args.features_path,
        specificity=args.specificity,
        n_reps=args.n_reps,
        consistency_cutoff=args.consistency_cutoff,
    )
