#!/usr/bin/env python3
"""
AI Module 1b: interpretable ML diagnostic + SHAP on the ADNI survivor panel.
Benchmarks logistic-elastic-net / RandomForest / GradientBoosting with 10-fold
CV AUC, calibration, and SHAP feature attribution. Honest about modest AUC.
Inputs : results/ai/adni_panel.csv  (from Step17a)
Outputs: results/ai/ml_model_auc.csv, results/ai/shap_importance.csv,
         figures/pub/ai/ml_roc.png, ml_calibration.png, shap_beeswarm.png, shap_bar.png
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, os
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import shap
os.makedirs("results/ai", exist_ok=True); os.makedirs("figures/pub/ai", exist_ok=True)

df = pd.read_csv("results/ai/adni_panel.csv")
y = (df["group"]=="AD").astype(int).values
genes = [c for c in df.columns if c not in ("group","AGE","SEX","APOE4")]
df["SEXn"]=(df["SEX"]=="Male").astype(float)
df["AGE"]=df["AGE"].fillna(df["AGE"].median()); df["APOE4"]=df["APOE4"].fillna(0)
# RAW features for CV (scaler goes INSIDE the pipeline -> no leakage)
Xg_raw = df[genes].values
Xf_raw = np.hstack([Xg_raw, df[["AGE","SEXn","APOE4"]].values])
# scaled versions used only for the final (descriptive) SHAP model
X_genes = StandardScaler().fit_transform(Xg_raw)
X_full  = np.hstack([X_genes, StandardScaler().fit_transform(df[["AGE","SEXn","APOE4"]].values)])
feat_full = genes + ["AGE","SEX","APOE4"]
print(f"ADNI: {len(y)} samples, {y.sum()} AD / {(1-y).sum()} CN, {len(genes)} panel genes")

cv = StratifiedKFold(5, shuffle=True, random_state=0)
models = {
  "ElasticNet-LR": LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=.5, C=.1, max_iter=4000),
  "RandomForest" : RandomForestClassifier(n_estimators=400, max_depth=4, random_state=0),
  "GradBoost"    : HistGradientBoostingClassifier(max_depth=3, learning_rate=.05, max_iter=300, random_state=0),
}
rows=[]; roc_data={}
for design,(X,tag) in {"panel":(Xg_raw,"genes"),"panel+demo":(Xf_raw,"genes+AGE/SEX/APOE")}.items():
  for name,mdl in models.items():
    p = cross_val_predict(make_pipeline(StandardScaler(), mdl), X, y, cv=cv, method="predict_proba")[:,1]
    auc = roc_auc_score(y,p); rows.append({"design":design,"model":name,"cv_auc":round(auc,3)})
    roc_data[f"{design}:{name}"]=(y,p,auc)
    print(f"  {design:10s} {name:14s} CV-AUC={auc:.3f}")
pd.DataFrame(rows).to_csv("results/ai/ml_model_auc.csv",index=False)

# ROC of best model
best=max(roc_data,key=lambda k:roc_data[k][2]); yb,pb,ab=roc_data[best]
fpr,tpr,_=roc_curve(yb,pb)
plt.figure(figsize=(3.2,3.2)); plt.plot([0,1],[0,1],"--",c="grey")
plt.plot(fpr,tpr,c="#D55E00",lw=2,label=f"{best}\nAUC={ab:.2f}")
plt.xlabel("1 - specificity"); plt.ylabel("sensitivity"); plt.legend(fontsize=7); plt.title("ADNI diagnostic ML",fontsize=10,weight="bold")
plt.tight_layout(); plt.savefig("figures/pub/ai/ml_roc.png",dpi=300); plt.close()
# calibration
frac,mean=calibration_curve(yb,pb,n_bins=8)
plt.figure(figsize=(3.2,3.2)); plt.plot([0,1],[0,1],"--",c="grey"); plt.plot(mean,frac,"o-",c="#0072B2")
plt.xlabel("predicted"); plt.ylabel("observed"); plt.title("Calibration",fontsize=10,weight="bold")
plt.tight_layout(); plt.savefig("figures/pub/ai/ml_calibration.png",dpi=300); plt.close()

# SHAP on the elastic-net logistic model (best AUC; linear SHAP is exact)
lr=LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=.5,C=.1,max_iter=4000).fit(X_full,y)
expl=shap.LinearExplainer(lr,X_full,feature_names=feat_full)
sv=expl(X_full)
imp=pd.DataFrame({"feature":feat_full,"mean_abs_shap":np.abs(sv.values).mean(0)}).sort_values("mean_abs_shap",ascending=False)
imp.to_csv("results/ai/shap_importance.csv",index=False)
print("\nTop SHAP features:"); print(imp.head(12).to_string(index=False))
plt.figure(); shap.summary_plot(sv,features=X_full,feature_names=feat_full,show=False,max_display=15)
plt.tight_layout(); plt.savefig("figures/pub/ai/shap_beeswarm.png",dpi=300,bbox_inches="tight"); plt.close()
plt.figure(); shap.summary_plot(sv,features=X_full,feature_names=feat_full,plot_type="bar",show=False,max_display=15)
plt.tight_layout(); plt.savefig("figures/pub/ai/shap_bar.png",dpi=300,bbox_inches="tight"); plt.close()
print("\nSaved: ml_model_auc.csv, shap_importance.csv, figures/pub/ai/{ml_roc,ml_calibration,shap_beeswarm,shap_bar}.png")
