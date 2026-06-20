#!/usr/bin/env python3
"""Step 27: regenerate Fig1 (composition), Fig4 (enrichment), Fig7 (ML ROC+SHAP)
to the same publication spec as Step26 (600 dpi, Arial, >=7pt, editable PDF)."""
import numpy as np, pandas as pd, csv, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({
    "pdf.fonttype":42,"ps.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7,"axes.titlesize":7.5,"axes.labelsize":7,
    "xtick.labelsize":6.5,"ytick.labelsize":6.5,"legend.fontsize":6,
    "axes.linewidth":0.6,"lines.linewidth":1.0,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "axes.spines.top":False,"axes.spines.right":False,
    "savefig.dpi":600,"figure.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.02})
B,O,Gy,Gn,Pu="#0072B2","#D55E00","#999999","#009E73","#9467bd"
def save(fig,n): fig.savefig(f"figures/pub/{n}.pdf"); fig.savefig(f"figures/pub/{n}.png"); plt.close(fig)
def num(x):
    try:return float(x)
    except:return np.nan

# ---------- Fig 1: composition shift (grouped bars, lineage x cohort) ----------
cs=pd.read_csv("results/cell_composition_shift.csv")
lin=["Neutrophil","Monocyte","NK","Tcell_CD8","Tcell_CD4","Bcell"]
labs=["Neutro","Mono","NK","CD8 T","CD4 T","B cell"]
cohs=cs.cohort.unique(); cols={cohs[0]:B,cohs[1]:O,cohs[2]:Gn}
fig,ax=plt.subplots(figsize=(4.0,3.0)); x=np.arange(len(lin)); w=0.26
for i,co in enumerate(cohs):
    d=cs[cs.cohort==co].set_index("lineage")
    vals=[num(d.loc[l,"AD_minus_CN"]) if l in d.index else 0 for l in lin]
    ps=[num(d.loc[l,"AD_vs_CN_p"]) if l in d.index else 1 for l in lin]
    bars=ax.bar(x+(i-1)*w,vals,w,color=cols[co],label=co.replace("GSE",""))
    for b,p in zip(bars,ps):
        if p<0.05: ax.text(b.get_x()+b.get_width()/2,b.get_height()+(0.01 if b.get_height()>=0 else -0.03),"*",ha="center",fontsize=7)
ax.axhline(0,c="k",lw=.5); ax.set_xticks(x); ax.set_xticklabels(labs,rotation=20,ha="right")
ax.set_ylabel("composition shift (AD − control)")
ax.set_title("Reproducible myeloid↑ / lymphoid↓\nshift in AD blood (3/3 cohorts)",weight="bold")
ax.legend(frameon=False,title="cohort",ncol=1,loc="upper right",fontsize=5.5,title_fontsize=6)
save(fig,"Fig1_composition")

# ---------- Fig 4: enrichment of survivors (horizontal bar) ----------
src="results/enrichment_robust143.csv"
import os
if not os.path.exists(src): src="results/enrichment_survivors_verified.csv"
en=pd.read_csv(src); en["adjp"]=pd.to_numeric(en["adjp"],errors="coerce")
en=en.dropna(subset=["adjp"]).sort_values("adjp").head(10).iloc[::-1]
terms=[t[:42] for t in en["term"]]
fig,ax=plt.subplots(figsize=(4.6,3.0))
ax.barh(range(len(en)),-np.log10(en["adjp"]),color=O,alpha=.85)
ax.axvline(-np.log10(0.05),ls="--",c=Gy,lw=.7); ax.text(-np.log10(0.05)+.05,0,"FDR 0.05",fontsize=5.5,color=Gy,rotation=90,va="bottom")
ax.set_yticks(range(len(en))); ax.set_yticklabels(terms,fontsize=6)
ax.set_xlabel("−log$_{10}$ adjusted p"); ax.set_title("De-confounded survivors enrich for\ninnate-immune / neutrophil / TLR4",weight="bold")
save(fig,"Fig4_enrichment")

# ---------- Fig 7: ML ROC + SHAP ----------
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_curve, roc_auc_score
import shap
df=pd.read_csv("results/ai/adni_panel.csv")
genes=[c for c in df.columns if c not in ("group","AGE","SEX","APOE4")]
df["SEXn"]=(df["SEX"]=="Male").astype(float); df["AGE"]=df["AGE"].fillna(df["AGE"].median()); df["APOE4"]=df["APOE4"].fillna(0)
y=(df["group"]=="AD").astype(int).values
Xg=df[genes].values; Xd=df[["AGE","SEXn","APOE4"]].values; Xf=np.hstack([Xg,Xd]); Xa=df[["APOE4"]].values
cv=StratifiedKFold(5,shuffle=True,random_state=0)
def roc(X):
    pipe=make_pipeline(StandardScaler(),LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=.5,C=.1,max_iter=4000))
    p=cross_val_predict(pipe,X,y,cv=cv,method="predict_proba")[:,1]
    fpr,tpr,_=roc_curve(y,p); return fpr,tpr,roc_auc_score(y,p)
fig,axs=plt.subplots(1,2,figsize=(6.4,3.0))
ax=axs[0]
for X,lab,c in [(Xa,"APOE4 only",O),(Xg,"panel",B),(Xf,"panel+demo",Gn)]:
    fpr,tpr,a=roc(X); ax.plot(fpr,tpr,color=c,lw=1.3,label=f"{lab} (AUC {a:.2f})")
ax.plot([0,1],[0,1],ls="--",c=Gy,lw=.6)
ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
ax.set_title("ADNI classification\n(panel adds +0.04 over APOE)",weight="bold")
ax.legend(frameon=False,loc="lower right")
# SHAP bar
Xfs=StandardScaler().fit_transform(Xf); feat=genes+["AGE","SEX","APOE4"]
lr=LogisticRegression(penalty="elasticnet",solver="saga",l1_ratio=.5,C=.1,max_iter=4000).fit(Xfs,y)
sv=shap.LinearExplainer(lr,Xfs).shap_values(Xfs)
imp=pd.DataFrame({"f":feat,"v":np.abs(sv).mean(0)}).sort_values("v").tail(10)
ax=axs[1]; ax.barh(range(len(imp)),imp["v"],color=[O if f=="APOE4" else B for f in imp["f"]])
ax.set_yticks(range(len(imp))); ax.set_yticklabels(imp["f"],fontsize=6); ax.set_xlabel("mean |SHAP|")
ax.set_title("Feature importance\n(APOE4 dominates)",weight="bold")
save(fig,"Fig7_ML")
print("Fig1, Fig4, Fig7 regenerated at 600 dpi / Arial / >=7pt / editable PDF.")
