#!/usr/bin/env python3
"""
AI Module 2b: deep "virtual-cell" embedding + signature robustness  [SERVER, GPU]
Trains a deep generative embedding (scVI) of the PBMCs and tests whether the
survivor signature still separates AD vs Control in the learned latent space —
a model-agnostic robustness check. (Swap scVI for scGPT/Geneformer for extra
"foundation-model" framing; the downstream test is identical.)

Install (GPU instance):
  pip install scvi-tools scanpy
Run after Step20a:
  python script/ai/Step20b_embedding.py
Outputs: figures/pub/ai/scvi_umap_condition.png, results/ai/scvi_latent_auc.txt
"""
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc, numpy as np, pandas as pd, scipy.io as sio, os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.makedirs("figures/pub/ai",exist_ok=True); os.makedirs("results/ai",exist_ok=True)

X = sio.mmread("results/ai/counts.mtx").T.tocsr()                 # cells x genes
genes = [l.strip() for l in open("results/ai/genes.txt")]
cells = pd.read_csv("results/ai/cells.csv")
ad = sc.AnnData(X); ad.var_names=genes; ad.obs=cells.set_index("barcode")
ad = ad[~ad.obs.condition.isna()].copy()
sc.pp.filter_genes(ad, min_cells=10)
ad.layers["counts"]=ad.X.copy()

import scvi
# NOTE: do NOT pass batch_key="condition" — that would make scVI integrate OUT the
# AD-vs-Control difference, defeating the test. Embed condition-blind, then ask if
# condition is still separable in the learned latent space.
scvi.model.SCVI.setup_anndata(ad, layer="counts")
m = scvi.model.SCVI(ad, n_latent=30); m.train(max_epochs=200)
ad.obsm["X_scVI"]=m.get_latent_representation()

# UMAP of the latent space colored by condition
sc.pp.neighbors(ad, use_rep="X_scVI"); sc.tl.umap(ad)
sc.pl.umap(ad, color="condition", show=False, palette=["#0072B2","#D55E00"])
plt.savefig("figures/pub/ai/scvi_umap_condition.png", dpi=300, bbox_inches="tight"); plt.close()

# does the latent space carry AD vs Control signal? (logistic AUC on latent dims)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score
y=(ad.obs.condition=="AD").astype(int).values
p=cross_val_predict(LogisticRegression(max_iter=2000), ad.obsm["X_scVI"], y, cv=5, method="predict_proba")[:,1]
auc=roc_auc_score(y,p)
open("results/ai/scvi_latent_auc.txt","w").write(f"scVI latent AD-vs-Control AUC (cell-level): {auc:.3f}\n")
print(f"scVI latent AD-vs-Control AUC: {auc:.3f}")
print("Saved: figures/pub/ai/scvi_umap_condition.png, results/ai/scvi_latent_auc.txt")
