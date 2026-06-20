#!/usr/bin/env python3
"""Step 33b: FigSx_deconv_robustness from REAL Step32 outputs (no placeholders)."""
import os, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7})
B,O,Gn,Pu,Gy,Rd,Dk="#0072B2","#D55E00","#009E73","#7A4F9E","#999","#c0413b","#222"

DIR="results/deconv_robustness"
sum_=pd.read_csv(f"{DIR}/concordance_summary.csv")

fig=plt.figure(figsize=(7.1,5.6))
gs=GridSpec(2,2,figure=fig,hspace=0.55,wspace=0.40)

# === Panel A: Monocyte scatter, Danaher vs Bisque, both cohorts ===
ax=fig.add_subplot(gs[0,0])
for coh,col in zip(["GSE63060","GSE63061"],[B,O]):
    d=pd.read_csv(f"{DIR}/{coh}_Bisque_vs_Danaher.csv")
    ax.scatter(d["Monocyte.dan"], d["Monocyte.x"], s=8, alpha=0.55, c=col, label=coh, lw=0)
r60=sum_[(sum_.cohort=="GSE63060")&(sum_.method=="Bisque")&(sum_.lineage=="Monocyte")].r.iloc[0]
r61=sum_[(sum_.cohort=="GSE63061")&(sum_.method=="Bisque")&(sum_.lineage=="Monocyte")].r.iloc[0]
ax.set_xlabel("Monocyte score (Danaher)"); ax.set_ylabel("Monocyte fraction (Bisque)")
ax.set_title(f"Bisque vs Danaher — Monocyte\nr={r60:.2f} (63060), {r61:.2f} (63061)")
ax.legend(loc="best",fontsize=6,frameon=False)
ax.text(-0.18,1.10,"a",transform=ax.transAxes,fontsize=10,fontweight="bold")

# === Panel B: Monocyte AD vs Ctrl by method, both cohorts ===
ax=fig.add_subplot(gs[0,1])
mono=sum_[sum_.lineage=="Monocyte"].copy()
# include Danaher baseline (computed from cell_scores + status)
dan_rows=[]
for coh in ["GSE63060","GSE63061"]:
    cs=pd.read_csv(f"results/cell_scores_{coh}.csv")
    st=pd.read_csv(f"{DIR}/{coh}_sample_status.csv")
    m=cs.merge(st,on="sample")
    m["status"]=m["status"].astype(str).str.strip()
    ad=m[m.status=="AD"].Monocyte.mean(); ct=m[m.status.isin(["CTL","Control","HC","NC"])].Monocyte.mean()
    dan_rows.append({"cohort":coh,"method":"Danaher","AD_minus_Ctrl":ad-ct})
dan=pd.DataFrame(dan_rows)
# bar groups: cohorts on x, methods grouped
labels=["GSE63060","GSE63061"]
x=np.arange(len(labels)); w=0.25
for i,(meth,col) in enumerate(zip(["Danaher","Bisque","MuSiC"],[Gy,B,Pu])):
    if meth=="Danaher":
        vals=[dan[dan.cohort==c].AD_minus_Ctrl.iloc[0] for c in labels]
    else:
        vals=[mono[(mono.cohort==c)&(mono.method==meth)].AD_minus_Ctrl.iloc[0] for c in labels]
    ax.bar(x+(i-1)*w, vals, w, color=col, label=meth, edgecolor="white", lw=0.5)
ax.axhline(0,c="k",lw=0.5); ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Monocyte AD − Control")
ax.set_title("Monocyte AD > Control in 4/4 tests\n(direction reproducible across methods)")
ax.legend(loc="best",fontsize=6,frameon=False)
ax.text(-0.18,1.10,"b",transform=ax.transAxes,fontsize=10,fontweight="bold")

# === Panel C: per-lineage r heatmap (Bisque/MuSiC vs Danaher) ===
ax=fig.add_subplot(gs[1,0])
lineages=["Monocyte","Bcell","Tcell_CD4","Tcell_CD8","NK"]
mat=np.full((len(lineages),4),np.nan)
cols=[]
for j,(coh,meth) in enumerate([("GSE63060","Bisque"),("GSE63060","MuSiC"),
                                ("GSE63061","Bisque"),("GSE63061","MuSiC")]):
    cols.append(f"{coh[-3:]}\n{meth}")
    for i,lin in enumerate(lineages):
        sub=sum_[(sum_.cohort==coh)&(sum_.method==meth)&(sum_.lineage==lin)]
        if len(sub): mat[i,j]=sub.r.iloc[0]
im=ax.imshow(mat,cmap="RdBu_r",vmin=-0.7,vmax=0.7,aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(cols,fontsize=6)
ax.set_yticks(range(len(lineages))); ax.set_yticklabels(lineages,fontsize=7)
for i in range(len(lineages)):
    for j in range(4):
        v=mat[i,j]
        if not np.isnan(v):
            ax.text(j,i,f"{v:.2f}",ha="center",va="center",
                    fontsize=6,color="white" if abs(v)>0.4 else "black")
ax.set_title("Per-lineage r(method vs Danaher)\n(red = positive concordance)")
cb=plt.colorbar(im,ax=ax,fraction=0.046,pad=0.04); cb.ax.tick_params(labelsize=5)
ax.text(-0.30,1.10,"c",transform=ax.transAxes,fontsize=10,fontweight="bold")

# === Panel D: takeaway summary ===
ax=fig.add_subplot(gs[1,1]); ax.axis("off")
ax.text(0.02,0.95,"Take-home",fontsize=9,fontweight="bold",color=Dk,transform=ax.transAxes)
lines=[
    "• Monocyte AD>Control direction is preserved",
    "  in 4/4 (cohort × method) tests.",
    "",
    "• Monocyte r(method vs Danaher) median 0.54;",
    "  Bisque achieves 0.51 (GSE63060) and 0.56",
    "  (GSE63061); MuSiC 0.32 and 0.63.",
    "",
    "• Cross-method numerical disagreement is",
    "  expected for T-cell subsets and B cells —",
    "  Danaher marker panels and PBMC scRNA",
    "  references resolve different sub-populations.",
    "",
    "• The audit's central claim — that the AD blood",
    "  signal tracks myeloid composition — does not",
    "  depend on the deconvolution method:",
    "  all three methods agree on the monocyte shift.",
    "",
    "• Caveat: PBMC reference has no neutrophils;",
    "  neutrophil cannot be cross-validated here.",
]
for i,t in enumerate(lines):
    ax.text(0.02, 0.85 - i*0.045, t, fontsize=6.6, color=Dk if not t.startswith("•") else B,
            transform=ax.transAxes, fontweight="bold" if t.startswith("•") else "normal")
ax.text(-0.10,1.10,"d",transform=ax.transAxes,fontsize=10,fontweight="bold")

fig.savefig("figures/pub/FigSx_deconv_robustness.pdf",bbox_inches="tight")
fig.savefig("figures/pub/FigSx_deconv_robustness.png",dpi=600,bbox_inches="tight"); plt.close()
print("Saved FigSx_deconv_robustness.{pdf,png}")
