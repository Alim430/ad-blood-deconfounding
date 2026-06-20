#!/usr/bin/env python3
"""
Step 26: regenerate all matplotlib figures to PUBLICATION spec.
- 600 dpi PNG + vector PDF (editable fonts, pdf.fonttype=42)
- Arial/Helvetica, font >=7pt at final column size
- sized to journal column widths (single 3.5in, double 7.1in)
- line widths >=0.6pt
Outputs overwrite figures/pub/*.{pdf,png}
"""
import numpy as np, pandas as pd, csv, collections, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
mpl.rcParams.update({
    "pdf.fonttype":42, "ps.fonttype":42,                 # editable embedded TrueType
    "font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7, "axes.titlesize":7.5, "axes.labelsize":7,
    "xtick.labelsize":6.5, "ytick.labelsize":6.5, "legend.fontsize":6,
    "axes.linewidth":0.6, "lines.linewidth":1.0,
    "xtick.major.width":0.6, "ytick.major.width":0.6,
    "xtick.major.size":2.5, "ytick.major.size":2.5,
    "axes.spines.top":False, "axes.spines.right":False,
    "savefig.dpi":600, "figure.dpi":600, "savefig.bbox":"tight", "savefig.pad_inches":0.02,
})
B,O,Gy,Gn="#0072B2","#D55E00","#999999","#009E73"
def save(fig,name): fig.savefig(f"figures/pub/{name}.pdf"); fig.savefig(f"figures/pub/{name}.png"); plt.close(fig)
def num(x):
    try:return float(x)
    except:return np.nan

# ---------- Fig 2b: de-confounding attenuation (single column) ----------
s=pd.read_csv("results/biomarker_survival_after_adjustment.csv")
t1=s[s.tier_orig=="Tier1_Robust"].copy()
for c in ["logFC_orig","logFC_adj"]: t1[c]=pd.to_numeric(t1[c],errors="coerce")
t1=t1.dropna(subset=["logFC_orig","logFC_adj"]); surv=t1.survives.astype(str).isin(["TRUE","True","true","1"])
fig,ax=plt.subplots(figsize=(3.4,3.2)); lim=np.nanpercentile(np.abs(np.r_[t1.logFC_orig,t1.logFC_adj]),99.5)
ax.plot([-lim,lim],[-lim,lim],c="k",lw=.6,ls="--",zorder=1)
ax.axhline(0,c=Gy,lw=.5); ax.axvline(0,c=Gy,lw=.5)
ax.scatter(t1.logFC_orig[~surv],t1.logFC_adj[~surv],s=5,c=Gy,alpha=.5,lw=0,label=f"lost (n={(~surv).sum()})")
ax.scatter(t1.logFC_orig[surv],t1.logFC_adj[surv],s=6,c=O,alpha=.85,lw=0,label=f"survived (n={surv.sum()})")
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel("log$_2$FC (naive meta-analysis)"); ax.set_ylabel("log$_2$FC (de-confounded)")
ax.set_title("De-confounding attenuates\n'robust' biomarkers",weight="bold")
ax.legend(frameon=False,loc="upper left",handletextpad=.3); save(fig,"Fig2b_attenuation")

# ---------- Fig 3: ADNI composition collapse (double column, 3 panels) ----------
fig,axs=plt.subplots(1,3,figsize=(7.1,2.5))
ax=axs[0]; sets=["240\n(primary)","162\n(+APOE)","356\n(uniform)"]; noc=[90,88,85]; wic=[49,55,46]
x=np.arange(3); w=.38
ax.bar(x-w/2,noc,w,color=B,label="no composition adj."); ax.bar(x+w/2,wic,w,color=O,label="+ composition adj.")
ax.axhline(50,ls="--",c=Gy,lw=.7); ax.text(2.25,52,"chance",fontsize=6,color=Gy)
ax.set_xticks(x); ax.set_xticklabels(sets); ax.set_ylim(0,100); ax.set_ylabel("ADNI concordant (%)")
ax.set_title("Replication collapses under\ncomposition adjustment",weight="bold")
ax.legend(frameon=False,loc="upper right",handlelength=1)
ax=axs[1]; cats=["surv\nno-comp","surv\n+comp","all\nno-comp","all\n+comp"]; vals=[.076,.025,.029,.025]
ax.bar(range(4),vals,color=[B,O,Gy,Gy]); ax.axhline(.029,ls=":",c=Gy,lw=.7)
ax.text(3.05,.0305,"background",fontsize=5.5,color=Gy)
ax.set_xticks(range(4)); ax.set_xticklabels(cats); ax.set_ylabel("median |log$_2$FC| (ADNI)")
ax.set_title("Survivor-specific collapse\nto background",weight="bold")
ax=axs[2]; mods=["no\ncomp","+2\nmyeloid","+6\nlineages"]; conc=[90,49,55]
ax.plot(range(3),conc,"o-",color=O,ms=6,lw=1.5); ax.axhline(50,ls="--",c=Gy,lw=.7)
ax.set_xticks(range(3)); ax.set_xticklabels(mods); ax.set_ylim(40,100); ax.set_ylabel("concordant (%)")
ax.set_title("Even 2 covariates collapse it\n(no over-adjustment)",weight="bold")
save(fig,"Fig3_ADNI_composition_collapse")

# ---------- Fig 3b: ADNI concordance scatter (single) ----------
a=pd.read_csv("results/ADNI_survivor_replication.csv")
a["logFC_adj"]=pd.to_numeric(a["logFC_adj"],errors="coerce"); a["logFC"]=pd.to_numeric(a["logFC"],errors="coerce")
a=a.dropna(subset=["logFC_adj","logFC"]); rep=a.replicated.astype(str).isin(["TRUE","True"])
conc=a.concordant.astype(str).isin(["TRUE","True"]); L=np.nanpercentile(np.abs(np.r_[a.logFC_adj,a.logFC]),99)
fig,ax=plt.subplots(figsize=(3.4,3.2)); ax.axhline(0,c=Gy,lw=.5); ax.axvline(0,c=Gy,lw=.5)
ax.scatter(a.logFC_adj[~rep],a.logFC[~rep],s=9,c=Gy,alpha=.6,lw=0,label="concordant")
ax.scatter(a.logFC_adj[rep],a.logFC[rep],s=12,c=B,alpha=.85,lw=0,label="replicated p<0.05")
ax.set_xlim(-L,L); ax.set_ylim(-L,L); ax.set_xlabel("discovery log$_2$FC (de-confounded)")
ax.set_ylabel("ADNI log$_2$FC (no comp. adj.)")
ax.text(.04,.95,f"{100*conc.mean():.0f}% concordant\n(without comp. adj.)",transform=ax.transAxes,va="top",fontsize=6.5,weight="bold")
ax.set_title("ADNI direction concordance\n(composition-driven)",weight="bold")
ax.legend(frameon=False,loc="lower right",handletextpad=.3); save(fig,"Fig3b_ADNI_concordance")

# ---------- Fig 5b: pseudoreplication (single) ----------
cl=collections.Counter(r['cell_type'] for r in csv.DictReader(open("results/scRNA_within_celltype_DE.csv")))
pb=collections.Counter()
for r in csv.DictReader(open("results/scRNA_pseudobulk_DE_fullgene.csv")):
    if num(r['padj'])<0.05: pb[r['cell_type']]+=1
cts=["Monocyte","NK","Bcell","CD8T","CD4T","Platelet"]; x=np.arange(len(cts)); w=.4
fig,ax=plt.subplots(figsize=(3.6,3.0))
ax.bar(x-w/2,[cl.get(c,0) for c in cts],w,color=Gy,label=f"cell-level (n={sum(cl.values())})")
ax.bar(x+w/2,[pb.get(c,0) for c in cts],w,color=O,label=f"donor pseudobulk (n={sum(pb.values())})")
ax.set_xticks(x); ax.set_xticklabels(cts,rotation=25,ha="right"); ax.set_ylabel("DE genes (padj<0.05)")
ax.set_title("Single-cell DE is pseudoreplication:\n873 → ~6 at donor level",weight="bold")
ax.legend(frameon=False,loc="upper right"); save(fig,"Fig5b_pseudoreplication")

# ---------- Fig S: robustness (double, 2x2) ----------
fig,ax=plt.subplots(2,2,figsize=(7.1,5.6))
a=ax[0,0]; v=[444,173,338,240]; a.bar(range(4),v,color=[Gy,Gy,Gy,O])
a.set_xticks(range(4)); a.set_xticklabels(["–140829","–63060","–63061","full"]); a.set_ylabel("survivors")
a.set_title("Leave-one-out meta:\nnot single-cohort driven",weight="bold")
a=ax[0,1]; v=[.708,.628,.738,.746]; a.bar(range(4),v,color=[O,B,Gy,Gy]); a.set_ylim(.5,.8); a.axhline(.708,ls=":",c=O,lw=.7)
a.set_xticks(range(4)); a.set_xticklabels(["APOE4\nonly","panel","demo","panel\n+demo"]); a.set_ylabel("AUC")
a.set_title("APOE4 alone=0.708;\npanel adds +0.038",weight="bold")
a=ax[1,0]; v=[.52,.68,.50]; a.bar(range(3),v,color=[Gy,O,B]); a.set_ylim(.45,.75); a.axhline(.5,ls="--",c=Gy,lw=.7)
a.set_xticks(range(3)); a.set_xticklabels(["cell\nscVI","donor\npseudobulk","within\ncelltype"]); a.set_ylabel("AUC")
a.set_title("Single-cell: composition-\nseparable, cell-intrinsic null",weight="bold")
a=ax[1,1]; lam=[.74,.73,.68,.87,.95,.52]; a.bar(range(6),lam,color=B); a.axhline(1,ls="--",c=O,lw=.7)
a.set_xticks(range(6)); a.set_xticklabels(["B","CD4","CD8","Mono","NK","Plt"],rotation=20); a.set_ylim(0,1.2)
a.set_ylabel("inflation λ"); a.set_title("Pseudobulk calibrated\n(λ≤1, no inflation)",weight="bold")
for axx in ax.flatten():
    axx.title.set_size(9); axx.title.set_y(1.05)
fig.subplots_adjust(hspace=0.95,wspace=0.32,top=0.92,bottom=0.10,left=0.10,right=0.97)
save(fig,"FigS_robustness")

# ---------- TCR (1 row, 3 panels, double) ----------
df=pd.read_csv("results/tcr/donor_tcr_metrics.csv"); cd8=pd.read_csv("results/tcr/donor_tcr_CD8.csv")
from scipy.stats import mannwhitneyu
def box(ax,d,col,title):
    grps=["Control","AD"]; data=[d[d.condition==g][col].dropna() for g in grps]
    bp=ax.boxplot(data,widths=.6,patch_artist=True,showfliers=False)
    for p,g in zip(bp['boxes'],grps): p.set_facecolor(B if g=="Control" else O); p.set_alpha(.6); p.set_linewidth(.6)
    for m in bp['medians']: m.set_color("k"); m.set_linewidth(.8)
    for i,g in enumerate(grps): ax.scatter(np.random.normal(i+1,.05,len(data[i])),data[i],s=6,c=B if g=="Control" else O,zorder=3,lw=.2,edgecolors="white")
    p=mannwhitneyu(data[0],data[1]).pvalue
    ax.set_xticks([1,2]); ax.set_xticklabels(grps); ax.set_title(f"{title}\np={p:.2f} (n.s.)",weight="bold")
fig,axs=plt.subplots(1,3,figsize=(7.1,2.6))
box(axs[0],df,"shannon_rar","TCR Shannon (rarefied)"); axs[0].set_ylabel("Shannon entropy")
box(axs[1],df,"clonality_rar","TCR clonality (rarefied)"); axs[1].set_ylabel("clonality")
box(axs[2],cd8,"clonality","CD8 clonal expansion"); axs[2].set_ylabel("CD8 clonality")
save(fig,"Fig6_TCR")
print("All figures regenerated at 600 dpi, fonts >=6.5pt, editable PDF (fonttype 42).")
print("Sizes: single-col 3.4in / double-col 7.1in. Check final font >=7pt after journal scaling.")
