#!/usr/bin/env python3
"""Step 29: rich multi-panel (CNS-style) main figures from extracted panel data.
Each figure = one 6-panel plate (2x3), publication spec (600dpi, Arial, >=7pt, editable PDF)."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
mpl.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.titlesize":7.5,
    "axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,"legend.fontsize":6,
    "axes.linewidth":0.6,"lines.linewidth":1.0,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600,"figure.dpi":600,
    "savefig.bbox":"tight","savefig.pad_inches":0.03})
B,O,Gy,Gn,Pu="#0072B2","#D55E00","#999999","#009E73","#9467bd"
def save(fig,n): fig.savefig(f"figures/pub/{n}.pdf"); fig.savefig(f"figures/pub/{n}.png"); plt.close(fig)
def lab(ax,L): ax.text(-0.18,1.06,L,transform=ax.transAxes,fontsize=10,fontweight="bold",va="top")
def bx(ax,d,col,gcol,title,ylab):
    grp=["Control","AD"]; data=[d[d[gcol]==g][col].dropna().values for g in grp]
    bp=ax.boxplot(data,widths=.6,patch_artist=True,showfliers=False)
    for p,g in zip(bp['boxes'],grp): p.set_facecolor(B if g=="Control" else O); p.set_alpha(.55); p.set_linewidth(.6)
    for m in bp['medians']: m.set_color("k"); m.set_linewidth(.8)
    for i in range(2): ax.scatter(np.random.normal(i+1,.05,len(data[i])),data[i],s=3,c=B if i==0 else O,alpha=.4,lw=0)
    try: p=mannwhitneyu(data[0],data[1]).pvalue; ax.set_title(f"{title}\np={p:.1e}",weight="bold")
    except: ax.set_title(title,weight="bold")
    ax.set_xticks([1,2]); ax.set_xticklabels(grp); ax.set_ylabel(ylab)

# ============ FIG 1: composition (6 panels) ============
cs=pd.read_csv("results/cell_composition_shift.csv")
cl=pd.read_csv("results/panels/composition_long.csv")
fig,ax=plt.subplots(2,3,figsize=(7.2,4.7))
# a: grouped bar AD-CN
a=ax[0,0]; lin=["Neutrophil","Monocyte","NK","Tcell_CD8","Tcell_CD4","Bcell"]; labs=["Neu","Mono","NK","CD8","CD4","B"]
cohs=cs.cohort.unique(); cols={cohs[0]:B,cohs[1]:O,cohs[2]:Gn}; x=np.arange(len(lin)); w=.26
for i,co in enumerate(cohs):
    d=cs[cs.cohort==co].set_index("lineage"); vals=[float(d.loc[l,"AD_minus_CN"]) if l in d.index else 0 for l in lin]
    a.bar(x+(i-1)*w,vals,w,color=cols[co],label=co.replace("GSE",""))
a.axhline(0,c="k",lw=.5); a.set_xticks(x); a.set_xticklabels(labs,rotation=15); a.set_ylabel("AD − control")
a.set_title("Composition shift (all lineages)",weight="bold"); a.legend(frameon=False,fontsize=5,title="cohort",title_fontsize=5.5); lab(a,"a")
# b-e: per-lineage boxplots pooled
for (axx,lin_name,t,L) in [(ax[0,1],"Monocyte","Monocytes","b"),(ax[0,2],"Neutrophil","Neutrophils","c"),
                            (ax[1,0],"Bcell","B cells","d"),(ax[1,1],"Tcell_CD4","CD4 T cells","e")]:
    bx(axx, cl[cl.lineage==lin_name],"score","group",t,"composition score"); lab(axx,L)
# f: NK
bx(ax[1,2], cl[cl.lineage=="NK"],"score","group","NK cells","composition score"); lab(ax[1,2],"f")
save(fig,"Fig1_composition")

# ============ FIG 2: de-confounding (6 panels) ============
ac=pd.read_csv("results/panels/attenuation_compcor.csv"); i2=pd.read_csv("results/panels/I2.csv")
surv=ac.survives.astype(str).isin(["TRUE","True","true"])
fig,ax=plt.subplots(2,3,figsize=(7.2,4.7))
a=ax[0,0]; a.bar([0,1,2],[723,240,162],color=[Gy,O,B]); a.set_xticks([0,1,2]); a.set_xticklabels(["Tier-1\nnaive","survive\n(no-APOE)","survive\n(+APOE)"])
a.set_ylabel("genes"); a.set_title("De-confounding attrition\n723 → 240 (~67% lost)",weight="bold"); lab(a,"a")
a=ax[0,1]; L=np.nanpercentile(np.abs(np.r_[ac.logFC_orig,ac.logFC_adj]),99.5)
a.plot([-L,L],[-L,L],"k--",lw=.6); a.scatter(ac.logFC_orig[~surv],ac.logFC_adj[~surv],s=4,c=Gy,alpha=.4,lw=0)
a.scatter(ac.logFC_orig[surv],ac.logFC_adj[surv],s=5,c=O,alpha=.7,lw=0); a.set_xlim(-L,L); a.set_ylim(-L,L)
a.set_xlabel("naive log₂FC"); a.set_ylabel("de-confounded log₂FC"); a.set_title("Effect attenuation",weight="bold"); lab(a,"b")
a=ax[0,2]; a.scatter(ac.comp_cor,ac.att,s=4,c=Gy,alpha=.4,lw=0)
z=np.polyfit(ac.comp_cor,ac.att,1); xx=np.linspace(ac.comp_cor.min(),ac.comp_cor.max(),50); a.plot(xx,np.polyval(z,xx),c=O,lw=1.2)
r=np.corrcoef(ac.comp_cor,ac.att)[0,1]; a.set_xlabel("composition correlation"); a.set_ylabel("effect attenuation")
a.set_title(f"Attenuation ∝ composition\nr={r:.2f}",weight="bold"); lab(a,"c")
a=ax[1,0]; a.bar(range(4),[444,173,338,240],color=[Gy,Gy,Gy,O]); a.set_xticks(range(4)); a.set_xticklabels(["−140829","−63060","−63061","full"])
a.set_ylabel("survivors"); a.set_title("Leave-one-out meta\n(not single-cohort)",weight="bold"); lab(a,"d")
a=ax[1,1]; a.hist(i2.I2[~i2.is_surv].dropna(),bins=30,density=True,color=Gy,alpha=.6,label="all genes")
a.hist(i2.I2[i2.is_surv].dropna(),bins=30,density=True,color=O,alpha=.7,label="survivors")
a.set_xlabel("heterogeneity I² (%)"); a.set_ylabel("density"); a.set_title("Survivors are low-I²\n(concordant)",weight="bold"); a.legend(frameon=False); lab(a,"e")
a=ax[1,2]; a.bar(range(3),[1.03,1.08,1.05],color=B); a.axhline(5,ls="--",c=O,lw=.7); a.text(0.1,5.2,"VIF=5 threshold",fontsize=5.5,color=O)
a.set_xticks(range(3)); a.set_xticklabels(["140829","63060","63061"]); a.set_ylim(0,6); a.set_ylabel("VIF (disease term)")
a.set_title("No collinearity\n(VIF≈1)",weight="bold"); lab(a,"f")
save(fig,"Fig2_deconfounding")

# ============ FIG 3: ADNI collapse (6 panels) ============
pg=pd.read_csv("results/panels/adni_pergene.csv"); nn=pd.read_csv("results/panels/negctrl_null.csv")
scv=float(open("results/panels/surv_collapse.txt").read().strip())
fig,ax=plt.subplots(2,3,figsize=(7.2,4.7))
a=ax[0,0]; sets=["240","162","356"]; noc=[90,88,85]; wic=[49,55,46]; x=np.arange(3); w=.38
a.bar(x-w/2,noc,w,color=B,label="no comp."); a.bar(x+w/2,wic,w,color=O,label="+comp."); a.axhline(50,ls="--",c=Gy,lw=.7)
a.set_xticks(x); a.set_xticklabels(sets); a.set_ylim(0,100); a.set_ylabel("concordant (%)"); a.set_title("Replication collapses\nunder composition adj.",weight="bold"); a.legend(frameon=False); lab(a,"a")
a=ax[0,1]; a.bar(range(4),[.076,.025,.029,.025],color=[B,O,Gy,Gy]); a.axhline(.029,ls=":",c=Gy,lw=.7)
a.set_xticks(range(4)); a.set_xticklabels(["surv\nno-c","surv\n+c","all\nno-c","all\n+c"]); a.set_ylabel("median |log₂FC|"); a.set_title("Survivor-specific\ncollapse to background",weight="bold"); lab(a,"b")
a=ax[0,2]; a.plot(range(3),[90,49,55],"o-",color=O,ms=6,lw=1.5); a.axhline(50,ls="--",c=Gy,lw=.7)
a.set_xticks(range(3)); a.set_xticklabels(["no\ncomp","+2\nmyeloid","+6\nlin"]); a.set_ylim(40,100); a.set_ylabel("concordant (%)"); a.set_title("Even 2 covariates\ncollapse it",weight="bold"); lab(a,"c")
a=ax[1,0]; sv=pg.is_surv.astype(bool); LL=np.nanpercentile(np.abs(np.r_[pg.lfc_nocomp,pg.lfc_comp]),99)
a.plot([-LL,LL],[-LL,LL],"k--",lw=.5); a.scatter(pg.lfc_nocomp[sv],pg.lfc_comp[sv],s=5,c=O,alpha=.6,lw=0)
a.set_xlim(-LL,LL); a.set_ylim(-LL,LL); a.set_xlabel("ADNI log₂FC (no comp.)"); a.set_ylabel("ADNI log₂FC (+comp.)"); a.set_title("Survivor effects shrink\nto axis under adj.",weight="bold"); lab(a,"d")
a=ax[1,1]; a.hist(nn.null_collapse,bins=40,color=Gy,alpha=.7); a.axvline(scv,c=O,lw=1.5); a.text(scv*0.95,a.get_ylim()[1]*0.8,f"survivors\n{scv:.3f}",color=O,fontsize=6,ha="right")
a.set_xlabel("effect-size collapse"); a.set_ylabel("random gene sets"); a.set_title("Negative control:\n~40× random (p<0.001)",weight="bold"); lab(a,"e")
a=ax[1,2]; a.hist(ac.comp_cor[~surv].dropna(),bins=30,density=True,color=Gy,alpha=.6,label="all")
a.hist(ac.comp_cor[surv].dropna(),bins=30,density=True,color=O,alpha=.7,label="survivors"); a.set_xlabel("composition correlation"); a.set_ylabel("density")
a.set_title("Survivors are composition-\ncorrelated (0.47 vs 0.16)",weight="bold"); a.legend(frameon=False); lab(a,"f")
save(fig,"Fig3_ADNI_composition_collapse")

# ============ FIG 6: TCR (6 panels) ============
df=pd.read_csv("results/tcr/donor_tcr_metrics.csv"); cd8=pd.read_csv("results/tcr/donor_tcr_CD8.csv")
fig,ax=plt.subplots(2,3,figsize=(7.2,4.7))
bx(ax[0,0],df,"shannon_rar","condition","TCR Shannon (rarefied)","Shannon"); lab(ax[0,0],"a")
bx(ax[0,1],df,"clonality_rar","condition","TCR clonality (rarefied)","clonality"); lab(ax[0,1],"b")
bx(ax[0,2],df,"top10","condition","Top-10 clonotype fraction","fraction"); lab(ax[0,2],"c")
bx(ax[1,0],cd8,"clonality","condition","CD8 clonal expansion","clonality"); lab(ax[1,0],"d")
bx(ax[1,1],cd8,"shannon","condition","CD8 Shannon","Shannon"); lab(ax[1,1],"e")
# f: richness depth artifact
a=ax[1,2]; a.scatter(df.n_cells,df.richness,s=10,c=[O if c=="AD" else B for c in df.condition],lw=0)
a.set_xlabel("T cells captured"); a.set_ylabel("unique clonotypes"); a.set_title("Richness = depth artifact\n(r=0.92)",weight="bold"); lab(a,"f")
save(fig,"Fig6_TCR")
print("Saved multi-panel Fig1, Fig2, Fig3, Fig6 (6 panels each).")
