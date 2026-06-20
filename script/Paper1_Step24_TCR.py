#!/usr/bin/env python3
"""
Step 24: TCR repertoire analysis (GSE226602) — the adaptive-immune layer.
Donor-level (n=50) AD vs Control clonal expansion + diversity, with rarefaction
(depth-equalised) to remove capture-depth bias. CD8-specific (Gate 2020 Nature).
Inputs : data/tcr/tcr_paired_clonotypes.csv.gz, data/tcr/tcr_contigs.csv.gz,
         results/ai/cells.csv
Outputs: results/tcr/donor_tcr_metrics.csv, figures/pub/tcr/*.png
"""
import csv, gzip, re, numpy as np, pandas as pd, os, collections
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.makedirs("results/tcr", exist_ok=True); os.makedirs("figures/pub/tcr", exist_ok=True)
np.random.seed(0)
OI={"Control":"#0072B2","AD":"#D55E00"}

# donor -> condition (strip barcode suffix + leading G)
cond={}
ct={}
for r in csv.DictReader(open("results/ai/cells.csv")):
    ct[r['barcode']]=r['cell_type']
    m=re.match(r'(.+?)_[ACGTN]+-1$', r['barcode'])
    if m: cond.setdefault(re.sub(r'^G','',m.group(1)), r['condition'])

def metrics(counts):
    c=np.asarray(counts,float); N=c.sum(); p=c/N; rich=len(c)
    H=-(p*np.log(p)).sum()
    return dict(n_cells=int(N), richness=rich, shannon=H,
                clonality=(1-H/np.log(rich)) if rich>1 else np.nan,
                simpson=(p**2).sum(), top10=np.sort(p)[::-1][:10].sum())

def rarefy(counts, Nsub, B=100):
    c=np.asarray(counts,float); p=c/c.sum(); Hs=[]; Cs=[]
    for _ in range(B):
        s=np.random.multinomial(Nsub,p); s=s[s>0]; pp=s/s.sum(); H=-(pp*np.log(pp)).sum()
        Hs.append(H); Cs.append(1-H/np.log(len(s)) if len(s)>1 else np.nan)
    return np.mean(Hs), np.nanmean(Cs)

# ---- overall repertoire per donor (from paired clonotypes) ----
pc=pd.read_csv("data/tcr/tcr_paired_clonotypes.csv.gz")
rows=[]
for donor,g in pc.groupby("id"):
    if donor not in cond: continue
    m=metrics(g["frequency"].values); m.update(donor=donor, condition=cond[donor]); rows.append(m)
df=pd.DataFrame(rows)
Nmin=int(df[df.n_cells>=50].n_cells.min())
print(f"Donors: {len(df)} | rarefaction depth Nmin={Nmin} (donors >=50 cells: {(df.n_cells>=50).sum()})")
keep=df.n_cells>=Nmin
rar=[rarefy(pc[pc.id==d]["frequency"].values, Nmin) for d in df.donor]
df["shannon_rar"]=[x[0] for x in rar]; df["clonality_rar"]=[x[1] for x in rar]

# ---- CD8-specific clonality (Gate 2020) ----
bc2cl={}
for r in csv.DictReader(gzip.open("data/tcr/tcr_contigs.csv.gz","rt")):
    cl=r.get('raw_clonotype_id','')
    if cl and cl not in ('None','NA'): bc2cl["G"+r['barcode']]=cl   # G-prefixed to match cells.csv
cd8=collections.defaultdict(list)
for bc,cl in bc2cl.items():
    if ct.get(bc)=="CD8T":
        d=re.match(r'(.+?)_[ACGTN]+-1$',bc);
        if d: cd8[re.sub(r'^G','',d.group(1))].append(cl)
cd8rows=[]
for donor,cls in cd8.items():
    if donor not in cond or len(cls)<30: continue
    c=pd.Series(cls).value_counts().values
    m=metrics(c); m.update(donor=donor, condition=cond[donor], cd8_cells=len(cls)); cd8rows.append(m)
cd8df=pd.DataFrame(cd8rows)

# ---- stats ----
def test(d,col):
    a=d[d.condition=="AD"][col].dropna(); b=d[d.condition=="Control"][col].dropna()
    if len(a)<3 or len(b)<3: return np.nan,np.nan,np.nan
    p=mannwhitneyu(a,b,alternative="two-sided").pvalue
    return a.median(), b.median(), p
print("\n=== Overall repertoire (rarefied), AD vs Control ===")
for col in ["shannon_rar","clonality_rar","top10","richness"]:
    am,bm,p=test(df,col); print(f"  {col:14s} AD={am:.3f} Ctrl={bm:.3f}  p={p:.3f}")
print(f"\n=== CD8-specific (Gate 2020 axis), n donors={len(cd8df)} ===")
for col in ["clonality","shannon","top10"]:
    am,bm,p=test(cd8df,col); print(f"  CD8 {col:10s} AD={am:.3f} Ctrl={bm:.3f}  p={p:.3f}")

df.to_csv("results/tcr/donor_tcr_metrics.csv",index=False)
cd8df.to_csv("results/tcr/donor_tcr_CD8.csv",index=False)

# ---- figures ----
def box(d,col,title,fn):
    fig,ax=plt.subplots(figsize=(3.2,3.6))
    grps=["Control","AD"]; data=[d[d.condition==g][col].dropna() for g in grps]
    bp=ax.boxplot(data,labels=grps,widths=.6,patch_artist=True,showfliers=False)
    for patch,g in zip(bp['boxes'],grps): patch.set_facecolor(OI[g]); patch.set_alpha(.6)
    for i,g in enumerate(grps): ax.scatter(np.random.normal(i+1,.06,len(data[i])),data[i],s=14,c=OI[g],zorder=3,edgecolors='white',linewidths=.3)
    am,bm,p=test(d,col); ax.set_title(f"{title}\nMann-Whitney p={p:.3f}",fontsize=10,weight="bold")
    ax.set_ylabel(col); [ax.spines[s].set_visible(False) for s in ["top","right"]]
    fig.tight_layout(); fig.savefig(fn,dpi=300); plt.close()
box(df,"clonality_rar","TCR clonality (rarefied)","figures/pub/tcr/tcr_clonality.png")
box(df,"shannon_rar","TCR Shannon diversity (rarefied)","figures/pub/tcr/tcr_shannon.png")
if len(cd8df): box(cd8df,"clonality","CD8 T clonal expansion","figures/pub/tcr/tcr_cd8_clonality.png")
print("\nSaved: results/tcr/*.csv, figures/pub/tcr/*.png")
