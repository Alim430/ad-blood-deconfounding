#!/usr/bin/env python3
"""Step 30: graphical abstract — the whole story as one clear visual narrative."""
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"]})
B,O,Gy,Gn,Rd,Dk="#0072B2","#D55E00","#bdbdbd","#009E73","#c0413b","#2b2b2b"
fig,ax=plt.subplots(figsize=(7.4,4.3)); ax.set_xlim(0,100); ax.set_ylim(0,62); ax.axis("off")

def box(x,y,w,h,fc,ec,r=2.2,lw=1.0,a=1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.2,rounding_size={r}",
                fc=fc,ec=ec,lw=lw,alpha=a,zorder=2))
def arrow(x1,y1,x2,y2,c=Dk):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=13,lw=1.6,color=c,zorder=1))
def T(x,y,s,fs=8,w="normal",c=Dk,ha="center",va="center",it=False):
    ax.text(x,y,s,fontsize=fs,fontweight=w,color=c,ha=ha,va=va,style="italic" if it else "normal",zorder=3)

# ---- title / hook ----
box(2,53,96,7.5,"#FFF3EC",O,lw=1.4)
T(50,56.8,"Blood-RNA “biomarkers” of Alzheimer’s disease are myeloid cell counts in disguise",10.0,"bold",O)

# ---- Stage 1: the claim ----
box(2,33,20,16,"#eef4fb",B)
T(12,46.5,"The claim",8.5,"bold",B)
T(12,42,"Bulk blood RNA:",7.2,c=Dk)
T(12,38.8,"723 “robust” AD",7.2,c=Dk)
T(12,35.8,"gene biomarkers",7.2,c=Dk)
T(12,31.0,"(AUC > 0.9 reported)",6.3,c=Gy,it=True)

# ---- Stage 2: de-confound ----
arrow(22.5,41,27.5,41)
T(25,44.3,"adjust cell\ncomposition",6.0,c=Dk)
box(28,33,20,16,"#fff7f0",O)
T(38,46.5,"The filter",8.5,"bold",O)
T(38,41.8,"67% vanish",9.5,"bold",Rd)
T(38,38.2,"240 survive",7.4,c=Dk)
T(38,34.8,"— the rest were\ncomposition",6.2,c=Gy,it=True)

# ---- Stage 3: the reveal (validation collapse) ----
arrow(48.5,41,53.5,41)
T(51,44.6,"validate —\nadjust BOTH\nsides",5.8,c=Dk)
box(54,30.5,20,18.5,"#fdeceb",Rd)
T(64,46.7,"The reveal",8.5,"bold",Rd)
# mini bar inside
bx=[60,67]; bv=[88,48]; bc=[B,O]
for x,v,c in zip(bx,bv,bc):
    ax.add_patch(plt.Rectangle((x-2.3,33.2),4.6,v*0.10,fc=c,ec="none",zorder=3))
    T(x,33.0+v*0.10+1.1,f"{v}%",6.3,"bold",c)
ax.plot([57.5,71.5],[33.2+5.0,33.2+5.0],"--",c=Gy,lw=0.8,zorder=4)  # chance line
T(72.4,33.2+5.0,"chance",5.3,c=Gy,ha="left")
T(60,31.2,"no\ncomp.",5.0,c=B); T(67,31.2,"+comp.",5.0,c=O)
T(64,29.3,"signal collapses to background",5.8,c=Rd,it=True)

# ---- Stage 4: multi-modal confirmation ----
arrow(74.5,41,79.5,41)
box(80,33,18,16,"#eefaf4",Gn)
T(89,46.5,"Same null,",8.0,"bold",Gn)
T(89,43.8,"every modality",8.0,"bold",Gn)
for i,(s) in enumerate(["bulk RNA: collapses","single-cell: 873 → 6","deep learning: AUC 0.52","TCR: no signal"]):
    T(89,41.2-i*2.5,"•  "+s,5.8,c=Dk,ha="center")

# ---- bottom takeaway bar ----
box(2,3,96,12,"#f3f3f3","#999",lw=1.0)
T(50,11.3,"The AD blood gene signal is predominantly explained by myeloid leukocyte composition.",9.2,"bold",Dk)
T(50,7.0,"Limited evidence for a composition-independent signature — and APOE genotype alone (AUC 0.71) beats the whole transcriptomic panel.",7.6,c=Dk)
T(50,4.2,"Lesson: adjust composition on discovery AND validation, and use donor-level statistics, or you will rediscover cell counts.",6.8,c=O,it=True)

fig.savefig("figures/pub/Graphical_Abstract.pdf",bbox_inches="tight")
fig.savefig("figures/pub/Graphical_Abstract.png",dpi=600,bbox_inches="tight"); plt.close()
print("Saved Graphical_Abstract.{pdf,png}")
