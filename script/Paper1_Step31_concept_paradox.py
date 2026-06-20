#!/usr/bin/env python3
"""Step 31: 'Reproducibility paradox' concept figure — old vs new model of AD blood-RNA signal."""
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"]})
B,O,Gy,Gn,Rd,Pu,Dk="#0072B2","#D55E00","#999","#009E73","#c0413b","#7A4F9E","#222"

fig,ax=plt.subplots(figsize=(7.4,4.6)); ax.set_xlim(0,100); ax.set_ylim(0,64); ax.axis("off")

def box(x,y,w,h,fc,ec,r=1.8,lw=1.0):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.2,rounding_size={r}",fc=fc,ec=ec,lw=lw,zorder=2))
def arrow(x1,y1,x2,y2,c=Dk,lw=1.4,head=12):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=head,lw=lw,color=c,zorder=1))
def T(x,y,s,fs=8,w="normal",c=Dk,ha="center",va="center",it=False):
    ax.text(x,y,s,fontsize=fs,fontweight=w,color=c,ha=ha,va=va,
            style="italic" if it else "normal",zorder=3)

# title
T(50,60.5,"The reproducibility paradox of AD blood transcriptomics, and its resolution",11,"bold",Dk)

# ---- LEFT: paradox ----
box(2,7,42,48,"#fafafa","#bbb",lw=0.8)
T(23,52.0,"The paradox",9,"bold",Dk)
T(23,49.2,"(observed for >10 years)",6.4,c=Gy,it=True)
# three studies, different genes
study_y=[42,34,26]; gene_sets=[["GeneA","GeneB","GeneC"],["GeneD","GeneE","GeneF"],["GeneG","GeneH","GeneI"]]
study_aucs=["AUC 0.91","AUC 0.88","AUC 0.93"]
for i,(yy,gs,au) in enumerate(zip(study_y,gene_sets,study_aucs)):
    box(5,yy-2.6,16,5.2,"#eef4fb",B,r=1.2,lw=0.6)
    T(13,yy+0.2,f"Study {i+1}",7,"bold",B)
    T(13,yy-1.6,au,5.6,c=B)
    box(25,yy-2.6,16,5.2,"#fff3ec",O,r=1.2,lw=0.6)
    T(33,yy+0.2,", ".join(gs),6,c=Dk)
    T(33,yy-1.6,"reported markers",5.4,c=Gy,it=True)
# arrow study->markers
for yy in study_y:
    arrow(21.2,yy,24.8,yy,c=Gy,lw=0.8,head=8)

# paradox box at bottom
box(5,9.5,36,7,"#fdeceb",Rd,r=1.4)
T(23,14.2,"High accuracy in every study",6.5,"bold",Rd)
T(23,11.5,"...yet gene sets barely overlap, and none replicates clinically.",5.9,c=Rd)

# ---- middle arrow / question ----
T(48.5,32,"Why?",10,"bold",Dk)
arrow(45,32,52,32,c=Dk,lw=1.6,head=14)

# ---- RIGHT: resolution ----
box(55,7,43,48,"#fafafa","#bbb",lw=0.8)
T(76.5,52.0,"The cell-composition account",9,"bold",Gn)
T(76.5,49.2,"(this work)",6.4,c=Gy,it=True)

# AD -> myeloid expansion -> cell counts -> bulk RNA -> any myeloid panel classifies
nodes=[
    (76.5,44.5,"AD",B,7.0),
    (76.5,38.5,"Myeloid-up / lymphoid-down shift\n(reproducible across cohorts)",O,5.8),
    (76.5,30.0,"Bulk blood RNA = composition-weighted\naverage of leukocyte states",Pu,5.8),
    (76.5,21.0,"ANY sufficiently large set of\nmyeloid-marker genes classifies AD\n(by indirectly counting cells)",Gn,5.8),
    (76.5,11.5,"Gene lists differ between studies;\nunderlying signal does NOT.",Dk,6.0),
]
prev=None
for i,(x,y,s,c,fs) in enumerate(nodes):
    fc="#eef4fb" if c==B else "#fff3ec" if c==O else "#f3eef9" if c==Pu else "#eafaf4" if c==Gn else "#f3f3f3"
    h=4.0 if i==0 else 5.8 if "\n" in s else 4.6
    box(60,y-h/2,33,h,fc,c,r=1.2,lw=0.8)
    T(x,y,s,fs,"bold" if i in (0,4) else "normal",c=c)
    if prev is not None:
        arrow(76.5,prev-1.4,76.5,y+(h/2)+0.2,c=Dk,lw=1.0,head=10)
    prev=y-h/2

# bottom takeaway band
box(2,0.5,96,4.3,"#ffeed8","#d4881f",r=1.0,lw=0.8)
T(50,2.7,"Different studies select different gene sets — but they are all reading out the same low-dimensional myeloid composition phenotype.",
  7.4,"bold",Dk)

fig.savefig("figures/pub/FigCM_concept_paradox.pdf",bbox_inches="tight")
fig.savefig("figures/pub/FigCM_concept_paradox.png",dpi=600,bbox_inches="tight")
plt.close()
print("Saved FigCM_concept_paradox.{pdf,png}")
