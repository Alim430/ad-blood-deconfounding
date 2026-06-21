#!/usr/bin/env python3
"""Step 30: graphical abstract — clean four-stage narrative, text constrained inside boxes."""
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"]})
B,O,Gy,Gn,Rd,Dk="#0072B2","#D55E00","#bdbdbd","#009E73","#c0413b","#2b2b2b"

# wider canvas so each stage box is roomy
fig,ax=plt.subplots(figsize=(9.2,4.8)); ax.set_xlim(0,120); ax.set_ylim(0,66); ax.axis("off")

def box(x,y,w,h,fc,ec,r=2.0,lw=1.1):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.2,rounding_size={r}",
                fc=fc,ec=ec,lw=lw,zorder=2))
def arrow(x1,x2,y,c=Dk):
    ax.add_patch(FancyArrowPatch((x1,y),(x2,y),arrowstyle="-|>",mutation_scale=14,lw=1.7,color=c,zorder=1))
def T(x,y,s,fs=8,w="normal",c=Dk,ha="center",va="center",it=False):
    ax.text(x,y,s,fontsize=fs,fontweight=w,color=c,ha=ha,va=va,
            style="italic" if it else "normal",zorder=4)

# ---- title ----
box(3,57.5,114,6.5,"#FFF3EC",O,lw=1.4)
T(60,60.8,"Blood-RNA “biomarkers” of Alzheimer’s disease are myeloid cell counts in disguise",12,"bold",O)

# ---- four stage boxes: centers 16, 46, 76, 106 ; width 26 ; y 34..51 ----
BY,BH=34,17
boxes=[(3,"#eef4fb",B,"The claim"),(33,"#fff7f0",O,"The filter"),
       (63,"#fdeceb",Rd,"The reveal"),(93,"#eefaf4",Gn,"Same null")]
for x,fc,ec,lab in boxes:
    box(x,BY,24,BH,fc,ec)
cx=[15,45,75,105]   # box centres (x+12)

# arrows between boxes (gaps: 27-33, 57-63, 87-93)
for x1,x2 in [(27.5,32.5),(57.5,62.5),(87.5,92.5)]:
    arrow(x1,x2,BY+BH/2)
T(30,45.5,"adjust cell\ncomposition",5.8,c=Dk)
T(60,46.0,"validate:\nadjust BOTH",5.8,c=Dk)
T(90,45.6,"check every\nmodality",5.8,c=Dk)

# Stage 1 — claim
T(15,48.5,"The claim",8.5,"bold",B)
T(15,43.6,"Bulk blood RNA",7,c=Dk)
T(15,40.6,"723 “robust” AD",7,c=Dk)
T(15,37.8,"gene biomarkers",7,c=Dk)
T(15,35.4,"(AUC > 0.9 reported)",5.8,c=Gy,it=True)

# Stage 2 — filter
T(45,48.5,"The filter",8.5,"bold",O)
T(45,43.8,"67% vanish",10,"bold",Rd)
T(45,40.0,"240 of 723 survive",7,c=Dk)
T(45,36.4,"the rest were\njust composition",5.8,c=Gy,it=True)

# Stage 3 — reveal (mini bar fully inside box 63..87)
T(75,48.5,"The reveal",8.5,"bold",Rd)
bx=[71,79]; bv=[88,48]; bc=[B,O]; base=36.5; scale=0.085
ax.plot([67.5,82.5],[base+50*scale,base+50*scale],"--",c=Gy,lw=0.8,zorder=3)  # ~chance
T(83.0,base+50*scale,"chance",5.0,c=Gy,ha="left")
for x,v,c in zip(bx,bv,bc):
    ax.add_patch(plt.Rectangle((x-2.6,base),5.2,v*scale,fc=c,ec="none",zorder=3))
    T(x,base+v*scale+0.9,f"{v}%",6.2,"bold",c)
T(71,35.2,"no comp.",5.0,c=B); T(79,35.2,"+comp.",5.0,c=O)

# Stage 4 — same null (left-aligned bullets inside box 93..117)
T(105,48.5,"Same null,",8,"bold",Gn)
T(105,45.6,"every modality",8,"bold",Gn)
for i,s in enumerate(["bulk RNA","single-cell (873→6)","deep learning (0.52)","TCR repertoire"]):
    T(95,42.4-i*2.4,"• "+s,5.8,c=Dk,ha="left")

# ---- bottom takeaway ----
box(3,3,114,12,"#f3f3f3","#999",lw=1.0)
T(60,11.4,"The AD blood gene signal is predominantly explained by myeloid leukocyte composition.",10,"bold",Dk)
T(60,7.4,"Limited evidence for a composition-independent signature; APOE genotype alone (AUC 0.71) beats the whole panel.",7.2,c=Dk)
T(60,4.3,"Lesson: adjust composition on discovery AND validation, and use donor-level statistics — or you rediscover cell counts.",6.6,c=O,it=True)

fig.savefig("figures/pub/Graphical_Abstract.pdf",bbox_inches="tight")
fig.savefig("figures/pub/Graphical_Abstract.png",dpi=600,bbox_inches="tight"); plt.close()
print("Saved Graphical_Abstract.{pdf,png}")
