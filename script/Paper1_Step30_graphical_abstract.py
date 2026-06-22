#!/usr/bin/env python3
"""Step 30: graphical abstract — 70% visual, 30% text.
   Genome Medicine style: icons, arrows, short phrases, minimal text blocks."""
import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
mpl.rcParams.update({"pdf.fonttype":42,"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","DejaVu Sans"]})

# ── Palette ──
B,O,Gy,Gn,Rd,Pu,Dk,WH="#0072B2","#D55E00","#9a9a9a","#009E73","#c0413b","#7A4F9E","#2b2b2b","white"
LtB,LtO,LtPu,LtGn="#d6eaf8","#fdebd0","#e8daef","#d5f5e3"

fig,ax=plt.subplots(figsize=(11.5,5.4)); ax.set_xlim(0,145); ax.set_ylim(0,75); ax.axis("off")

# ── Helpers ──
def box(x,y,w,h,fc,ec,r=2.5,lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.4,rounding_size={r}",
                fc=fc,ec=ec,lw=lw,zorder=2))
def arrow(x1,x2,y,c=Dk,lw=2.0,ms=20):
    ax.add_patch(FancyArrowPatch((x1,y),(x2,y),arrowstyle="-|>",mutation_scale=ms,
                lw=lw,color=c,zorder=1))
def varrow(x,y1,y2,c=Dk,lw=1.5,ms=14):
    ax.add_patch(FancyArrowPatch((x,y1),(x,y2),arrowstyle="-|>",mutation_scale=ms,
                lw=lw,color=c,zorder=1))
def T(x,y,s,fs=8,w="normal",c=Dk,ha="center",va="center",it=False):
    ax.text(x,y,s,fontsize=fs,fontweight=w,color=c,ha=ha,va=va,
            style="italic" if it else "normal",zorder=4)
def circ(x,y,r,fc,ec=Dk,lw=0.8):
    ax.add_patch(Circle((x,y),r,fc=fc,ec=ec,lw=lw,zorder=3))
def rect(x,y,w,h,fc,ec=Dk,lw=0.8):
    ax.add_patch(Rectangle((x,y),w,h,fc=fc,ec=ec,lw=lw,zorder=3))

# ═══════════════════════════════════════════════════════
# HEADER — title + subtitle
# ═══════════════════════════════════════════════════════
box(1,65,143,9,LtO,O,lw=1.6)
T(72,70.5,"A multi-modal audit of Alzheimer's disease blood transcriptomics",
   13,"bold",Dk)
T(72,67.0,"A cell-composition account of the blood-RNA reproducibility paradox",
   9.5,"normal",O,it=True)

# ═══════════════════════════════════════════════════════
# FOUR PANELS
# ═══════════════════════════════════════════════════════
BY,BH=26,36
panels = [
    (1,  LtB, B,  "1"),
    (36, LtO, O,  "2"),
    (71, LtPu,Pu, "3"),
    (106,LtGn,Gn, "4"),
]
for x,fc,ec,num in panels:
    box(x,BY,33,BH,fc,ec)

# arrows between panels
for x1,x2 in [(34.5,35.5),(69.5,70.5),(104.5,105.5)]:
    arrow(x1,x2,BY+BH/2,Dk,1.8,18)

# ───────────────────────────────────────────────────────
# PANEL 1 — PARADOX  (visual: 3 study columns with gene-bars)
# ───────────────────────────────────────────────────────
T(17.5,59.5,"Paradox",11,"bold",B)

# Three "studies" — show as colored boxes with DOTS representing genes (no tiny text)
sx = [7.5, 17.5, 27.5]
colors_study = [B,O,Pu]
np.random.seed(42)
for i,xi in enumerate(sx):
    col = colors_study[i]
    # study label
    T(xi, 57, f"Study {chr(65+i)}", 7.5, "bold", col)
    # mini gene list as colored dots (different patterns per study)
    box(xi-4.5, 44, 9, 11, WH, col, r=1.5, lw=1.0)
    for j in range(5):
        # horizontal "gene bar" - different lengths per study (visual only)
        blen = [6.5,5.0,7.0][i] + np.random.uniform(-1,1)
        rect(xi-blen/2, 51.5-j*1.5, blen, 0.9, col, WH, 0.3)
    # different dot positions to show "different genes"
    for j in range(3):
        circ(xi + np.random.uniform(-2,2), 43 - j*0.1, 0.25, col, WH, 0.2)

T(17.5,41.5,"different genes",7.5,c=Dk,it=True)

# big downward arrow to "same outcome"
varrow(17.5,40.5,37.5,Rd,2.0,18)

T(17.5,35,"similar AD",8.5,"bold",Dk)
T(17.5,32.5,"classification",8.5,"bold",Dk)

# Why? — the question
box(11.5,27,12,4,LtB,Rd,r=1.5,lw=1.0)
T(17.5,29,"Why?",15,"bold",Rd,it=True)

# ───────────────────────────────────────────────────────
# PANEL 2 — AUDIT  (visual: stacked modality icons)
# ───────────────────────────────────────────────────────
T(52.5,59.5,"Audit",11,"bold",O)
T(52.5,56,"~1,350 participants",9.5,"bold",Dk)

# Five modalities as compact colored chips
modalities = [
    ("3 bulk cohorts",    B),
    ("ADNI validation",   O),
    ("scRNA-seq",         Pu),
    ("TCR repertoire",    "#2E86C1"),
    ("Interpretable ML",  Gn),
]
for i,(lab,col) in enumerate(modalities):
    yy = 51.5 - i*3.8
    box(41.5, yy-1.3, 22, 3.0, WH, col, r=1.5, lw=1.0)
    T(52.5, yy+0.2, lab, 8, "bold", col)

T(52.5,31,"5 modalities",8.5,c=O,it=True)
T(52.5,28,"→ same answer",9,"bold",O)

# ───────────────────────────────────────────────────────
# PANEL 3 — CELL-COMPOSITION ACCOUNT  (visual centerpiece)
# ───────────────────────────────────────────────────────
T(87.5,59.5,"Cell-composition",11,"bold",Pu)
T(87.5,56.5,"account",11,"bold",Pu)

# Central equation — LARGE and visual
T(87.5,52,"Blood RNA",13,"bold",Dk)
T(87.5,47.5,"≈",36,"bold",Pu)
T(87.5,43,"Immune-cell",13,"bold",Dk)
T(87.5,40,"composition",13,"bold",Dk)

# AD composition shift — visual cells with arrows
# Myeloid UP (red circles, bigger)
circ(80,34,2.2,"#E74C3C",Dk,0.7); T(80,34,"M",7,"bold",WH)
T(80,30.5,"monocyte ↑",6.5,"bold","#E74C3C")
circ(95,34,2.0,"#E67E22",Dk,0.7); T(95,34,"N",7,"bold",WH)
T(95,30.5,"neutrophil ↑",6.5,"bold","#E67E22")

# Lymphoid DOWN (blue circle, smaller)
circ(87.5,34,1.5,"#3498DB",Dk,0.7); T(87.5,34,"L",6.5,"bold",WH)
T(87.5,37.5,"lymphocyte ↓",6.5,"bold","#3498DB")

# Mini-flow: different panels → same phenotype
T(87.5,26.5,"Different panels",7,c=Dk)
varrow(87.5,26,24.5,Dk,1.2,12)
T(87.5,23,"same phenotype",7.5,"bold",Pu)

# ───────────────────────────────────────────────────────
# PANEL 4 — IMPLICATIONS / RULES
# ───────────────────────────────────────────────────────
T(122.5,59.5,"Implications",11,"bold",Gn)
T(122.5,55.5,"Requirements for",8.5,c=Dk)
T(122.5,53,"credible discovery",8.5,c=Dk)

# Three rules — two-line phrases, centred so they stay inside the panel
rules = [
    ("Symmetric composition", "adjustment"),
    ("Donor-level", "single-cell statistics"),
    ("Composition-aware", "benchmarking"),
]
for i,(l1,l2) in enumerate(rules):
    yy = 47 - i*6.0
    circ(110, yy-0.8, 1.0, Gn, Dk, 0.5)
    T(122.5, yy,    l1, 8.0, "bold", Dk)
    T(122.5, yy-2.4,l2, 7.4, c=Dk)

# ═══════════════════════════════════════════════════════
# BOTTOM — Key numbers (visual, bold)
# ═══════════════════════════════════════════════════════
box(1,2,143,18,"#f8f8f8","#bbb",lw=1.0)
T(72,18,"Key findings",8,"bold",Gy)

nums = [
    (21,  "67%",          "biomarkers lost\nafter de-confounding"),
    (52,  "88%→48%",      "replication collapses\n(symmetric adjustment)"),
    (83,  "873→6",        "scRNA DE genes\n(donor-level stats)"),
    (117, "0.63≈0.63",   "composition vs panel\n(AUC, head-to-head)"),
]
for x,big,sub in nums:
    box(x-9, 6.5, 18, 8.5, WH, "#ddd", r=1.5, lw=0.6)
    T(x, 12, big, 17, "bold", Rd)
    lines = sub.split("\n")
    for li,ln in enumerate(lines):
        T(x, 8.8 - li*1.5, ln, 6, c=Dk)

# ── Save ──
fig.savefig("figures/pub/Graphical_Abstract.pdf",bbox_inches="tight")
fig.savefig("figures/pub/Graphical_Abstract.png",dpi=600,bbox_inches="tight")
plt.close()
print("Saved Graphical_Abstract.{pdf,png}")
