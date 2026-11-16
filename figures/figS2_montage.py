"""Supplementary Figure S2: observed versus matched-depth overlap measures.
2x2 montage, labelled A-D, of the four panels written by pipeline Step R2
(publication_reporting_v1/supplementary_figures/FigureS2A-D_*.png)."""
import glob, os
from _style import *
from _common import load_config
import matplotlib.image as mpimg

cfg = load_config(description=__doc__); OUT = cfg.out_dir("figures")
D = cfg.result("publication_reporting_v1", "supplementary_figures")
panels = [glob.glob(os.path.join(D, f"FigureS2{L}_*.png"))[0] for L in "ABCD"]
fig, axs = plt.subplots(2, 2, figsize=(5.3, 5.3))
fig.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.03, wspace=0.16, hspace=0.08)
for ax, f, L in zip(axs.flat, panels, "ABCD"):
    ax.imshow(mpimg.imread(f)); ax.axis("off")
    ax.text(-0.06, 1.02, L, transform=ax.transAxes, fontsize=12, fontweight="bold", color=INK, va="top", ha="left")
save(fig, OUT, "Figure_S2")
