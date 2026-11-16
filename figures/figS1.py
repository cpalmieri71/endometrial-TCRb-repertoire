"""Supplementary Figure S1: sampling depth, coverage and rarefaction. Requires analysis A1 and A2."""
import json, os
from _style import *
from _common import load_config
cfg = load_config(description=__doc__); OUT = cfg.out_dir("figures"); ADIR = cfg.out_dir("analysis")
m=pd.read_csv(os.path.join(ADIR, "recomputed_metrics.csv"))
w=m.pivot(index='subject',columns='compartment',values=['depth','goods_coverage','observed_richness'])
cv=pd.read_csv(os.path.join(ADIR, "rarefaction_curves.csv"))
A=pd.read_csv(os.path.join(ADIR, "paired_rarefied_perpair.csv"))
COMMON=json.load(open(os.path.join(ADIR, "rarefaction_results.json")))["cohort"]["11 pairs"]["common_n"]  # cohort-wide depth, 11 pairs
fig=plt.figure(figsize=(7.09,7.4))
gs=fig.add_gridspec(3,4,height_ratios=[1.0,1.25,1.0],hspace=0.80,wspace=0.60,
                    left=0.085,right=0.985,top=0.955,bottom=0.065)
def lab(ax,L,dx=-0.33,dy=1.20):
    ax.text(dx,dy,L,transform=ax.transAxes,fontsize=9.5,fontweight='bold',color=INK,va='top',ha='left')
def paired_s1(ax,p,b,ylab,logy=True):
    for x,y in zip(p,b): ax.plot([0,1],[x,y],color="#c5c5c2",lw=0.8,zorder=2)
    ax.scatter(np.zeros_like(p),p,s=16,color=BLU,zorder=4,lw=0.5,edgecolor='white')
    ax.scatter(np.ones_like(b),b,s=16,color=RED,zorder=4,lw=0.5,edgecolor='white')
    for x,v in [(0,p),(1,b)]:
        ax.plot([x-0.18,x+0.18],[np.median(v)]*2,color=INK,lw=1.6,zorder=5,solid_capstyle='butt')
    if logy: ax.set_yscale('log')
    ax.set_xlim(-0.45,1.45); ax.set_xticks([0,1]); ax.set_xticklabels(["PBMC","Tumor\nbiopsy"],fontsize=6.6)
    ax.tick_params(axis='x',length=0,pad=3); ax.spines['bottom'].set_visible(False)
    ax.set_ylabel(ylab); ax.grid(True,axis='y',color=GRID,lw=0.5,zorder=0); ax.set_axisbelow(True)
# A,B,C
axA=fig.add_subplot(gs[0,0]); paired_s1(axA,w[('depth','PBMC')].values,w[('depth','BIOPSY')].values,"Sequencing depth (reads)")
lab(axA,"A")
axB=fig.add_subplot(gs[0,1]); paired_s1(axB,w[('goods_coverage','PBMC')].values,w[('goods_coverage','BIOPSY')].values,"Good's coverage",logy=False)
axB.set_ylim(0.80,1.01); lab(axB,"B")
axC=fig.add_subplot(gs[0,2:])
dr=w[('depth','BIOPSY')]/w[('depth','PBMC')]; rr=w[('observed_richness','BIOPSY')]/w[('observed_richness','PBMC')]
axC.axhline(1,color=MUT,lw=0.7,ls=(0,(3,2)),zorder=1); axC.axvline(1,color=MUT,lw=0.7,ls=(0,(3,2)),zorder=1)
axC.scatter(dr,rr,s=26,color=BLU,lw=0.6,edgecolor='white',zorder=4)
axC.set_xscale('log'); axC.set_yscale('log')
axC.set_xlabel("Biopsy / PBMC sequencing-depth ratio"); axC.set_ylabel("Biopsy / PBMC\nobserved-richness ratio",linespacing=1.3)
axC.grid(True,color=GRID,lw=0.5,zorder=0); axC.set_axisbelow(True)
axC.text(0.03,0.06,f"{int((dr>1).sum())} pairs with deeper biopsy;\nrichness ratio < 1 in all of them",
         transform=axC.transAxes,ha='left',va='bottom',fontsize=6.2,color=INK,linespacing=1.3)
lab(axC,"C",dx=-0.145)
# D rarefaction
axD=fig.add_subplot(gs[1,:])
for (sub,comp),g in cv.groupby(['subject','compartment']):
    g=g.sort_values('n_reads')
    axD.plot(g.n_reads,g.richness,lw=0.9,alpha=0.82,color=RED if comp=="BIOPSY" else BLU,
             ls='-' if comp=="BIOPSY" else (0,(4,1.4)),solid_capstyle='round',zorder=3)
axD.axvline(COMMON,color=MUT,lw=0.7,ls=(0,(1.6,1.6)),zorder=2)
axD.annotate(f"common depth\n{COMMON:,} reads",xy=(COMMON,70),xytext=(COMMON*1.45,72),fontsize=6.2,
             color=INK2,ha='left',va='bottom',linespacing=1.25)
axD.set_xscale('log'); axD.set_yscale('log'); axD.set_xlim(90,1.5e6); axD.set_ylim(40,1.6e5)
axD.set_xlabel("Sequencing depth (reads, subsampled)"); axD.set_ylabel("Expected observed clonotypes")
axD.grid(True,color=GRID,lw=0.5,zorder=0); axD.set_axisbelow(True)
axD.legend(handles=[Line2D([],[],color=BLU,lw=1.3,ls=(0,(4,1.4)),label="PBMC"),
                    Line2D([],[],color=RED,lw=1.3,ls='-',label="Tumor biopsy")],
           loc='upper left',frameon=False,handlelength=2.4,borderpad=0.2)
lab(axD,"D",dx=-0.062,dy=1.10)
# E paired standardised
specs=[("S","Observed clonotypes",True),("D1","Hill D$_1$",True),("D2","Hill D$_2$",True),("T10","Top-10 frequency",False)]
for j,(mk,ttl,logy) in enumerate(specs):
    ax=fig.add_subplot(gs[2,j])
    b=A[f"{mk}_B"].values; p=A[f"{mk}_P"].values
    for bi,pi in zip(b,p): ax.plot([0,1],[pi,bi],color="#c5c5c2",lw=0.8,zorder=2)
    ax.scatter(np.zeros_like(p),p,s=15,color=BLU,zorder=4,lw=0.5,edgecolor='white')
    ax.scatter(np.ones_like(b),b,s=15,color=RED,zorder=4,lw=0.5,edgecolor='white')
    if logy: ax.set_yscale('log')
    else: ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6]); ax.set_ylim(0.03,0.63); ax.yaxis.set_minor_locator(mt.NullLocator())
    if mk=="D2":
        ax.set_yticks([20,50,100,200,500]); ax.yaxis.set_minor_locator(mt.NullLocator())
        ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v,_:f"{v:g}"))
    if mk in ("S","D1"):
        ax.yaxis.set_major_locator(mt.LogLocator(base=10,numticks=5)); ax.yaxis.set_minor_locator(mt.NullLocator())
    ax.set_xlim(-0.42,1.42); ax.set_xticks([0,1]); ax.set_xticklabels(["PBMC","Biopsy"],fontsize=6.5)
    ax.tick_params(axis='x',length=0,pad=3); ax.spines['bottom'].set_visible(False)
    ax.grid(True,axis='y',color=GRID,lw=0.5,zorder=0); ax.set_axisbelow(True)
    wt=stats.wilcoxon(b,p); n=int((b<p).sum())
    ax.text(0.0,1.21,ttl,transform=ax.transAxes,ha='left',va='bottom',fontsize=7.2,color=INK)
    ax.text(0.0,1.025,f"{n}/{len(b)} lower in biopsy · P = {wt.pvalue:.4f}",transform=ax.transAxes,ha='left',
            va='bottom',fontsize=6.0,color=INK if wt.pvalue<0.05 else MUT)
    if j==0: lab(ax,"E",dx=-0.52,dy=1.42)
fig.text(0.085,0.004,"Panel E: within each pair the deeper repertoire was subsampled without replacement to the depth of the shallower "
 "member (mean of 100 subsamples). Wilcoxon signed-rank, two-sided, n = 13.",fontsize=5.8,color=INK2)
save(fig, OUT, "Figure_S1")
