"""Figures 1 and 2. Requires analysis A1 and A3 outputs and pipeline Step 2."""
import json, os
from _style import *
from _common import load_config
cfg = load_config(description=__doc__); OUT = cfg.out_dir("figures"); ADIR = cfg.out_dir("analysis")
OV = cfg.result("02_exact_clonotype_overlap_500") + "/"
m=pd.read_csv(os.path.join(ADIR, "recomputed_metrics.csv"))
w=m.pivot(index='subject',columns='compartment',values=['observed_richness','hill_D1','hill_D2','top10_frequency'])
res=pd.read_csv(OV+"within_subject_resampling_summary.csv")
st=pd.read_csv(OV+"within_subject_biopsy_abundance_strata_resampling_summary.csv")
st=st[st.metric=='pbmc_detection_fraction']
raw=pd.read_csv(OV+"within_subject_overlap_summary_raw.csv")

# ============ Figure 1: compartmental repertoire architecture ============
fig,axs=plt.subplots(1,4,figsize=(7.09,2.5))
fig.subplots_adjust(left=0.085,right=0.985,top=0.80,bottom=0.17,wspace=0.56)
specs=[('observed_richness',"Observed clonotypes",True,"A"),('hill_D1',"Hill D$_1$",True,"B"),
       ('hill_D2',"Hill D$_2$",True,"C"),('top10_frequency',"Top-10 cumulative read frequency (%)",False,"D")]
for ax,(col,ylab,logy,L) in zip(axs,specs):
    p=w[(col,'PBMC')].values; b=w[(col,'BIOPSY')].values
    if col=='top10_frequency': p,b=p*100,b*100
    paired(ax,p,b,["PBMC","Tumor\nbiopsy"],ylab,logy=logy)
    wt=stats.wilcoxon(b,p); n=int((b<p).sum())
    stat_note(ax,f"{n}/{len(b)} lower in biopsy\nP = {wt.pvalue:.4f}",wt.pvalue<0.05)
    lab(ax,L,dx=-0.44,dy=1.30)
    if col=='top10_frequency': ax.set_ylim(0,52)
    if col=='hill_D2':
        ax.set_yticks([20,50,100,200,500]); ax.yaxis.set_minor_locator(mt.NullLocator())
        ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v,_:f"{v:g}"))
    if col in ('observed_richness','hill_D1'):
        ax.yaxis.set_major_locator(mt.LogLocator(base=10,numticks=5)); ax.yaxis.set_minor_locator(mt.NullLocator())
save(fig, OUT, "Figure_1")

# ============ Figure 2: tumor-blood clonotype sharing ============
cl=pd.read_csv(OV+"within_subject_clonotype_table.csv.gz",low_memory=False)
cl=cl[cl.in_biopsy].copy(); cl['y']=cl.in_pbmc.astype(int); cl['lf']=np.log10(cl.biopsy_frequency)
cl=cl[np.isfinite(cl.lf)]
fig=plt.figure(figsize=(7.09,4.6))
gs=fig.add_gridspec(2,6,hspace=0.72,wspace=1.5,left=0.075,right=0.985,top=0.90,bottom=0.095)
# A: shared clonotypes vs shared read abundance
axA=fig.add_subplot(gs[0,0:2])
cf=res[res.metric=='biopsy_clone_fraction_shared'].set_index('subject').std_median*100
af=res[res.metric=='biopsy_abundance_fraction_shared'].set_index('subject').std_median*100
paired(axA,cf.values,af.values,["Shared\nclonotypes","Shared read\nabundance"],"Fraction of biopsy repertoire (%)",logy=False,cl=BLU,cr=BLU)
axA.set_ylim(0,58); lab(axA,"A",dx=-0.40,dy=1.17)
# B: strata gradient
axB=fig.add_subplot(gs[0,2:4])
order=['top_1pct','1_to_5pct','5_to_10pct','10_to_25pct','25_to_50pct','50_to_100pct']
piv=st.pivot(index='subject',columns='stratum',values='std_median')[order]*100
x=np.arange(6)
for s in piv.index: axB.plot(x,piv.loc[s],color=LINE,lw=0.8,zorder=2)
axB.plot(x,piv.median(),color=INK,lw=1.7,marker='o',ms=4,zorder=5,solid_capstyle='round')
axB.set_xticks(x); axB.set_xticklabels(["Top 1%","1–5%","5–10%","10–25%","25–50%","50–100%"],
                                        rotation=38,ha='right',rotation_mode='anchor',fontsize=6.2)
axB.set_xlabel("Biopsy abundance-rank stratum",labelpad=1)
axB.set_ylabel("Biopsy clonotypes detected\nin PBMC (%)",linespacing=1.3)
axB.set_ylim(-3,103); axB.grid(True,axis='y',color=GRID,lw=0.5,zorder=0); axB.set_axisbelow(True)
lab(axB,"B",dx=-0.36,dy=1.17)
# C: top1% vs least 50%
axC=fig.add_subplot(gs[0,4:6])
hi=piv['top_1pct'].values; lo=piv['50_to_100pct'].values
paired(axC,hi,lo,["Top 1%","Least-abundant\n50%"],"Biopsy clonotypes detected\nin PBMC (%)",logy=False,cl=BLU,cr=BLU)
axC.set_ylim(-3,103)
wt=stats.wilcoxon(hi,lo); stat_note(axC,f"{int((hi>lo).sum())}/{len(hi)} patients · P = {wt.pvalue:.4f}")
lab(axC,"C",dx=-0.40,dy=1.17)
# D: dose-response
axD=fig.add_subplot(gs[1,0:3])
grid=np.linspace(cl.lf.min(),cl.lf.max(),120)
for s,g in cl.groupby('subject'):
    if g.y.nunique()<2 or len(g)<100: continue
    X=np.column_stack([np.ones(len(g)),g.lf.values]); b,_=irls(X,g.y.values.astype(float))
    l1,h1=np.percentile(g.lf,[1,99]); mk=(grid>=l1)&(grid<=h1)
    axD.plot(grid[mk],1/(1+np.exp(-(b[0]+b[1]*grid[mk]))),color=MUT,lw=0.7,alpha=0.75,zorder=2)
subs=sorted(cl.subject.unique())
Dm=np.column_stack([(cl.subject==s).astype(float).values for s in subs])
X=np.column_stack([Dm,cl.lf.values]); bb,_=irls(X,cl.y.values.astype(float))
axD.plot(grid,1/(1+np.exp(-(np.mean(bb[:len(subs)])+bb[-1]*grid))),color=RED,lw=2.0,zorder=4)
axD.set_xlabel("Intratumoral clonotype frequency (log$_{10}$)")
axD.set_ylabel("Probability of PBMC detection")
axD.set_ylim(-0.03,1.03); axD.grid(True,color=GRID,lw=0.5,zorder=0); axD.set_axisbelow(True)
mdl=json.load(open(os.path.join(ADIR, "clonotype_detection_model.json")))  # written by analysis/a3
axD.text(0.97,0.06,f"OR {mdl['or_per_decade']:.2f} per decade\n(95% CI {mdl['ci95'][0]:.2f}–{mdl['ci95'][1]:.2f})",transform=axD.transAxes,ha='right',va='bottom',
         fontsize=6.3,color=INK,linespacing=1.35)
axD.legend(handles=[Line2D([],[],color=MUT,lw=0.9,label="individual patients"),
                    Line2D([],[],color=RED,lw=1.8,label="within-patient estimate")],
           loc='upper left',frameon=False,handlelength=1.8,borderpad=0.15)
lab(axD,"D",dx=-0.165,dy=1.14)
# E: enrichment
axE=fig.add_subplot(gs[1,3:6])
bi=raw.shared_fraction_biopsy_enriched_ge2x.values*100; pb=raw.shared_fraction_pbmc_enriched_ge2x.values*100
paired(axE,bi,pb,["Enriched\nin tumor","Enriched\nin blood"],"Shared clonotypes ≥2-fold\nenriched (%)",logy=False,cl=RED,cr=BLU)
axE.set_ylim(-4,100)
wt=stats.wilcoxon(bi,pb); stat_note(axE,f"{int((bi>pb).sum())}/{len(bi)} patients · P = {wt.pvalue:.4f}")
lab(axE,"E",dx=-0.165,dy=1.14)
save(fig, OUT, "Figure_2")
