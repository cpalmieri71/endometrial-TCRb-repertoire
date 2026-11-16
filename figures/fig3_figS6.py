"""Figure 3 and Supplementary Figure S6. Requires pipeline Steps 6 and 8, and the private
patient key and curated annotations (patient codes and POLE-mutant tumors)."""
import json, os
from _style import *
from _common import load_config, patient_codes, curated_annotations
cfg = load_config(description=__doc__); OUT = cfg.out_dir("figures"); ADIR = cfg.out_dir("analysis")
R = cfg.results + "/"
codes = patient_codes(cfg)
POLE_MUT = [k for k, a in curated_annotations(cfg).items() if a["molecular_class"] == "POLEmut"]
assoc=pd.read_csv(R+"06_immune_genomic_associations/01_primary_associations.csv")
mat=pd.read_csv(R+"06_immune_genomic_associations/12_subject_analysis_matrix.csv")
loso=pd.read_csv(R+"06_immune_genomic_associations/06_primary_loso_rho_long.csv")
hist=pd.read_csv(R+"08_histology_associations/17_all_histology_associations.csv")

# ============ Figure 3: immune-genomic associations ============
fig=plt.figure(figsize=(7.09,4.9))
gs=fig.add_gridspec(2,2,hspace=0.70,wspace=0.52,left=0.115,right=0.975,top=0.93,bottom=0.095)
IM=['biopsy_hill_D1','biopsy_top10_frequency','within_top_5pct_pbmc_detection_fraction_std_median','within_morisita_horn_std_median']
GE=['oca_tmb','oca_msi_score','oca_gim','oca_loh']
IML=["Hill D$_1$","Top-10 freq.","Top-5% recovery","Morisita–Horn"]; GEL=["TMB","MSI","GIM","LOH"]
axA=fig.add_subplot(gs[0,0])
Mx=assoc.pivot(index='immune_feature',columns='genomic_feature',values='spearman_rho').loc[IM,GE].values
Q=assoc.pivot(index='immune_feature',columns='genomic_feature',values='bh_q_within_family').loc[IM,GE].values
F=assoc.pivot(index='immune_feature',columns='genomic_feature',values='maxT_fwer_p').loc[IM,GE].values
im=axA.imshow(Mx,cmap=CMAP,norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1),aspect='auto')
for i in range(4):
    for j in range(4):
        mk="*" if (Q[i,j]<0.05 and F[i,j]<0.05) else ""
        axA.text(j,i,f"{Mx[i,j]:.2f}{mk}",ha='center',va='center',fontsize=6.2,
                 color='white' if abs(Mx[i,j])>0.55 else INK)
axA.set_xticks(range(4)); axA.set_xticklabels(GEL,fontsize=6.7)
axA.set_yticks(range(4)); axA.set_yticklabels(IML,fontsize=6.7); axA.tick_params(length=0)
for sp in axA.spines.values(): sp.set_visible(False)
cb=fig.colorbar(im,ax=axA,fraction=0.046,pad=0.03); cb.set_label("Spearman r",fontsize=6.5)
cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
axA.text(0,-0.21,"* supported by BH-FDR and max-r FWER",transform=axA.transAxes,fontsize=6,color=INK2,va='top')
lab(axA,"A",dx=-0.42,dy=1.16)
for ax,col,ylab,L,pos in [(fig.add_subplot(gs[0,1]),'biopsy_hill_D1',"Intratumoral Hill D$_1$","B",None),
                          (fig.add_subplot(gs[1,0]),'biopsy_top10_frequency',"Intratumoral Top-10\ncumulative read frequency","C",None)]:
    x=mat.oca_tmb.values; y=mat[col].values; pole=mat.subject.isin(POLE_MUT).values  # POLE-mutant: pathogenic exonuclease-domain mutation
    ax.scatter(x[~pole],y[~pole],s=26,color=BLU,lw=0.6,edgecolor='white',zorder=4,label="POLE wild-type")
    ax.scatter(x[pole],y[pole],s=26,color=RED,lw=0.6,edgecolor='white',zorder=4,marker='D',label="POLE-mutant")
    ax.set_xscale('log'); ax.set_xticks([4,6,10,20,40,70])
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v,_:f"{v:g}")); ax.xaxis.set_minor_locator(mt.NullLocator())
    ax.set_xlim(3.2,90)
    if col=='biopsy_hill_D1': ax.set_yscale('log'); ax.set_ylim(22,4200)
    r,_=stats.spearmanr(x,y)
    pp=assoc[(assoc.immune_feature==col)&(assoc.genomic_feature=='oca_tmb')].permutation_p_two_sided.iloc[0]
    ax.set_xlabel("TMB (mutations/Mb)"); ax.set_ylabel(ylab,linespacing=1.3)
    ax.grid(True,color=GRID,lw=0.5,zorder=0); ax.set_axisbelow(True)
    ax.text(0.97,0.95,f"r = {r:.3f}\nP = {pp:.5f}",transform=ax.transAxes,ha='right',va='top',
            fontsize=6.3,color=INK,linespacing=1.3)
    if L=="B": ax.legend(loc='lower left',frameon=False,handlelength=1.1,borderpad=0.15,fontsize=6.2)
    lab(ax,L,dx=-0.26,dy=1.16)
axD=fig.add_subplot(gs[1,1])
lo=loso[(loso.immune_feature=='biopsy_hill_D1')&(loso.genomic_feature=='oca_tmb')].sort_values('rho')
full=assoc[(assoc.immune_feature=='biopsy_hill_D1')&(assoc.genomic_feature=='oca_tmb')].spearman_rho.iloc[0]
yy=np.arange(len(lo))
axD.axvline(full,color=RED,lw=1.4,zorder=2)
axD.scatter(lo.rho,yy,s=20,color=BLU,zorder=4,lw=0.6,edgecolor='white')
axD.set_yticks(yy); 
axD.set_yticklabels([codes[int(s)] for s in lo.left_out_subject],fontsize=5.8)
axD.set_xlabel("Spearman r, TMB vs Hill D$_1$"); axD.set_ylabel("Patient left out",fontsize=7)
axD.set_xlim(-0.92,-0.73); axD.grid(True,axis='x',color=GRID,lw=0.5,zorder=0); axD.set_axisbelow(True)
axD.text(0.97,0.97,f"full cohort\nr = {full:.3f}",transform=axD.transAxes,fontsize=6.3,color=RED,
         ha='right',va='top',linespacing=1.3)
lab(axD,"D",dx=-0.26,dy=1.16)
save(fig, OUT, "Figure_3")

# ============ Supplementary Figure S6: histopathological endpoints ============
short={'biopsy_hill_D1':'Hill D$_1$','biopsy_top10_frequency':'Top-10 freq.',
 'within_top_5pct_pbmc_detection_fraction_std_median':'Top-5% recovery',
 'within_morisita_horn_std_median':'Morisita–Horn','oca_tmb':'TMB','oca_msi_score':'MSI',
 'oca_gim':'GIM','oca_loh':'LOH'}
fams=[('grade','A','Endometrioid grade (n = 11)',"Spearman r (95% bootstrap interval)"),
      ('myometrial_invasion','B','Myometrial invasion >50% (n = 12)',"Rank-biserial r (95% bootstrap interval)"),
      ('lvsi','C','LVSI positive (n = 11)',"Rank-biserial r (95% bootstrap interval)")]
fig,axs=plt.subplots(1,3,figsize=(7.09,2.9))
fig.subplots_adjust(left=0.105,right=0.985,top=0.82,bottom=0.21,wspace=0.80)
for ax,(fam,L,ttl,xl) in zip(axs,fams):
    d=hist[hist.analysis_family==fam].sort_values('effect_size')
    forest(ax,[short.get(f,f) for f in d.descriptor_feature],d.effect_size.values,
           d.bootstrap_ci_low.values,d.bootstrap_ci_high.values,d.family_bh_q.values,xl)
    ax.set_title(ttl,loc='left',pad=8,fontsize=7.2,color=INK)
    ax.set_xlabel(xl,fontsize=6.6)
    lab(ax,L,dx=-0.56,dy=1.26)
fig.text(0.105,0.015,"No comparison reached within-family Benjamini–Hochberg q < 0.05, and none reached the "
 "permutation max-|effect| family-wise threshold.",fontsize=6,color=INK2)
save(fig, OUT, "Figure_S6")
