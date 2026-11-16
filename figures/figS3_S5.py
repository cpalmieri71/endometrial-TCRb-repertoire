"""Supplementary Figures S3-S5: immune panel, genomic panel, association effect sizes.
Requires pipeline Steps 4-6."""
import json, os
from _style import *
from _common import load_config
cfg = load_config(description=__doc__); OUT = cfg.out_dir("figures"); ADIR = cfg.out_dir("analysis")
# ---------------- S3: immune panel ----------------
red=pd.read_csv(cfg.result("04_immune_feature_panel/09_cross_domain_redundancy_pairs.csv"))
stab=pd.read_csv(cfg.result("04_immune_feature_panel/07_domain_selection_stability.csv"))
F=['biopsy_hill_D1','biopsy_top10_frequency','within_top_5pct_pbmc_detection_fraction_std_median','within_morisita_horn_std_median']
L=["Hill D$_1$","Top-10 freq.","Top-5% recovery","Morisita–Horn"]
M=np.full((4,4),np.nan); np.fill_diagonal(M,1.0)
for _,r in red.iterrows():
    if r.feature_a in F and r.feature_b in F:
        i,j=F.index(r.feature_a),F.index(r.feature_b); M[i,j]=M[j,i]=r.spearman_rho
fig,axs=plt.subplots(1,2,figsize=(7.09,3.1),gridspec_kw=dict(width_ratios=[1,0.85],wspace=1.05))
fig.subplots_adjust(left=0.135,right=0.965,top=0.86,bottom=0.34)
im=heat(axs[0],M,L,L,"A  Correlation among frozen immune representatives",
        "No pair reaches the prespecified redundancy threshold of |r| ≥ 0.80.")
cb=fig.colorbar(im,ax=axs[0],fraction=0.046,pad=0.03); cb.set_label("Spearman r",fontsize=6.5)
cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
core=stab[stab.panel_type=='core']
_short={'tumor_effective_diversity':'Effective diversity','tumor_dominance':'Clonal dominance',
 'abundance_conditioned_tumor_to_blood_recovery':'Blood recovery','tumor_blood_abundance_similarity':'Blood similarity'}
dom=[_short.get(d,d.replace('_',' ')) for d in core.domain_id]
axs[1].barh(np.arange(len(core))[::-1],core.most_selected_fraction_all_runs*100,height=0.55,
            color=BLU,zorder=3)
axs[1].set_yticks(np.arange(len(core))[::-1]); axs[1].set_yticklabels(dom,fontsize=6.2)
axs[1].set_xlabel("LOSO runs selecting the same representative (%)")
axs[1].set_xlim(0,105); axs[1].set_title("B  Stability of outcome-blind selection",loc='left',pad=7,fontsize=7.4)
axs[1].grid(True,axis='x',color=GRID,lw=0.5,zorder=0); axs[1].set_axisbelow(True)
for i,v in enumerate(core.most_selected_fraction_all_runs*100):
    axs[1].text(v+2,np.arange(len(core))[::-1][i],f"{v:.0f}%",va='center',fontsize=6,color=INK2)
save(fig, OUT, "Figure_S3")

# ---------------- S4: genomic ----------------
G=pd.read_csv(cfg.result("05_genomic_feature_panel/04_genomic_correlation_spearman.csv"),index_col=0)
GL=["TMB","MSI","GIM","LOH"]
fig,ax=plt.subplots(figsize=(3.6,3.4)); fig.subplots_adjust(left=0.19,right=0.84,top=0.87,bottom=0.33)
im=heat(ax,G.values,GL,GL,"Correlation among frozen genomic measures",
        "Largest pair GIM–LOH (r = 0.769), below the 0.80 redundancy threshold.")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04); cb.set_label("Spearman r",fontsize=6.5)
cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
save(fig, OUT, "Figure_S4")

# ---------------- S5: effect sizes ----------------
pri=pd.read_csv(cfg.result("06_immune_genomic_associations/01_primary_associations.csv"))
sec=pd.read_csv(cfg.result("06_immune_genomic_associations/08_secondary_recovery_sensitivity_associations.csv"))
pol=pd.read_csv(cfg.result("06_immune_genomic_associations/10_exploratory_pole_associations.csv"))
short={'biopsy_hill_D1':'Hill D$_1$','biopsy_top10_frequency':'Top-10 freq.',
 'within_top_5pct_pbmc_detection_fraction_std_median':'Top-5% recovery',
 'within_morisita_horn_std_median':'Morisita–Horn',
 'within_biopsy_clone_fraction_shared_std_median':'Global recovery',
 'within_top_10pct_pbmc_detection_fraction_std_median':'Top-10% recovery',
 'within_biopsy_abundance_fraction_shared_std_median':'Shared read abund.'}
gshort={'oca_tmb':'TMB','oca_msi_score':'MSI','oca_gim':'GIM','oca_loh':'LOH'}
fig=plt.figure(figsize=(7.09,7.0))
gs=fig.add_gridspec(3,1,height_ratios=[16,12,4],hspace=0.34,left=0.30,right=0.97,top=0.965,bottom=0.085)
for ax,d,title in [(fig.add_subplot(gs[0]),pri.sort_values('spearman_rho'),"A  Primary family (16 prespecified tests)"),
                   (fig.add_subplot(gs[1]),sec.sort_values('spearman_rho'),"B  Secondary recovery sensitivity family (12 tests)")]:
    lbl=[f"{short.get(a,a)} × {gshort.get(b,b)}" for a,b in zip(d.immune_feature,d.genomic_feature)]
    forest(ax,lbl,d.spearman_rho.values,d.bootstrap_ci_low.values,d.bootstrap_ci_high.values,
           d.bh_q_within_family.values,"Spearman r (95% bootstrap interval)",title=title)
axp=fig.add_subplot(gs[2])
ec=[c for c in pol.columns if 'rank_biserial' in c][0]
qc=[c for c in pol.columns if 'bh' in c.lower()][0]
locol=[c for c in pol.columns if 'ci_low' in c][0]; hicol=[c for c in pol.columns if 'ci_high' in c][0]
d=pol.sort_values(ec)
forest(axp,[short.get(f,f) for f in d.immune_feature],d[ec].values,d[locol].values,d[hicol].values,
       d[qc].values,"Rank-biserial correlation (95% bootstrap interval)",title="C  Exploratory POLE source flag (4 tests)")
fig.text(0.30,0.018,"Red points mark associations with within-family BH q < 0.05.",fontsize=6,color=INK2)
save(fig, OUT, "Figure_S5")

