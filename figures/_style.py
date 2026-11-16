"""Shared figure style (7.09-in two-column width, 600 dpi) and plotting helpers."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.ticker as mt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import stats
BLU="#4C72B0"; RED="#C44E52"; INK="#1a1a1a"; INK2="#555555"; MUT="#8a8a8a"; GRID="#e4e4e2"; LINE="#c5c5c2"
plt.rcParams.update({"font.family":"FreeSans","font.size":7,"axes.labelsize":7.5,
 "xtick.labelsize":6.5,"ytick.labelsize":6.5,"legend.fontsize":6.5,"axes.edgecolor":"#9a9a98",
 "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2.5,
 "ytick.major.size":2.5,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":600,
 "figure.facecolor":"white","axes.facecolor":"white","pdf.fonttype":42,"ps.fonttype":42})
def lab(ax,L,dx=-0.30,dy=1.18):
    ax.text(dx,dy,L,transform=ax.transAxes,fontsize=9.5,fontweight='bold',color=INK,va='top',ha='left')
def paired(ax,left,right,labels,ylab,logy=True,cl=BLU,cr=RED,medians=True):
    for a,b in zip(left,right): ax.plot([0,1],[a,b],color=LINE,lw=0.8,zorder=2)
    ax.scatter(np.zeros_like(left),left,s=17,color=cl,zorder=4,lw=0.5,edgecolor='white')
    ax.scatter(np.ones_like(right),right,s=17,color=cr,zorder=4,lw=0.5,edgecolor='white')
    if medians:
        for x,v in [(0,left),(1,right)]:
            ax.plot([x-0.19,x+0.19],[np.median(v)]*2,color=INK,lw=1.7,zorder=5,solid_capstyle='butt')
    if logy: ax.set_yscale('log')
    ax.set_xlim(-0.46,1.46); ax.set_xticks([0,1]); ax.set_xticklabels(labels,fontsize=6.7)
    ax.tick_params(axis='x',length=0,pad=3); ax.spines['bottom'].set_visible(False)
    ax.set_ylabel(ylab); ax.grid(True,axis='y',color=GRID,lw=0.5,zorder=0); ax.set_axisbelow(True)
def stat_note(ax,txt,sig=True):
    ax.text(0.5,1.015,txt,transform=ax.transAxes,ha='center',va='bottom',fontsize=6.1,
            color=INK if sig else MUT)
def forest(ax,labels,eff,lo,hi,q,xlabel,title=None,xlim=(-1.05,1.05)):
    y=np.arange(len(labels))[::-1]; sig=np.asarray(q)<0.05
    ax.axvline(0,color=MUT,lw=0.7,zorder=1)
    for i,(l,h) in enumerate(zip(lo,hi)):
        ax.plot([l,h],[y[i],y[i]],color="#b9b9b6",lw=1.4,zorder=2,solid_capstyle='round')
    ax.scatter(eff,y,s=20,color=[RED if s else BLU for s in sig],zorder=4,lw=0.6,edgecolor='white')
    ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=6.2)
    ax.set_xlabel(xlabel); ax.set_xlim(*xlim); ax.set_ylim(-0.7,len(labels)-0.3)
    ax.grid(True,axis='x',color=GRID,lw=0.5,zorder=0); ax.set_axisbelow(True)
    if title: ax.set_title(title,loc='left',pad=7,fontsize=7.4,color=INK)
CMAP=LinearSegmentedColormap.from_list("bwr",[BLU,"#ffffff",RED])
from _common import irls as _irls


def irls(X, y):
    b, _, H = _irls(X, y); return b, H


def heat(ax, M, xl, yl, title, note=None, fs=6.4):
    im = ax.imshow(M, cmap=CMAP, norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha='center', va='center', fontsize=fs, color='white' if abs(M[i, j]) > 0.55 else INK)
    ax.set_xticks(range(len(xl))); ax.set_xticklabels(xl, fontsize=6.6, rotation=28, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl, fontsize=6.6); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, loc='left', pad=7, fontsize=7.4, color=INK)
    if note: ax.text(0, -0.44, note, transform=ax.transAxes, fontsize=6, color=INK2, va='top')
    return im


def save(fig, out, name):
    fig.savefig(os.path.join(out, f"{name}.pdf"), bbox_inches='tight')
    fig.savefig(os.path.join(out, f"{name}.png"), bbox_inches='tight', dpi=600)
    print("wrote", name)
