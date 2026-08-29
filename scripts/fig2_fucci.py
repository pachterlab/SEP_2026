# Figure 2
# Reproduces panels A-B: FUCCI top-100 clrcycle analysis and embedding benchmark.
# Inputs: data/raw/GSE121265/{GSE121265_fucci-counts.txt.gz,GSE121265_eset-final.rds.gz,fit.quant.rds}
# Outputs: results/fig2/*.csv and figures/Figure2.{pdf,svg,png}
from pathlib import Path
import sys
import argparse
import subprocess
from io import StringIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"scripts"))
from lib.fucci import TWO_PI,extract_source_tables,load_inputs,run_clrcycle,run_embedding_benchmark

RAW=ROOT/"data/raw/GSE121265"; PROCESSED=ROOT/"data/processed"; OUT=ROOT/"results/fig2"
required=[RAW/"GSE121265_fucci-counts.txt.gz",RAW/"GSE121265_eset-final.rds.gz",RAW/"fit.quant.rds"]
missing=[str(p) for p in required if not p.exists()]
if missing: raise FileNotFoundError("Missing source inputs:\n"+"\n".join(missing))
metadata,features,ranking=extract_source_tables(required[1],required[2],PROCESSED)
expression,panel,phenotype=load_inputs(required[0],metadata,features,ranking)
analysis=run_clrcycle(expression,panel,phenotype,OUT)
embeddings,benchmark=run_embedding_benchmark(analysis,OUT)

# Plotting is opt-in so numerical validation can gate rendering.
parser=argparse.ArgumentParser(); parser.add_argument("--plot",action="store_true")
args=parser.parse_args()
if not args.plot:
 print(benchmark.to_string(index=False)); raise SystemExit(0)
b=np.array([1.045301,2.680390,5.271979])
cmap=LinearSegmentedColormap.from_list("fucci_cycle",[(0,"#18A83B"),(b[0]/TWO_PI-.015,"#18A83B"),(b[0]/TWO_PI+.015,"#FFE119"),(b[1]/TWO_PI-.015,"#FFE119"),(b[1]/TWO_PI+.015,"#E31A1C"),(b[2]/TWO_PI-.012,"#E31A1C"),(b[2]/TWO_PI,"#FFE119"),(b[2]/TWO_PI+.012,"#18A83B"),(1,"#18A83B")],N=256)
theta=analysis["theta"]; cells=analysis["cells"]; order=analysis["order"]

def published_peak_phases():
 code='''a<-commandArgs(TRUE);fit<-readRDS(a[1]);m<-read.delim(a[2]);ts<-sort(m$theta);cat("ENSG\\ttheta_g\\n");for(g in a[-c(1,2)]){v<-fit[[g]]$trend.yy;cat(g,"\\t",if(length(v)==length(ts)) ts[which.max(v)] else NA,"\\n",sep="")}'''
 done=subprocess.run(["Rscript","-e",code,str(RAW/"fit.quant.rds"),str(PROCESSED/"GSE121265_metadata.tsv"),*order.feature.astype(str)],capture_output=True,text=True,check=True)
 return pd.read_csv(StringIO(done.stdout),sep="\t").theta_g.to_numpy(float)

BASE_SIZE=12; TITLE_SIZE=13; PANEL_SIZE=16; GENE_SIZE=10.5; METRIC_SIZE=11
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":BASE_SIZE,
 "axes.titlesize":TITLE_SIZE,"axes.titleweight":"bold","axes.labelsize":BASE_SIZE,
 "xtick.labelsize":BASE_SIZE,"ytick.labelsize":BASE_SIZE,
 "legend.fontsize":BASE_SIZE,"svg.fonttype":"none","pdf.fonttype":42})

# A/B occupy the tall upper-left block; B is twice the width of A.  The three
# expression examples form a deliberately shallow, shared-y row underneath.
fig=plt.figure(figsize=(18,10.5),constrained_layout=False)
outer=fig.add_gridspec(1,2,width_ratios=(2.82,1.28),wspace=.03,
                       left=.055,right=.975,top=.88,bottom=.19)
fig.suptitle(
 "FUCCI cell cycle projection, learned feature geometry (GSE121265), and cell embedding benchmark",
 fontsize=PANEL_SIZE,fontweight="bold",x=.055,ha="left",y=.975,
)
left=outer[0].subgridspec(2,1,height_ratios=(1.0,.42),hspace=.36)
top=left[0].subgridspec(1,3,width_ratios=(1,2,0.06),wspace=.34)
gene_ax=fig.add_subplot(top[0,0]); cell_ax=fig.add_subplot(top[0,1])
examples=left[1].subgridspec(1,4,width_ratios=(1,1,1,.40),wspace=.10)
example_axes=[]
for i in range(3):
 example_axes.append(fig.add_subplot(examples[0,i],sharey=example_axes[0] if i else None))

def panel_label(ax,label,x=-.10,y=1.04):
 ax.text(x,y,label,transform=ax.transAxes,fontsize=PANEL_SIZE,fontweight="bold",
         ha="right",va="bottom",clip_on=False)

cell_ax.scatter(cells.clrcycle_cosine_aligned,cells.clrcycle_sine_aligned,
 c=theta,cmap=cmap,norm=Normalize(0,TWO_PI),s=11,alpha=.82,linewidths=0,rasterized=True)
_clrcycle_row=benchmark.set_index("method").loc["CLRcycle first Fourier"]
_rho=float(analysis["metrics"]["rho_clrcycle"])
_b_metrics=(rf"$\boldsymbol{{\rho}}_{{\mathbf{{clrcycle}}}}$={_rho:.3f}, "
            f"R={_clrcycle_row.circular_agreement:.3f}, "
            f"distance stress={_clrcycle_row.distance_stress:.3f}")
median=float(np.median(cells.clrcycle_radius))
cell_ax.add_patch(plt.Circle((0,0),median,fill=False,ls="--",lw=1,color=".45",
                            label="Median cell radius"))
cell_ax.axhline(0,color=".85",lw=.6); cell_ax.axvline(0,color=".85",lw=.6)
cell_ax.set_aspect("equal",adjustable="datalim"); cell_ax.set_box_aspect(1)
cell_ax.set(xlabel="clrcycle cosine coordinate",ylabel="clrcycle sine coordinate")
cell_ax.set_title("FUCCI cells projected by clrcycle",y=1.10,pad=0)
cell_ax.text(.5,1.025,_b_metrics,transform=cell_ax.transAxes,ha="center",va="bottom",
             fontsize=METRIC_SIZE,fontweight="bold")
cell_ax.legend(loc="upper left",fontsize=8.5,frameon=True,framealpha=.88,
               facecolor="white",edgecolor=".8",handlelength=2.0)
panel_label(cell_ax,"B",x=-.08,y=1.08)

ga=order.clrcycle_angle_aligned.to_numpy(float); peaks=published_peak_phases()
gene_ax.scatter(np.cos(ga),np.sin(ga),c=peaks,cmap=cmap,norm=Normalize(0,TWO_PI),
 s=22,linewidths=0,zorder=3)
gene_ax.add_patch(plt.Circle((0,0),1,fill=False,color=".65",lw=1))
gene_ax.set(xlim=(-1.52,1.52),ylim=(-1.48,1.65)); gene_ax.set_aspect("equal"); gene_ax.axis("off")
# Keep the entire panel-A header above the radial gene labels.
gene_ax.set_title("Learnt feature circle",y=1.14,pad=0)
panel_label(gene_ax,"A",x=-.10,y=1.12)
# Fewer, substantially larger labels are staggered radially to remain legible.
_label_table=order.loc[(order["rank"]<=15)|order.gene.isin({"DIABLO","CDC6","TNPO1"})]
_label_rows=list(_label_table.itertuples())
_label_rows.sort(key=lambda row: float(row.clrcycle_angle_aligned)%TWO_PI)
for _idx,row in enumerate(_label_rows):
 a=float(row.clrcycle_angle_aligned); deg=np.degrees(a)%360
 radius=1.10+.105*(_idx%3); flip=90<deg<270
 gene_ax.text(radius*np.cos(a),radius*np.sin(a),row.gene,fontsize=GENE_SIZE,
  fontweight="bold" if row.gene in {"DIABLO","CDC6","TNPO1"} else "normal",
  rotation=deg+(180 if flip else 0),ha="right" if flip else "left",va="center",
  rotation_mode="anchor",clip_on=False)

grid=np.linspace(0,TWO_PI,400); design=np.c_[np.ones(len(theta)),np.cos(theta),np.sin(theta)]; gd=np.c_[np.ones(len(grid)),np.cos(grid),np.sin(grid)]
_example_values=[]
for _idx,(ax,gene) in enumerate(zip(example_axes,("DIABLO","CDC6","TNPO1"))):
 feature=order.loc[order.gene.eq(gene),"feature"].iloc[0]; j=list(expression.columns).index(feature); y=analysis["centered_clr"][:,j]; _example_values.append(y)
 coef=np.linalg.lstsq(design,y,rcond=None)[0]; peak=np.mod(np.arctan2(coef[2],coef[1]),TWO_PI)
 ax.scatter(theta,y,c=theta,cmap=cmap,norm=Normalize(0,TWO_PI),s=5,alpha=.35,linewidths=0,rasterized=True)
 ax.plot(grid,gd@coef,color="black",lw=1.4); ax.axvline(peak,color=".4",ls="--",lw=.8,
                                                       label="Fitted peak")
 ax.set(xlim=(0,TWO_PI),xticks=[0,np.pi/2,np.pi,3*np.pi/2,TWO_PI],
        xticklabels=["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"],title=gene,
        xlabel=r"FUCCI phase $\theta$")
 if _idx==0: ax.set_ylabel("Centered CLR expression")
 else: ax.tick_params(axis="y",labelleft=False)
 ax.legend(loc="lower left",fontsize=8.5,frameon=True,framealpha=.88,
           facecolor="white",edgecolor=".8",handlelength=2.0)
 panel_label(ax,chr(ord("C")+_idx),x=-.03,y=1.04)
_ymin=min(np.nanmin(y) for y in _example_values); _ymax=max(np.nanmax(y) for y in _example_values)
_ypad=.04*(_ymax-_ymin); example_axes[0].set_ylim(_ymin-_ypad,_ymax+_ypad)

right=outer[1].subgridspec(3,1,hspace=.62); metric=benchmark.set_index("method")
method_order=["PCA (2 PCs)","UMAP (2D)","Spectral embedding (2D)"]
for k,name in enumerate(method_order):
 ax=fig.add_subplot(right[k,0]); xy=embeddings[name]
 ax.scatter(xy[:,0],xy[:,1],c=theta,cmap=cmap,norm=Normalize(0,TWO_PI),s=11,alpha=.8,linewidths=0,rasterized=True)
 row=metric.loc[name]; display_name=name.replace(" (2D)","").replace(" (2 PCs)","")
 ax.set_title(display_name,y=1.12,pad=0)
 ax.text(.5,1.035,f"R={row.circular_agreement:.3f}, distance stress={row.distance_stress:.3f}",
         transform=ax.transAxes,ha="center",va="bottom",fontsize=METRIC_SIZE,
         fontweight="bold")
 ax.set_aspect("equal",adjustable="datalim"); ax.set_box_aspect(1); ax.set_anchor("W")
 ax.set_xlabel("PC1" if name.startswith("PCA") else "UMAP1" if name.startswith("UMAP") else "Spectral 1")
 ax.set_ylabel("PC2" if name.startswith("PCA") else "UMAP2" if name.startswith("UMAP") else "Spectral 2")
 panel_label(ax,chr(ord("F")+k),x=-.16,y=1.17)

# Shared continuous legend with the biological labels used in the source figure.
sm=plt.cm.ScalarMappable(norm=Normalize(0,TWO_PI),cmap=cmap)
cb_ax=fig.add_axes([.065,.075,.74,.025]); cb=fig.colorbar(sm,cax=cb_ax,orientation="horizontal")
cb.set_ticks([0,np.pi/2,np.pi,3*np.pi/2,TWO_PI],labels=["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
bio=cb.ax.secondary_xaxis("top")
bio.set_xticks([.42,1.86,3.98,5.80],labels=["G2/M (GMNN high)","S","G1 (CDT1 high)","G2/M"])
bio.tick_params(length=0,pad=5,labelsize=BASE_SIZE)
cb.set_label(r"FUCCI phase: measured cell $\theta$ and gene $\theta_g$ (radians, modulo $2\pi$)",labelpad=6)
for x in b: cb.ax.axvline(x,color="black",lw=.7,alpha=.55)
for ext in ("pdf","svg","png"):
 fig.savefig(ROOT/"figures"/f"Figure2.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
plt.close(fig); print(benchmark.to_string(index=False))
