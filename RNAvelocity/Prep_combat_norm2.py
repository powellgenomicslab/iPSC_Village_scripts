import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import anndata as ad
import matplotlib as plt
import numba
import umap
import itertools
from combat.pycombat import pycombat
import math
from scipy.sparse import csr_matrix
import skmisc
from plotnine import *


datadir = '/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/seurat/'
genedir = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap/"


# scvelo settings
scv.logging.print_version()
# Running scvelo 0.2.3 (python 3.8.5) on 2021-05-31 10:45.
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo', figsize=(40,40))  # for beautified visualization
​
# scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
# scanpy==1.7.2 anndata==0.7.6 numpy==1.20.3 scipy==1.6.3 pandas==1.2.4 scikit-learn==0.24.2
sc.settings.set_figure_params(dpi=200)
​

### Load in sample data info
samples = pd.read_table('/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/velocyto_files.tsv', delim_whitespace=True, header=0)

# ## And Genes to use
# gene_df = pd.read_csv(genedir + "velocity_genes.csv")

### Read in data
adata = list()

for i in range(0,8):
	adata.append(scv.read(samples['Directory'][i], cache=True))


### Update droplet names to match seurat names
for i in range(0,8):
	adata[i].obs_names = [row.replace("x", "") for row in [row.replace(":", "_") for row in adata[i].obs_names]]
	adata[i].var_names_make_unique()


#### Combine the data together
adata_concat = ad.concat(adata)


### Get metadata
meta = pd.read_csv(datadir + "metadata.csv")
# ## PHATE
# umap_cord = pd.read_csv(datadir + "cell_embeddings_phate.csv")
# # ## Seurat
# umap_cord = pd.read_csv(datadir + "cell_embeddings.csv")


### Filter anndata object for droplet ids
sample_obs = pd.read_csv(datadir + "cellID_obs.csv")
adata_concat = adata_concat[np.isin(adata_concat.obs.index,sample_obs["x"])]


### Order the metadata droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Cell ID'})

meta = meta.rename(columns = {'Unnamed: 0':'Cell ID'})
meta_ordered = adata_barcodes.merge(meta, on = "Cell ID")


# ### Order the umap droplets in same orderas adata
# umap = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
# umap_ordered = adata_barcodes.merge(umap, on = "Cell ID")


### Update metadata columns
meta_ordered['Site'] = meta_ordered['MULTI_ID'].str.replace('\d','')

## Add cryopreserved to site for necessary syndey samples
meta_ordered['Site'][["Thawed" in i for i in meta_ordered['Time']]] = meta_ordered['Site']+ "_cryopreserved"

meta_ordered['Time'][["Village Day 4" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 7" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 0" in i for i in meta_ordered['Time']]] = "Baseline"

meta_ordered['Site_Individual'] = meta_ordered['Site'] + "-" + meta_ordered['Final_Assignment']



##### Separate adata_concat into location-time conditions for separate processing
adata_sep = list()
# umap_sep = list()
meta_sep = list()
conditions = list()

for cond in meta_ordered['Site_Individual'].unique():
	print(cond)
	conditions.append(cond)
	### Add UMAP info to anndata
	# umap_sep.append(umap_ordered[meta_ordered['Site_Individual'].isin([cond])].iloc[:,1:])
	adata_sep.append(adata_concat[meta_ordered['Site_Individual'].isin([cond])])
	meta_sep.append(meta_ordered[meta_ordered['Site_Individual'].isin([cond])].iloc[:,1:])


for i in range(0,12):
	# adata_sep[i].obsm['X_umap'] = umap_sep[i].values
	adata_sep[i].obs['clusters'] = meta_sep[i]['integrated_snn_res.0.28'].values
	adata_sep[i].obs['individual'] = meta_sep[i]['Final_Assignment'].values
	adata_sep[i].obs['cell_cycle'] = meta_sep[i]['phases'].values
	adata_sep[i].obs['site_rep'] = meta_sep[i]['Site_rep'].values
	adata_sep[i].obs['Location'] = meta_sep[i]['Site'].values
	adata_sep[i].obs['Time'] = meta_sep[i]['Time'].values
	adata_sep[i].obs['Site_Individual'] = meta_sep[i]['Site_Individual'].values



### Filter
adata_sep_norm = list()
for i in range(0,12):
	# Also filter cells with too few unspliced counts
	adata_sep[i].obs['n_unspliced_counts'] = adata_sep[i].layers['unspliced'].sum(axis=1).A1
	adata_sep[i] = adata_sep[i][adata_sep[i].obs['n_unspliced_counts'] >= 1000]
	# Filter genes
	sc.pp.filter_genes(adata_sep[i], min_cells=20)
	# Filter genes that don't have unspliced counts detected in at least 10 cells
	adata_sep[i].var['n_unspliced_counts'] = np.count_nonzero(adata_sep[i].layers['unspliced'].toarray(), axis=0)
	adata_sep[i] = adata_sep[i][:,adata_sep[i].var['n_unspliced_counts']>10]
	scv.pp.filter_and_normalize(adata_sep[i])
	sc.pp.highly_variable_genes(adata_sep[i], min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = "Site_Individual")



pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['matrix'])
pd.DataFrame.sparse.from_spmatrix(adata_sep[i].X)
pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['spliced'])
pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['unspliced'])


Mtmp = list()
Rtmp = list()
M = list()
R = list()
numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']

for i in range(0,12):
	Mtmp.append(pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['spliced'] + adata_sep[i].layers['unspliced']))
	M.append(Mtmp[i].T)
	Rtmp.append(pd.DataFrame(adata_sep[i].layers['spliced']/(adata_sep[i].layers['spliced'] + adata_sep[i].layers['unspliced'])))
	R.append(Rtmp[i].T)
	# for c in [c for c in M[i].columns]:
	# 	M[i][c] = np.log10(M[i][c] +1)



##### Batch Normalize the data #####
Mb = list()
for i in range(0,12):
	Mb.append(pycombat(M[i], adata_sep[i].obs['Time']))



##### Convert back to S and U matrices
Sb = list()
Ub = list()

for i in range(0,12):
	Sb.append(Mb[i] * R[i])
	Sb[i] = Sb[i].replace(np.nan, 0).T
	Ub.append(Mb[i] * (1 - R[i]))
	Ub[i] = Ub[i].replace(np.nan, 0).T


pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['spliced'])
Sb[i]
pd.DataFrame.sparse.from_spmatrix(adata_sep[i].layers['unspliced'])
Ub[i]
pd.DataFrame.sparse.from_spmatrix(adata_sep[i].X)
Mb[i]

##### Add the updated dataframes to the adata
for i in range(0,12):
	adata_sep[i].layers['spliced'] = csr_matrix(Sb[i])
	adata_sep[i].layers['unspliced'] = csr_matrix(Ub[i])
	# adata_sep[i].layers['matrix'] = csr_matrix(Mb[i].T)



### Get abundances of spliced to unspliced
for i in range(0,12):
	scv.utils.show_proportions(adata_sep[i])


for i in range(0,12):
	print(i)
	##### RNA Velocity #####
	scv.pp.moments(adata_sep[i])
	scv.tl.velocity(adata_sep[i], mode = "stochastic")
	scv.tl.velocity_graph(adata_sep[i])
	scv.tl.velocity_pseudotime(adata_sep[i])



# adata_sep[i].layers['Ms'] = csr_matrix(adata_sep[i].layers['Ms'])
# adata_sep[i].layers['Mu'] = csr_matrix(adata_sep[i].layers['Mu'])

 



### Try with scvelo umap
outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_combat_corrected2/"
os.makedirs(outdir)
os.chdir(outdir)

# i = 1

for i in range(0,12):
	name = conditions[i]
	## UMAP + clustering
	# Identify highly-variable genes.
	sc.pp.highly_variable_genes(adata_sep[i], min_mean=0.0125, max_mean=3, min_disp=0.5)
	sc.pl.highly_variable_genes(adata_sep[i], save= name + "_variable_genes.png")
	# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
	sc.pp.scale(adata_sep[i], max_value=1)
	# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
	sc.tl.pca(adata_sep[i], svd_solver='arpack', n_comps = 50)
	# compute the neighborhood graph of cells using the PCA representation of the data matrix. 
	# You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
	sc.pp.neighbors(adata_sep[i], n_neighbors=10)
	# Embedding the neighborhood graph - reccomend embedding the graph in 2 dimensions using UMAP (McInnes et al., 2018).
	#  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
	#  In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
	#tl.paga(adata)
	#pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
	#tl.umap(adata, init_pos='paga')
	sc.tl.umap(adata_sep[i])
	# Clustering the neighborhood graph - recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018).
	# Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
	sc.tl.leiden(adata_sep[i])
	sc.pl.umap(adata_sep[i], color=['leiden'], save= name + "_umap.pdf")
	sc.pl.umap(adata_sep[i], color=['clusters'], save= name + "_umap_seurat_clusters.pdf")
	sc.pl.umap(adata_sep[i], color=['clusters'], save= name + "_umap_seurat_clusters.png")
	sc.pl.umap(adata_sep[i], color=['Time'], save= name + "_umap_time.png")
	scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap', color=['leiden'], save= name + "_stream.png", size=500)
	# Creates a plot AND adds embedded velocity vectors 'velocity_umap' to the obsm slot
	scv.pl.velocity_embedding(adata_sep[i], basis='umap', arrow_length=2, arrow_size=1.5, dpi=150, color=['leiden'], save= name + "_arrow.png")
	# Just a different representation, nothing new is added to adata



top_genes = list()

for i in range(0,12):
	name = conditions[i]
	scv.tl.recover_dynamics(adata_sep[i], fit_basal_transcription=True, max_iter=20)
	scv.tl.velocity(adata_sep[i], mode='dynamical')
	scv.tl.velocity_graph(adata_sep[i])
	scv.tl.recover_latent_time(adata_sep[i])
	scv.pl.scatter(adata_sep[i], color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= name + "_latent_time.png")
	top_genes.append(adata_sep[i].var['fit_likelihood'].sort_values(ascending=False).index[:100])
	# top_genes = adata_sep[i].var['fit_likelihood'].sort_values(ascending=False).index[:100]
	scv.pl.heatmap(adata_sep[i], var_names = top_genes[i][1:100], sortby='latent_time', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= name + "_heatmap.png")
	scv.pl.heatmap(adata_sep[i], var_names = top_genes[i][1:100], sortby='velocity_pseudotime', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= name + "pseudotime_heatmap.png")
	# scv.pl.heatmap(adata_sep[i], var_names = top_genes[1:100], sortby='latent_time', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= name + "_heatmap.png")
	# scv.pl.heatmap(adata_sep[i], var_names = top_genes[1:100], sortby='velocity_pseudotime', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= name + "pseudotime_heatmap.png")
	scv.pl.scatter(adata_sep[i], color='clusters', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= name + "_seurat_clusters.png")
	scv.pl.scatter(adata_sep[i], color='velocity_pseudotime', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= name + "_velocity_pseudotime_umap.png")
	scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap', color=['latent_time'], save= name + "_stream_latent_colors.png", size=150)
	scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap', color=['clusters'], save= name + "_stream_seurat_clusters.png", size=50)
	plot = (ggplot(adata_sep[i].obs, aes(x='latent_time', color='Time', fill='Time')) + geom_density(alpha=0.1))
	plot.save(filename = name + "_density_plot.png")
	adata_sep[i].write(filename=outdir +  name + "_adata_concat.h5ad")
	adata_sep[i].obs.to_csv(outdir + name + "_metadata.csv")





