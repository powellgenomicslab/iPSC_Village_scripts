import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import os
import anndata as ad
import matplotlib as plt
import numba
import umap


datadir = '/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/seurat/'


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


### Read in seurat info
# sample_obs = pd.read_csv(datadir + "cellID_obs_phate.csv")
# umap_cord = pd.read_csv(datadir + "cell_embeddings_phate.csv")

## Seurat
sample_obs = pd.read_csv(datadir + "cellID_obs.csv")
umap_cord = pd.read_csv(datadir + "cell_embeddings.csv")


### Filter anndata object for droplet ids
adata_concat = adata_concat[np.isin(adata_concat.obs.index,sample_obs["x"])]


### Order the umap droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Cell ID'})

umap = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})


umap_ordered = adata_barcodes.merge(umap, on = "Cell ID")

### Add UMAP info to anndata
umap_ordered = umap_ordered.iloc[:,1:]
adata_concat.obsm['X_umap'] = umap_ordered.values




### Filter
# Also filter cells with too few unspliced counts
adata_concat.obs['n_unspliced_counts'] = adata_concat.layers['unspliced'].sum(axis=1).A1
adata_concat = adata_concat[adata_concat.obs['n_unspliced_counts'] >= 1000]

# Filter genes
sc.pp.filter_genes(adata_concat, min_cells=20)

# Filter genes that don't have unspliced counts detected in at least 10 cells
adata_concat.var['n_unspliced_counts'] = np.count_nonzero(adata_concat.layers['unspliced'].toarray(), axis=0)
adata_concat = adata_concat[:,adata_concat.var['n_unspliced_counts']>10]




##### RNA Velocity #####
scv.pp.filter_and_normalize(adata_concat)
scv.pp.moments(adata_concat)
scv.tl.velocity(adata_concat, mode = "stochastic")
scv.tl.velocity_graph(adata_concat)
scv.tl.velocity_pseudotime(adata_concat)



### PHATE Figures ###
outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/phate/"
os.chdir(outdir)
scv.pl.velocity_embedding(adata_concat, basis = 'umap',save="velocity_embedding.png")
scv.pl.velocity_embedding_stream(adata_concat, basis='umap',legend_fontsize=12, title='', smooth=0.7, min_mass=1, dpi=600, save="dynamical_stream.png")
scv.pl.velocity_embedding(adata_concat, basis='umap', arrow_length=4, arrow_size=4, dpi=600, save="dynamical_arrow.png")
scv.pl.proportions(adata_concat, save="proportions.png")
scv.pl.scatter(adata_concat, color='velocity_pseudotime', cmap='gnuplot', dpi=600, save="pseudotime.png")



### Seurat Figures ###
outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/seurat/"
os.chdir(outdir)
scv.pl.velocity_embedding(adata_concat, basis = 'umap',save="velocity_embedding.png")
scv.pl.velocity_embedding_stream(adata_concat, basis='umap',legend_fontsize=12, title='', smooth=0.7, min_mass=1, dpi=600, save="dynamical_stream.png")
scv.pl.velocity_embedding(adata_concat, basis='umap', arrow_length=4, arrow_size=4, dpi=600, save="dynamical_arrow.png")
scv.pl.scatter(adata_concat, color='velocity_pseudotime', cmap='gnuplot', dpi=600, rescale_color=[0,1], save="pseudotime.png")




# ### Try with scvelo umap
# outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap/"
# os.chdir(outdir)

# ## UMAP + clustering
# # Identify highly-variable genes.
# sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc.pl.highly_variable_genes(adata_concat, save= "variable_genes.png")
# # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
# sc.pp.scale(adata_concat, max_value=10)
# # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
# sc.tl.pca(adata_concat, svd_solver='arpack')
# # compute the neighborhood graph of cells using the PCA representation of the data matrix. 
# # You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
# sc.pp.neighbors(adata_concat, n_neighbors=10, n_pcs=40)
# # Embedding the neighborhood graph - reccomend embedding the graph in 2 dimensions using UMAP (McInnes et al., 2018).
# #  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
# #  In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
# #tl.paga(adata)
# #pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
# #tl.umap(adata, init_pos='paga')
# sc.tl.umap(adata_concat)
# # Clustering the neighborhood graph - recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018).
# # Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
# sc.tl.leiden(adata_concat)
# sc.pl.umap(adata_concat, color=['leiden'], save= "umap.pdf")
# scv.pl.velocity_embedding_stream(adata_concat, basis='umap', color=['leiden'], save= "stream.png", size=500)
# # Creates a plot AND adds embedded velocity vectors 'velocity_umap' to the obsm slot
# scv.pl.velocity_embedding(adata_concat, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150, color=['leiden'], save= "arrow.png")
# # Just a different representation, nothing new is added to adata


outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/seurat/"
os.chdir(outdir)

### Dynamical model to get latent time
scv.tl.recover_dynamics(adata_concat, max_iter=20)
scv.tl.velocity(adata_concat, mode='dynamical')
scv.tl.velocity_graph(adata_concat)
scv.tl.recover_latent_time(adata_concat)



adata_concat.write(filename=outdir + "adata_concat.h5ad")
adata_concat = ad.read_h5ad(filename=outdir + "adata_concat.h5ad")





scv.pl.scatter(adata_concat, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], dpi = 600, save="latent_time.png")
scv.pl.scatter(adata_concat, color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], dpi = 600, save="latent_time_norescale.png")

top_genes = adata_concat.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_concat, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100, save="heatmap.png")






