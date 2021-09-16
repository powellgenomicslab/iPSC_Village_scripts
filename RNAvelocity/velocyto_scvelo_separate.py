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


### Get metadata
meta = pd.read_csv(datadir + "metadata.csv")
## PHATE
umap_cord = pd.read_csv(datadir + "cell_embeddings_phate.csv")
# ## Seurat
umap_cord = pd.read_csv(datadir + "cell_embeddings.csv")



### Filter anndata object for droplet ids
sample_obs = pd.read_csv(datadir + "cellID_obs.csv")
adata_concat = adata_concat[np.isin(adata_concat.obs.index,sample_obs["x"])]


### Order the metadata droplets in same orderas adata
adata_barcodes = pd.DataFrame(adata_concat.obs.index)
adata_barcodes = adata_barcodes.rename(columns = {0:'Cell ID'})

meta = meta.rename(columns = {'Unnamed: 0':'Cell ID'})
meta_ordered = adata_barcodes.merge(meta, on = "Cell ID")


### Order the umap droplets in same orderas adata
umap = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = adata_barcodes.merge(umap, on = "Cell ID")


### Update metadata columns
meta_ordered['Site'] = meta_ordered['MULTI_ID'].str.replace('\d','')

## Add cryopreserved to site for necessary syndey samples
meta_ordered['Site'][["Thawed" in i for i in meta_ordered['Time']]] = meta_ordered['Site']+ "_cryopreserved"

meta_ordered['Time'][["Village Day 4" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 7" in i for i in meta_ordered['Time']]] = "Village"
meta_ordered['Time'][["Thawed Village Day 0" in i for i in meta_ordered['Time']]] = "Baseline"

meta_ordered['Time_Site_Individual'] = meta_ordered['Time'] + "-" + meta_ordered['Site'] + "-" + meta_ordered['Final_Assignment']


##### Separate adata_concat into location-time conditions for separate processing
adata_sep = list()
umap_sep = list()
meta_sep = list()

for cond in meta_ordered['Time_Site_Individual'].unique():
	print(cond)
	### Add UMAP info to anndata
	umap_sep.append(umap_ordered[meta_ordered['Time_Site_Individual'].isin([cond])].iloc[:,1:])
	adata_sep.append(adata_concat[meta_ordered['Time_Site_Individual'].isin([cond])])
	meta_sep.append(meta_ordered[meta_ordered['Time_Site_Individual'].isin([cond])].iloc[:,1:])


for i in range(0,24):
	adata_sep[i].obsm['X_umap'] = umap_sep[i].values
	adata_sep[i].obs['clusters'] = meta_sep[i]['integrated_snn_res.0.28'].values
	adata_sep[i].obs['individual'] = meta_sep[i]['Final_Assignment'].values
	adata_sep[i].obs['cell_cycle'] = meta_sep[i]['phases'].values
	adata_sep[i].obs['site_rep'] = meta_sep[i]['Site_rep'].values


### Filter
for i in range(0,24):
	# Also filter cells with too few unspliced counts
	adata_sep[i].obs['n_unspliced_counts'] = adata_sep[i].layers['unspliced'].sum(axis=1).A1
	adata_sep[i] = adata_sep[i][adata_sep[i].obs['n_unspliced_counts'] >= 1000]
	# Filter genes
	sc.pp.filter_genes(adata_sep[i], min_cells=20)
	# Filter genes that don't have unspliced counts detected in at least 10 cells
	adata_sep[i].var['n_unspliced_counts'] = np.count_nonzero(adata_sep[i].layers['unspliced'].toarray(), axis=0)
	adata_sep[i] = adata_sep[i][:,adata_sep[i].var['n_unspliced_counts']>10]


### Get abundances of spliced to unspliced
for i in range(0,24):
	scv.utils.show_proportions(adata_sep[i])
# Abundance of ['spliced', 'unspliced']: [0.79 0.21]
# Abundance of ['spliced', 'unspliced']: [0.78 0.22]
# Abundance of ['spliced', 'unspliced']: [0.77 0.23]
# Abundance of ['spliced', 'unspliced']: [0.8 0.2]
# Abundance of ['spliced', 'unspliced']: [0.8 0.2]
# Abundance of ['spliced', 'unspliced']: [0.81 0.19]
# Abundance of ['spliced', 'unspliced']: [0.83 0.17]
# Abundance of ['spliced', 'unspliced']: [0.76 0.24]


for i in range(0,24):
	print(i)
	##### RNA Velocity #####
	scv.pp.filter_and_normalize(adata_sep[i])
	scv.pp.moments(adata_sep[i])
	scv.tl.velocity(adata_sep[i], mode = "stochastic")
	scv.tl.velocity_graph(adata_sep[i])
	scv.tl.velocity_pseudotime(adata_sep[i])





# ### PHATE Figures ###
# outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/phate/"
# os.chdir(outdir)
# for i in range(0,8):
# 	scv.pl.velocity_embedding(adata_sep[i], basis = 'umap',save="velocity_embedding.png")
# 	scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap',legend_fontsize=12, title='', smooth=0.7, min_mass=1, dpi=600, save= str(i) + "_dynamical_stream.png")
# 	scv.pl.velocity_embedding(adata_sep[i], basis='umap', arrow_length=4, arrow_size=4, dpi=600, save= str(i) + "_dynamical_arrow.png")
# 	# scv.pl.proportions(adata_sep[i], save=i + "_proportions.png")
# 	scv.pl.scatter(adata_sep[i], color='velocity_pseudotime', cmap='gnuplot', dpi=600, save= str(i) + "_pseudotime.png")



# ### Seurat Figures ###
# outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/seurat/"
# os.chdir(outdir)
# for i in range(0,8):
# 	scv.pl.velocity_embedding(adata_sep[i], basis = 'umap',save="velocity_embedding.png")
# 	scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap',legend_fontsize=12, title='', smooth=0.7, min_mass=1, dpi=600, save= str(i) + "_dynamical_stream.png")
# 	scv.pl.velocity_embedding(adata_sep[i], basis='umap', arrow_length=4, arrow_size=4, dpi=600, save= str(i) + "_dynamical_arrow.png")
# 	# scv.pl.proportions(adata_sep[i], save=i + "_proportions.png")
# 	scv.pl.scatter(adata_sep[i], color='velocity_pseudotime', cmap='gnuplot', dpi=600, save= str(i) + "_pseudotime.png")



### Try with scvelo umap
outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/scvelo_umap/"
os.chdir(outdir)
for i in range(0,24):
	# ## UMAP + clustering
	# # Identify highly-variable genes.
	# sc.pp.highly_variable_genes(adata_sep[i], min_mean=0.0125, max_mean=3, min_disp=0.5)
	# sc.pl.highly_variable_genes(adata_sep[i], save= str(i) + "_variable_genes.png")
	# # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
	# sc.pp.scale(adata_sep[i], max_value=10)
	# # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
	# sc.tl.pca(adata_sep[i], svd_solver='arpack')
	# # compute the neighborhood graph of cells using the PCA representation of the data matrix. 
	# # You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
	# sc.pp.neighbors(adata_sep[i], n_neighbors=10, n_pcs=40)
	# # Embedding the neighborhood graph - reccomend embedding the graph in 2 dimensions using UMAP (McInnes et al., 2018).
	# #  It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
	# #  In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
	# #tl.paga(adata)
	# #pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
	# #tl.umap(adata, init_pos='paga')
	# sc.tl.umap(adata_sep[i])
	# # Clustering the neighborhood graph - recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018).
	# # Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
	# sc.tl.leiden(adata_sep[i])
	# sc.pl.umap(adata_sep[i], color=['leiden'], save= str(i) + "_umap.pdf")
	# sc.pl.umap(adata_sep[i], color=['clusters'], save= str(i) + "_umap_seurat_clusters.pdf")
	sc.pl.umap(adata_sep[i], color=['clusters'], save= str(i) + "_umap_seurat_clusters.png")
	# scv.pl.velocity_embedding_stream(adata_sep[i], basis='umap', color=['leiden'], save= str(i) + "_stream.png", size=500)
	# # Creates a plot AND adds embedded velocity vectors 'velocity_umap' to the obsm slot
	# scv.pl.velocity_embedding(adata_sep[i], basis='umap', arrow_length=2, arrow_size=1.5, dpi=150, color=['leiden'], save= str(i) + "_arrow.png")
	# # Just a different representation, nothing new is added to adata


# outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/velocyto/seurat/"
# os.chdir(outdir)

for i in range(0,24):
	scv.tl.recover_dynamics(adata_sep[i], fit_basal_transcription=True, max_iter=20)
	scv.tl.velocity(adata_sep[i], mode='dynamical')
	scv.tl.velocity_graph(adata_sep[i])
	scv.tl.recover_latent_time(adata_sep[i])


top_genes = list()

### Plot la tent times
for i in range(0,24):
	scv.pl.scatter(adata_sep[i], color='latent_time', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= str(i) + "_latent_time.png")
	top_genes.append(adata_sep[i].var['fit_likelihood'].sort_values(ascending=False).index[:100])

for i in range(0,24):
	scv.pl.heatmap(adata_sep[i], var_names=top_genes[i], sortby='latent_time', col_color='clusters', yticklabels=True, n_convolve=100, figsize=(40, 20), save= str(i) + "_heatmap.png")



for i in range(0,24):
	scv.pl.scatter(adata_sep[i], color='clusters', fontsize=24, size=100, color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1], save= str(i) + "_seurat_clusters.png")


for i in range(0,24):
	adata_sep[i].write(filename=outdir +  str(i) + "_adata_concat.h5ad")

for i in range(0,24):
	adata_sep[i].obs.to_csv(outdir + str(i) + "_metadata.csv")



##### Get the genes that define the pseudotime in all conditions #####
uniuqe_top_genes = pd.DataFrame(set(list(itertools.chain.from_iterable(top_genes))), columns =['Gene'])

for i in range(0,24):
	cond = meta_ordered['Time_Site_Individual'].unique()[i]
	print(cond)
	uniuqe_top_genes[cond] = 'NA'
		if gene in top_genes[i]:
			uniuqe_top_genes[cond][uniuqe_top_genes.Gene == gene] = 1
		else:
			uniuqe_top_genes[cond][uniuqe_top_genes.Gene == gene] = 0

uniuqe_top_genes["Sum"] = uniuqe_top_genes.iloc[:,1:25].sum(axis=1)
uniuqe_top_genes = uniuqe_top_genes.sort_values(by = "Sum", ascending= False)


top_genes_in_all = list(set.intersection(*map(set, top_genes)))


##### Write the table of all the genes to use for rerunning the velocity analysis #####
uniuqe_top_genes.to_csv(outdir + "velocity_genes.csv")



# The metadata column with the cluster information: integrated_snn_res.0.28


for i in range(0,8):
	any(["POU5F1" in i for i in top_genes[i]])