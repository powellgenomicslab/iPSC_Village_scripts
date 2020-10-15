import sys
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
 
print('Arguments:', len(sys.argv))
print('List:', str(sys.argv))
var_number = float(sys.argv[3])
print(var_number)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

## Basic run with scrublet
input_dir = os.path.join(sys.argv[1])
counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1)) # Use with the raw data
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15, sim_doublet_ratio = 2)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=150, 
                                                          min_gene_variability_pctl=var_number, 
                                                          n_prin_comps=30)

scrub.call_doublets(threshold=0.40)

outdir = sys.argv[2]
scrub.plot_histogram();
plt.savefig(os.path.join(outdir,'figure1.png'))
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')
scrub.plot_embedding('UMAP', order_points=True);
plt.savefig(os.path.join(outdir,'figure2.png'))

np.savetxt(os.path.join(outdir,'predicted_doublet_mask.txt'), scrub.predicted_doublets_, fmt='%s')
np.savetxt(os.path.join(outdir,'doublet_scores.txt'), scrub.doublet_scores_obs_, fmt='%.4f')

