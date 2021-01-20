import sys
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import read10x
 
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
counts_matrix = read10x.import_cellranger_mtx(input_dir + "matrix.mtx.gz")
barcodes_df = read10x.read_barcodes(input_dir + "barcodes.tsv.gz")


dbl_rate = counts_matrix.shape[0]/1000 * 0.008

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=150, 
                                                          min_gene_variability_pctl=var_number, 
                                                          n_prin_comps=30)



outdir = sys.argv[2]
scrublet_doublet_threshold = sys.argv[4]

try:
    scrublet_doublet_threshold = float(scrublet_doublet_threshold)
except:
    scrublet_doublet_threshold = None
print('scrublet doublet threshold: ')
print(scrublet_doublet_threshold)

np.savetxt(os.path.join(outdir,'predicted_doublet_mask.txt'), scrub.predicted_doublets_, fmt='%s')
np.savetxt(os.path.join(outdir,'doublet_scores.txt'), scrub.doublet_scores_obs_, fmt='%.4f')



if scrublet_doublet_threshold is None:
  ### Plotting and saving
  scrub.plot_histogram();
  plt.savefig(os.path.join(outdir,'doublet_score_histogram.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(outdir,'UMAP.png'))

  results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
  scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
  dataframe = pd.concat([barcodes_df, results, scores], axis=1)
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

  dataframe.to_csv(os.path.join(outdir,'scrublet_results.txt'), sep = "\t", index = False)

else:
  print(scrublet_doublet_threshold)
  scrub.call_doublets(threshold=scrublet_doublet_threshold)
  ### Plotting and saving
  scrub.plot_histogram(); 
  plt.savefig(os.path.join(outdir,'doublet_score_histogram_manual_threshold.png'))
  print('Running UMAP...')
  scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
  print('Done.')
  scrub.plot_embedding('UMAP', order_points=True);
  plt.savefig(os.path.join(outdir,'UMAP_manual_threshold.png'))

  results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
  scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
  dataframe = pd.concat([barcodes_df, results, scores], axis=1)
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
  dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")
  
  dataframe.to_csv(os.path.join(outdir,'scrublet_results.txt'), sep = "\t", index = False)