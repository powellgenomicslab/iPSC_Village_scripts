# Setup
import importlib
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sccoda.datasets as scd



### Set up directories ###
dir = '/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Expression_Boxplots/pluri_degs/'
outdir = '/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scCODA/'


# Load data
cell_counts = pd.read_csv(outdir + 'cell_line_numbers.tsv', sep = '\t')

print(cell_counts)


# Convert data to anndata object
data_all = dat.from_pandas(cell_counts, covariate_columns=["Group"])

print(data_all)


# Extract condition from mouse name and add it as an extra column to the covariates
data_all.obs["Location"] = data_all.obs["Group"].str.replace(r"_.+", "", regex=True)
data_all.obs["Village"] = data_all.obs["Group"].str.replace(r"Brisbane_", "", regex=True).str.replace(r"Melbourne_", "", regex=True).str.replace(r"Sydney_", "", regex=True).str.replace(r"_Replicate[0-9]_.+", "", regex=True)
data_all.obs["Replicate"] = data_all.obs["Group"].str.replace(r".+_Replicate", "Replicate", regex=True).str.replace("_.+", "", regex=True)
data_all.obs["Cryopreservation"] = data_all.obs["Group"].str.replace(r".+_", "", regex=True)


# Uni-culture vs village
## Brisbane
data_brisbane = data_all[data_all.obs["Location"].isin(["Brisbane"])]
viz.boxplots(data_brisbane, feature_name="Village")
plt.savefig(outdir + 'brisbane_numbers.png')

model_brisbane = mod.CompositionalAnalysis(data_brisbane, formula="Village", reference_cell_type="automatic")
brisbane_sim_results = model_brisbane.sample_hmc()

brisbane_sim_results.summary()
print(brisbane_sim_results.credible_effects())


# Run scCODA with each cell type as the reference
cell_types_brisbane = data_brisbane.var.index
results_cycle_brisbane = pd.DataFrame(index=cell_types_brisbane, columns=["times_credible"]).fillna(0)

models_brisbane = []
results_brisbane = []

for ct in cell_types_brisbane:
    print(f"Reference: {ct}")
    # Run inference
    model_temp = mod.CompositionalAnalysis(data_brisbane, formula="Village", reference_cell_type=ct)
    models_brisbane.append(model_temp)
    results_temp = model_temp.sample_hmc(num_results=20000)
    results_brisbane.append(results_temp)
    # Select credible effects
    cred_eff = results_temp.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)
    # add up credible effects
    results_cycle_brisbane["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle_brisbane["pct_credible"] = results_cycle_brisbane["times_credible"]/len(cell_types_brisbane)
results_cycle_brisbane["is_credible"] = results_cycle_brisbane["pct_credible"] > 0.5
print(results_cycle_brisbane)

# Get average difference
### couldn't work out how to pull from summary in python, will do manually - angry emoji
results_brisbane[0].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           0.000000      1735.548669          0.855849
#                    MBE1006          -0.665075      1502.677050         -0.103652
#                    TOB0421          -1.146676       899.107614         -0.798455

results_brisbane[1].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           0.672738      1738.916835          0.863100
#                    MBE1006           0.000000      1501.590528         -0.107455
#                    TOB0421          -0.481418       896.825971         -0.801995

results_brisbane[2].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           1.155125      1739.197970          0.864988
#                    MBE1006           0.480325      1500.238853         -0.108542
#                    TOB0421           0.000000       897.896510         -0.801505

# Get average difference
### couldn't work out how to pull from summary in python, will do manually - angry emoji


results_brisbane[0].summary_extended()
results_brisbane[1].summary_extended()
results_brisbane[2].summary_extended()


brisbane_logfc_results_dt = pd.DataFrame({'Cell Type': ['FSA0006', 'MBE1006', 'TOB0421'],
        'FSA0006_ref': [0.855849, -0.103652, -0.798455],
        'MBE1006_ref': [0.863100, -0.107455, -0.801995],
        'TOB0421_ref': [0.864988, -0.108542, -0.801505]})

brisbane_logfc_results_dt['mean'] = brisbane_logfc_results_dt.mean(axis=1)

brisbane_logfc_results_dt['is_credible'] = [True, True, True]

brisbane_logfc_results_dt.to_csv(outdir + "brisbane_logfc_results.tsv", sep = "\t")




## Melbourne
data_melbourne = data_all[data_all.obs["Location"].isin(["Melbourne"])]
viz.boxplots(data_melbourne, feature_name="Village")
plt.savefig(outdir + 'melbourne_numbers.png')

model_melbourne = mod.CompositionalAnalysis(data_melbourne, formula="Village", reference_cell_type="automatic")
melbourne_sim_results = model_melbourne.sample_hmc()

melbourne_sim_results.summary()
print(melbourne_sim_results.credible_effects())


# Run scCODA with each cell type as the reference
cell_types_melbourne = data_melbourne.var.index
results_cycle_melbourne = pd.DataFrame(index=cell_types_melbourne, columns=["times_credible"]).fillna(0)

models_melbourne = []
results_melbourne = []

for ct in cell_types_melbourne:
    print(f"Reference: {ct}")
    # Run inference
    model_temp = mod.CompositionalAnalysis(data_melbourne, formula="Village", reference_cell_type=ct)
    models_melbourne.append(model_temp)
    results_temp = model_temp.sample_hmc(num_results=20000)
    results_melbourne.append(results_temp)
    # Select credible effects
    cred_eff = results_temp.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)
    # add up credible effects
    results_cycle_melbourne["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle_melbourne["pct_credible"] = results_cycle_melbourne["times_credible"]/len(cell_types_melbourne)
results_cycle_melbourne["is_credible"] = results_cycle_melbourne["pct_credible"] > 0.5
print(results_cycle_melbourne)


# Get average difference
### couldn't work out how to pull from summary in python, will do manually - angry emoji
results_melbourne[0].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           0.000000      3471.687790          1.147833
#                    MBE1006          -1.930117       648.252082         -1.636737
#                    TOB0421          -1.203827      1065.893461         -0.588923

results_melbourne[1].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           1.932333      3478.980493          1.152990
#                    MBE1006           0.000000       650.122718         -1.634778
#                    TOB0421           0.716773      1056.730123         -0.600693

results_melbourne[2].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           1.139058      3411.776681          1.136354
#                    MBE1006          -0.780606       662.039516         -1.633135
#                    TOB0421           0.000000      1112.017136         -0.506959



# results_melbourne[0].summary_extended().to_csv(outdir + "melbourne_rep1_summary_extended.tsv", sep = "\t")
# results_melbourne[1].summary_extended().to_csv(outdir + "melbourne_rep2_summary_extended.tsv", sep = "\t")
# results_melbourne[2].summary_extended().to_csv(outdir + "melbourne_rep3_summary_extended.tsv", sep = "\t")


results_melbourne[0].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 0
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 211.069 sec. Acceptance rate: 72.1%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              7.019   5.145    8.814  1.027       958.971621
# MBE1006              7.539   5.694    9.373  1.028      1613.016781
# TOB0421              7.509   5.627    9.305  1.027      1565.344931


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006           0.000000   0.000    0.000  0.000                    0.0      1737.951492          0.857828
#                    MBE1006          -0.666177  -0.767   -0.564  0.053                    1.0      1501.597873         -0.103263
#                    TOB0421          -1.150533  -1.264   -1.043  0.059                    1.0       897.783968         -0.802040


results_melbourne[1].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 1
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 216.384 sec. Acceptance rate: 35.2%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              6.137   5.029    6.951  0.693       954.967911
# MBE1006              6.658   5.538    7.457  0.688      1607.889516
# TOB0421              6.637   5.508    7.467  0.703      1574.475907


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006           0.663964   0.565    0.783  0.056                    1.0      1727.811797          0.855422
#                    MBE1006           0.000000   0.000    0.000  0.000                    0.0      1497.642328         -0.102475
#                    TOB0421          -0.475140  -0.597   -0.376  0.055                    1.0       911.879208         -0.787957

results_melbourne[2].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 2
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 215.073 sec. Acceptance rate: 53.5%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              6.318   4.911    7.935  0.975       957.849632
# MBE1006              6.841   5.453    8.440  0.977      1615.970211
# TOB0421              6.808   5.423    8.421  0.979      1563.513490


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006           1.152184   1.044    1.260  0.055                    1.0      1738.891749          0.860297
#                    MBE1006           0.482503   0.388    0.577  0.052                    1.0      1501.653833         -0.105848
#                    TOB0421           0.000000   0.000    0.000  0.000                    0.0       896.787752         -0.801953



melbourne_logfc_results_dt = pd.DataFrame({'Cell Type': ['FSA0006', 'MBE1006', 'TOB0421'],
        'FSA0006_ref': [1.147833, -1.636737, -0.588923],
        'MBE1006_ref': [1.152990, -1.634778, -0.600693],
        'TOB0421_ref': [1.136354, -1.633135, -0.506959]})

melbourne_logfc_results_dt['mean'] = melbourne_logfc_results_dt.mean(axis=1)

melbourne_logfc_results_dt['is_credible'] = [True, True, True]

melbourne_logfc_results_dt.to_csv(outdir + "melbourne_logfc_results.tsv", sep = "\t")






## Sydney - Fresh
data_sydney = data_all[data_all.obs["Location"].isin(["Sydney"]) & data_all.obs["Cryopreservation"].isin(["Fresh"])]
viz.boxplots(data_sydney, feature_name="Village")
plt.savefig(outdir + 'sydney_numbers.png')

model_sydney = mod.CompositionalAnalysis(data_sydney, formula="Village", reference_cell_type="TOB0421")

sydney_sim_results = model_sydney.sample_hmc()
sydney_sim_results.summary()
print(sydney_sim_results.credible_effects())


# Run scCODA with each cell type as the reference
cell_types_sydney = data_sydney.var.index
results_cycle_sydney = pd.DataFrame(index=cell_types_sydney, columns=["times_credible"]).fillna(0)

models_sydney = []
results_sydney = []

for ct in cell_types_sydney:
    print(f"Reference: {ct}")
    # Run inference
    model_temp = mod.CompositionalAnalysis(data_sydney, formula="Village", reference_cell_type=ct)
    models_sydney.append(model_temp)
    results_temp = model_temp.sample_hmc(num_results=20000)
    results_sydney.append(results_temp)
    # Select credible effects
    cred_eff = results_temp.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)
    # add up credible effects
    results_cycle_sydney["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle_sydney["pct_credible"] = results_cycle_sydney["times_credible"]/len(cell_types_sydney)
results_cycle_sydney["is_credible"] = results_cycle_sydney["pct_credible"] > 0.5
print(results_cycle_sydney)

# Get average difference
### couldn't work out how to pull from summary in python, will do manually - angry emoji
results_sydney[0].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006                0.0      1130.524812               0.0
#                    MBE1006                0.0      1255.685903               0.0
#                    TOB0421                0.0       819.289285               0.0

results_sydney[1].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006                0.0      1132.029342               0.0
#                    MBE1006                0.0      1244.846090               0.0
#                    TOB0421                0.0       828.624568               0.0

results_sydney[2].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006                0.0      1125.134806               0.0
#                    MBE1006                0.0      1255.963310               0.0
#                    TOB0421                0.0       824.401884               0.0



results_sydney[0].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 0
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 201.450 sec. Acceptance rate: 58.3%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              5.127   3.915    6.303  0.625      1131.273416
# MBE1006              5.229   3.959    6.361  0.629      1252.753483
# TOB0421              4.807   3.663    6.004  0.621       821.473102


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006                0.0   0.000    0.000  0.000               0.000000      1131.273416               0.0
#                    MBE1006                0.0  -0.256    0.028  0.072               0.405067      1252.753483               0.0
#                    TOB0421                0.0  -0.170    0.150  0.043               0.274533       821.473102               0.0

results_sydney[1].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 1
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 201.570 sec. Acceptance rate: 55.7%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              5.004   3.969    6.005  0.546      1131.297199
# MBE1006              5.100   4.049    6.086  0.548      1245.285646
# TOB0421              4.693   3.652    5.684  0.547       828.917155


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006                0.0  -0.066    0.233  0.061               0.338667      1131.297199               0.0
#                    MBE1006                0.0   0.000    0.000  0.000               0.000000      1245.285646               0.0
#                    TOB0421                0.0  -0.131    0.190  0.046               0.284733       828.917155               0.0

results_sydney[2].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 2
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 197.450 sec. Acceptance rate: 52.1%

# Intercepts:
#            Final Parameter  HDI 3%  HDI 97%     SD  Expected Sample
# Cell Type                                                          
# FSA0006              5.077   4.100    6.166  0.587      1126.016304
# MBE1006              5.185   4.187    6.299  0.596      1254.435925
# TOB0421              4.766   3.769    5.855  0.591       825.047771


# Effects:
#                               Final Parameter  HDI 3%  HDI 97%     SD  Inclusion probability  Expected Sample  log2-fold change
# Covariate          Cell Type                                                                                                   
# Village[T.Village] FSA0006                0.0  -0.069    0.221  0.058               0.339867      1126.016304               0.0
#                    MBE1006                0.0  -0.242    0.059  0.068               0.392600      1254.435925               0.0
#                    TOB0421                0.0   0.000    0.000  0.000               0.000000       825.047771               0.0



sydney_logfc_results_dt = pd.DataFrame({'Cell Type': ['FSA0006', 'MBE1006', 'TOB0421'],
        'FSA0006_ref': [0, 0, 0],
        'MBE1006_ref': [0, 0, 0],
        'TOB0421_ref': [0, 0, 0]})

sydney_logfc_results_dt['mean'] = sydney_logfc_results_dt.mean(axis=1)

sydney_logfc_results_dt['is_credible'] = [False, False, False]

sydney_logfc_results_dt.to_csv(outdir + "sydney_logfc_results.tsv", sep = "\t")





## Sydney - Cryopreserved
data_sydney_cryo = data_all[data_all.obs["Location"].isin(["Sydney"]) & data_all.obs["Cryopreservation"].isin(["Cryopreserved"])]
viz.boxplots(data_sydney_cryo, feature_name="Village")
plt.savefig(outdir + 'sydney_numbers.png')

model_sydney_cryo = mod.CompositionalAnalysis(data_sydney_cryo, formula="Village", reference_cell_type="automatic")

sydney_sim_results_cryo = model_sydney_cryo.sample_hmc()
sydney_sim_results_cryo.summary()
print(sydney_sim_results_cryo.credible_effects())


# Run scCODA with each cell type as the reference
cell_types_sydney_cryo = data_sydney_cryo.var.index
results_cycle_sydney_cryo = pd.DataFrame(index=cell_types_sydney_cryo, columns=["times_credible"]).fillna(0)

models_sydney_cryo = []
results_sydney_cryo = []

for ct in cell_types_sydney_cryo:
    print(f"Reference: {ct}")
    # Run inference
    model_temp = mod.CompositionalAnalysis(data_sydney_cryo, formula="Village", reference_cell_type=ct)
    models_sydney_cryo.append(model_temp)
    results_temp = model_temp.sample_hmc(num_results=20000)
    results_sydney_cryo.append(results_temp)
    # Select credible effects
    cred_eff = results_temp.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)
    # add up credible effects
    results_cycle_sydney_cryo["times_credible"] += cred_eff.astype("int")

# Calculate percentages
results_cycle_sydney_cryo["pct_credible"] = results_cycle_sydney_cryo["times_credible"]/len(results_cycle_sydney_cryo)
results_cycle_sydney_cryo["is_credible"] = results_cycle_sydney_cryo["pct_credible"] > 0.5
print(results_cycle_sydney_cryo)

# Get average difference
### couldn't work out how to pull from summary in python, will do manually - angry emoji
results_sydney_cryo[0].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006            0.00000      1620.195457          1.368064
#                    MBE1006           -2.05755       294.934049         -1.600354
#                    TOB0421           -1.82538       280.037161         -1.265402

results_sydney_cryo[1].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           1.967652      1628.040102          1.372877
#                    MBE1006           0.000000       308.733991         -1.465845
#                    TOB0421           0.000000       258.392574         -1.465845

results_sydney_cryo[2].summary()
# Effects:
#                               Final Parameter  Expected Sample  log2-fold change
# Covariate          Cell Type                                                    
# Village[T.Village] FSA0006           1.925718      1608.624854          1.359406
#                    MBE1006           0.000000       320.030556         -1.418819
#                    TOB0421           0.000000       266.511256         -1.418819



results_sydney_cryo[0].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 0
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 211.440 sec. Acceptance rate: 74.1%

# Intercepts:
#            Final Parameter  ...  Expected Sample
# Cell Type                   ...                 
# FSA0006              4.711  ...       631.273016
# MBE1006              5.057  ...       892.242932
# TOB0421              4.773  ...       671.650718

# [3 rows x 5 columns]


# Effects:
#                               Final Parameter  ...  log2-fold change
# Covariate          Cell Type                   ...                  
# Village[T.Village] FSA0006           0.000000  ...          1.359483
#                    MBE1006          -2.051918  ...         -1.600809
#                    TOB0421          -1.813005  ...         -1.256130

# [3 rows x 7 columns]



results_sydney_cryo[1].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 cell types
# Reference index: 1
# Formula: Village
# Spike-and-slab threshold: 1.000

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 214.655 sec. Acceptance rate: 44.2%

# Intercepts:
#            Final Parameter  ...  Expected Sample
# Cell Type                   ...                 
# FSA0006              3.108  ...       629.681591
# MBE1006              3.412  ...       853.387975
# TOB0421              3.231  ...       712.097101

# [3 rows x 5 columns]


# Effects:
#                               Final Parameter  ...  log2-fold change
# Covariate          Cell Type                   ...                  
# Village[T.Village] FSA0006           1.966162  ...          1.370768
#                    MBE1006           0.000000  ...         -1.465804
#                    TOB0421           0.000000  ...         -1.465804

# [3 rows x 7 columns]

results_sydney_cryo[2].summary_extended()
# Compositional Analysis summary (extended):

# Data: 6 samples, 3 ce
# ll types
# Reference index: 2
# Formula: Village
# Spike-and-slab threshold: 0.997

# MCMC Sampling: Sampled 20000 chain states (5000 burnin samples) in 219.024 sec. Acceptance rate: 45.9%

# Intercepts:
#            Final Parameter  ...  Expected Sample
# Cell Type                   ...
# FSA0006              3.152  ...       629.950098
# MBE1006              3.464  ...       860.609283
# TOB0421              3.264  ...       704.607286

# [3 rows x 5 columns]


# Effects:
#                               Final Parameter  ...  log2-fold change
# Covariate          Cell Type                   ...
# Village[T.Village] FSA0006           1.909161  ...          1.348700
#                    MBE1006           0.000000  ...         -1.405636
#                    TOB0421           0.000000  ...         -1.405636

# [3 rows x 7 columns]



sydney_logfc_results_dt_cryo = pd.DataFrame({'Cell Type': ['FSA0006', 'MBE1006', 'TOB0421'],
        'FSA0006_ref': [1.368064, -1.600354, -1.265402],
        'MBE1006_ref': [1.372877, -1.465845, -1.465845],
        'TOB0421_ref': [1.359406, -1.418819, -1.418819]})

sydney_logfc_results_dt_cryo['mean'] = sydney_logfc_results_dt_cryo.mean(axis=1)

sydney_logfc_results_dt_cryo['is_credible'] = [True, False, False]

sydney_logfc_results_dt_cryo.to_csv(outdir + "sydney_cryo_logfc_results.tsv", sep = "\t")

