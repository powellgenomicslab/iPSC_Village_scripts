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

sydney_logfc_results_dt_cryo = pd.DataFrame({'Cell Type': ['FSA0006', 'MBE1006', 'TOB0421'],
        'FSA0006_ref': [1.368064, -1.600354, -1.265402],
        'MBE1006_ref': [1.372877, -1.465845, -1.465845],
        'TOB0421_ref': [1.359406, -1.418819, -1.418819]})

sydney_logfc_results_dt_cryo['mean'] = sydney_logfc_results_dt_cryo.mean(axis=1)

sydney_logfc_results_dt_cryo['is_credible'] = [True, False, False]

sydney_logfc_results_dt_cryo.to_csv(outdir + "sydney_cryo_logfc_results.tsv", sep = "\t")

