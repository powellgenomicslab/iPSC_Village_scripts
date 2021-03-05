# Prepare inputs
import os
import pandas as pd
from glob import glob
import re
import subprocess
import gzip

def matchFolders(x, dir_list = None):
    for folder in dir_list:
        if re.search(r'^' + x + "$", folder):
            return(folder)
        elif re.search(r'^' + x + '\D', folder):
            return(folder)
        elif re.search(x + "$", folder):
            return(folder)
        elif re.search(x + "\D", folder):
            return(folder)

def get_barcodes_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if f.endswith("barcodes.tsv.gz")]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))

def get_bam_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if f.endswith(".bam")]:
            return(os.path.join(dirpath, filename))

def get_matrix_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("matrix.mtx", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath, filename))

def get_matrix_dirs(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if re.search("matrix.mtx", f)]:
            if re.search(r'filtered', os.path.join(dirpath, filename)):
                return(os.path.join(dirpath))

def get_expected_Ndoublet(pool):
    if pool.endswith(".gz"):
        lines_in_file = gzip.open(pool, 'r').readlines()
        number_of_doublets = len(lines_in_file)*0.008
    else:
        lines_in_file = open(pool, 'r').readlines()
        number_of_doublets = ((len(lines_in_file)**2)*0.008)/1000
    return(number_of_doublets)


def get_scrnaseq_dirs(scrnaseq_dir, sample_file):
    ### Check that all the directories exist and can be accessed
    if not os.path.exists(scrnaseq_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(scrnaseq_dir))

    ### Read in samplesheet from the configfile specified from the command line
    samples = pd.read_csv(sample_file, sep = "\t")

    ### Expect first colunn to be pools
    pools = samples.iloc[:, 0]

    ### Match pools to scrna seq directories to make a list of each scRNA-seq dir
    scrna_seq_dirlist = os.listdir(scrnaseq_dir)
    try:
        scrnaseq_filelist = [os.path.join(scrnaseq_dir, matchFolders(pool, dir_list = scrna_seq_dirlist)) for pool in pools]
    except TypeError:
        print(error)
        print("Could not find a scRNA-seq directory for all of the pools in your pool list. Please check that they are spelled correctly and you do not have any additional pool names that are not in {}  ".format(individual_list_dir))
    scrnaseq_filedict = dict(zip(pools, scrnaseq_filelist))
    scrnaseq_libs = pd.Series(scrnaseq_filedict, name="scRNAseq_Directories")

    ### Get the directories of all the barcodes including the files themselves
    try:
        barcode_filelist = [get_barcodes_files(pool) for pool in scrnaseq_filelist]
    except:
        print(error)
        print("Could not find a barcode file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'barcodes.tsv' within the name.")
    barcode_filedict = dict(zip(pools, barcode_filelist))
    barcode_libs = pd.Series(barcode_filedict, name="Barcode_Files")

    ### Get the directories of all the bams including the files themselves
    try:
        bam_filelist = [get_bam_files(pool) for pool in scrnaseq_filelist]
    except:
        print(error)
        print("Could not find a bam file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain '.bam' within the name.")
    bam_filedict = dict(zip(pools, bam_filelist))
    bamlibs = pd.Series(bam_filedict, name="Bam_Files")

    ### Get the directories of all the matrix files including the files themselves
    try:
        matrix_filelist = [get_matrix_files(pool) for pool in scrnaseq_filelist]
    except:
        print(error)
        print("Could not find a matrix file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'matrix.mtx' within the name.")
    matrix_filedict = dict(zip(pools, matrix_filelist))
    matrix_libs = pd.Series(matrix_filedict, name="Matrix_Files")

    ### Get the directories of all the matrix EXCLUDING the matrix filenames themselves
    try:
        matrix_dirlist = [get_matrix_dirs(pool) for pool in scrnaseq_filelist]
    except:
        print(error)
        print("Could not find a matrix file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain 'matrix.mtx' within the name.")
    matrix_dirdict = dict(zip(pools, matrix_dirlist))
    matrix_dir_libs = pd.Series(matrix_dirdict, name="Matrix_Directories")

    number_doublets = [get_expected_Ndoublet(pool) for pool in barcode_libs]
    number_doubletsdict = dict(zip(pools, number_doublets))
    number_doublets_lib = pd.Series(number_doubletsdict, name="Number_Expected_Doublets")

    dataframe=pd.concat([scrnaseq_libs, barcode_libs, bamlibs, matrix_libs, matrix_dir_libs, number_doublets_lib.astype(int)], axis=1)
    return(dataframe)

