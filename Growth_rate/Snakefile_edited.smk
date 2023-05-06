import glob
from os import path


indir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/ratrack/data/"
outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/ratrack/"

# filenames = ['166_ratrack']

files = glob.glob(indir + "*.csv")
filenames = [os.path.basename(x).replace(".csv","") for x in files]
print(filenames)
# filenames = path.split(glob.glob(indir + "*.csv")).split('.')[0] 
# print(filenames)


rule all:
    input:
        expand(outdir + "results/{filename}.fit.csv", filename = filenames),
        expand(outdir + "results/{filename}.pdf", filename = filenames)
        # expand(outdir + "intermediate/{filename}.groups.csv", filename = filenames)


rule define_groups:
    input:
        obs = indir + "{filename}.csv",
        par = indir + "{filename}.toml"
    output:
        temp(outdir + "intermediate/{filename}.groups.csv")
        # outdir + "intermediate/{filename}.groups.csv"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/csvtools.py define-groups \
            -i {input.obs} \
            -p {input.par} \
            -o {output}
        """


rule zero_time:
    input:
        outdir + "intermediate/{filename}.groups.csv"
    output:
        temp(outdir + "intermediate/{filename}.groups.zero.csv")
        # "intermediate/{filename}.groups.zero.csv"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/csvtools.py zero-time-longform \
            -i {input} \
            -o {output}
        """


rule discard_deaths:
    input:
        outdir + "intermediate/{filename}.groups.zero.csv"
    output:
        temp(outdir + "intermediate/{filename}.groups.zero.deathless.csv")
        # "intermediate/{filename}.groups.zero.deathless.csv"
        # "intermediate/{filename}.groups.zero.csv"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/csvtools.py delete-column \
            -i {input} \
            -o {output} \
            -n dead
        """


rule split_by_group:
    input:
        obs = outdir + "intermediate/{filename}.groups.zero.deathless.csv",
        par = indir + "{filename}.toml"
    output:
        dynamic(outdir + "intermediate/{filename}.g{k}.data.csv"),
        dynamic(outdir + "intermediate/{filename}.g{k}.toml"),
        # groupinfo = "intermediate/{filename}.groupinfo.toml",
    params:
        outdir = outdir
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/csvtools.py split-by-group \
            -i {input.obs} \
            -p {input.par} \
            -z {params.outdir}
        """


rule reconstruct:
    output:
        outdir + "intermediate/{filename}.g{k}.db"
    input:
        obs = outdir + "intermediate/{filename}.g{k}.data.csv",
        par = outdir + "intermediate/{filename}.g{k}.toml"
    run:
        shell(" \
             /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/abc.py reconstruct \
                 -p {input.par} \
                 -o {input.obs} \
                 -d {output} \
        ")


rule abc_plots:
    input:
        db = outdir + "intermediate/{filename}.g{k}.db",
        obs = outdir + "intermediate/{filename}.g{k}.data.csv",
        par = outdir + "intermediate/{filename}.g{k}.toml"
    output:
        outdir + "intermediate/{filename}.g{k}.abc.pdf"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/plots.py abc-info \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            --save {output}
        """


rule fit_plots:
    input:
        db = outdir + "intermediate/{filename}.g{k}.db",
        obs = outdir + "intermediate/{filename}.g{k}.data.csv",
        par = outdir + "intermediate/{filename}.g{k}.toml"
    output:
        outdir + "intermediate/{filename}.g{k}.fit.pdf"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/plots.py result-single \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            --save {output}
        """


rule fit_tables:
    input:
        db = outdir + "intermediate/{filename}.g{k}.db",
        obs = outdir + "intermediate/{filename}.g{k}.data.csv",
        par = outdir + "intermediate/{filename}.g{k}.toml"
    output:
        outdir + "intermediate/{filename}.g{k}.fit.csv"
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/plots.py tabulate-single \
            -p {input.par} \
            -o {input.obs} \
            -d {input.db} \
            -c {output}
        """


def report_inputs(wildcards):
    import toml
    params = toml.load(indir + wildcards.filename + '.toml')
    if params['abc_params']['birthrate_coupling_sets'] == 'all':
        groups = [0]
    elif params['abc_params']['birthrate_coupling_sets'] == 'none':
        with open(outdir +  wildcards.filename + '.csv') as in_csv:
            names = []
            for row in in_csv:
                names.append(row.split(',')[0])
            groups = list(range(len(set(names)) - 1))
    else:
        groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
    print(params)
    print(groups)
    sources = {}
    for k in groups:
        sources['g' + str(k) + '.abc'] = outdir + 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.abc.pdf'
        sources['g' + str(k) + '.fit'] = outdir + 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.fit.pdf'
        sources['g' + str(k) + '.db'] = outdir + 'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.db'
    print("collected sources for report:", sources)
    return sources



def table_inputs(wildcards):
    import toml
    params = toml.load(indir + wildcards.filename + '.toml')
    # groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
    if params['abc_params']['birthrate_coupling_sets'] == 'all':
        groups = [0]
    elif params['abc_params']['birthrate_coupling_sets'] == 'none':
        with open(indir + wildcards.filename + '.csv') as in_csv:
            names = []
            for row in in_csv:
                names.append(row.split(',')[0])
            # print(set(names), len(set(names)))
            groups = list(range(len(set(names)) - 1))
    else:
        groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
    sources = {}
    for k in groups:
        # sources['g' + str(k) + '.db'] = 'intermediate/' \
        #     + wildcards.filename + '.g' + str(k) + '.db'
        sources['g' + str(k) + '.csv'] = outdir +'intermediate/' \
            + wildcards.filename + '.g' + str(k) + '.fit.csv'
    print("collected sources for table:", sources)
    return sources



def num_to_aa(n):
    """
    convert an integer to a base 25 letter number
    NOTE: this function isn't actually correct
    but it's good enough for this use-case
    """
    if n == 0:
        return 'A'
    aa = ''
    while n > 0:
        aa += chr(65 + n%25)
        n //= 25
    return aa


rule produce_table:
    input:
        unpack(table_inputs)
    output:
        outdir + "results/{filename}.fit.csv"
    run:
        print(input)
        tables = ' '.join(['-i ' + x for x in input])
        shell('/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/ratrack/bin/python /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/code/csvtools.py merge -o ' + str(output) + ' ' + tables)


rule produce_report:
    input:
        unpack(report_inputs)
    output:
        outdir + "results/{filename}.pdf"
    run:
        # identify most likely model to put into report
        import toml
        import re
        import numpy as np
        from pyabc import History
        params = toml.load(indir + wildcards.filename + '.toml')
        if params['abc_params']['birthrate_coupling_sets'] == 'all':
            groups = [0]
        else:
            groups = list(range(len(params['abc_params']['birthrate_coupling_sets'])))
        sources = {x: {} for x in groups}
        print(groups)
        for group in groups:
            print(group)
            sources[group]['names'] = params['abc_params']['birthrate_coupling_sets'][group]
            db_path = 'sqlite:///' + input['g' + str(group) + '.db']
            abc_history = History(db_path)
            axs = abc_history.get_model_probabilities()
            final = np.array(axs[-1:])
            final = final[final > 0]
            sources[group]['model_fraction'] = np.max(final)
            sources[group]['best_model'] = np.where(final == sources[group]['model_fraction'])[0][0]
            sources[group]['num_models'] = len(final)

        # collect sources into output file
        files = ' '.join([num_to_aa(x) + '=' + input['g' + str(x) + '.fit']
                          for x in groups])
        pages = ' '.join([' '.join([num_to_aa(y) \
                                    + str(x*sources[y]['num_models'] + sources[y]['best_model'] + 1)
                          for x in range(2)]) for y in groups])
        command = 'pdftk ' + files + ' cat ' + pages + ' output ' + str(output)
        shell(command)


