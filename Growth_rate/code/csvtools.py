"""
csv manipulation tools
"""

import csv
from copy import deepcopy
from os import path

import click
import toml


@click.group()
def main():
    """
    tools for manipulating csv files during data preprocessing
    """
    pass


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def zero_time_split(infile, outfile):
    """
    Set minimum time to 0 and adjust other times to match
    """
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            times = []
            for line in rdr:
                times.append(float(line['time']))
            min_time = min(times)
            # go back to start for copying with modification
            in_csv.seek(0)
            rdr.__next__() # skip header
            wtr.writeheader()
            for line in rdr:
                line['time'] = float(line['time']) - min_time
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def zero_time_longform(infile, outfile):
    """
    set minimum time to 0 (and adjust others accordingly)
    treats each set with identical name, well and generation as one dataset
    """
    min_times = {}
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            # times = []
            for line in rdr:
                # times.append(float(line['time']))
                # id_string = '.'.join([line[x] for x in ['name', 'generation', 'well']])
                id_string = line['name']
                if id_string not in min_times:
                    min_times[id_string] = []
                min_times[id_string].append(float(line['time']))
            # min_time = min(times)
            min_times = {k: min(v) for k, v in min_times.items()}
            # go back to start for copying with modification
            in_csv.seek(0)
            rdr.__next__() # skip header
            wtr.writeheader()
            for line in rdr:
                # id_string = '.'.join([line[x] for x in ['name', 'generation', 'well']])
                id_string = line['name']
                line['time'] = float(line['time']) - min_times[id_string]
                wtr.writerow(line)



@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfolder', type=click.Path())
def split_longform(infile, outfolder):
    """
    decompose a single longform csv into one small aptly named file for each cell colony
    """
    data = {}
    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for line in rdr:
            idx = line['name']
            if idx not in data:
                data[idx] = []
            data[idx].append(line)

    for dataset in data.values():
        filename = outfolder + '.'.join([dataset[0]['name'],
                                         'csv'])
        with open(filename, 'w') as out_csv:
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            wtr.writeheader()
            for line in dataset:
                del line['name']
                del line['well']
                del line['generation']
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
@click.option('-n', 'name', type=str)
@click.option('-p', 'position', type=int)
@click.option('-v', 'init_val', type=str)
def add_column(infile, outfile, name, position, init_val):
    """
    add a column named 'name' to a csv
    position if given will be the position of the new column
    """
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            fieldnames = rdr.fieldnames[:]
            if name in fieldnames:
                print("column with given name already in csv")
                return
            if position is not None:
                fieldnames.insert(position, name)
            else:
                fieldnames.append(name)
            wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
            wtr.writeheader()
            for line in rdr:
                line[name] = init_val
                wtr.writerow(line)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
@click.option('-n', 'name', type=str)
def delete_column(infile, outfile, name):
    """
    delete a column named 'name' from a csv
    """
    with open(infile, 'r') as in_csv:
        with open(outfile, 'w') as out_csv:
            rdr = csv.DictReader(in_csv)
            fieldnames = rdr.fieldnames[:]
            if name in fieldnames:
                exists = True
                fieldnames.remove(name)
            else:
                exists = False
                print("column with given name not in csv:", name)
            wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
            wtr.writeheader()
            for line in rdr:
                print(line)
                if exists:
                    del line[name]
                wtr.writerow(line)



@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-p', 'paramfile', type=click.Path())
@click.option('-o', 'outfile', type=click.Path())
def define_groups(infile, paramfile, outfile):
    """
    define the coupling groups
    """

    data = {}
    ignored = set()

    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for l in rdr:
            # id_string = '.'.join([l[x] for x in ['name', 'generation', 'well']])
            id_string = l['name']
            if id_string not in data:
                data[id_string] = {str(x): [] for x in rdr.fieldnames}
            for k, v in l.items():
                data[id_string][k].append(v)

    params = toml.load(paramfile)
    birthrate_coupling = params['abc_params']['birthrate_coupling_sets']


    running_birthrate_group = 0
    if not birthrate_coupling or birthrate_coupling == 'none':
        # subdivide completely by names
        for id_string, obs in data.items():
            obs['birthrate_group'] = running_birthrate_group
            running_birthrate_group += 1
            print(obs['birthrate_group'])
    elif birthrate_coupling == 'all':
        for id_string, obs in data.items():
            obs['birthrate_group'] = 0
    else:
        for id_string, obs in data.items():
            try:
                obs['birthrate_group'] = [id_string in x for x in birthrate_coupling].index(True)
                print('Creating group for name', id_string)
            except:
                print('Ignoring data rows', id_string, 'reason: not in coupling sets')
                ignored.add(id_string)

    with open(outfile, 'w') as out_csv:
        fieldnames = list(data[list(data.keys())[0]].keys())
        wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
        wtr.writeheader()
        for __, obs in data.items():
            # print(obs)
            if obs['name'][0] in ignored:
                continue
            for i, __ in enumerate(obs['time']):
                row = {k: obs[k][i] for k in fieldnames if k[-5:] != 'group'}
                row['birthrate_group'] = obs['birthrate_group']
                wtr.writerow(row)


@main.command()
@click.option('-i', 'infile', type=click.Path())
@click.option('-p', 'paramfile', type=click.Path())
@click.option('-z', 'outdir_base', type=click.Path())
def split_by_group(infile, paramfile, outdir_base):
    """
    decompose a single longform csv based on the groups
    also decomposes the paramfile to match
    """
    data = {}
    outfilebase = outdir_base + 'intermediate/' + path.split(paramfile)[1].split('.')[0]
    print(outfilebase)
    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for line in rdr:
            idx = line['birthrate_group']
            if idx not in data:
                data[idx] = []
            data[idx].append(line)

    params = toml.load(paramfile)


    sets = deepcopy(params['abc_params']['birthrate_coupling_sets'])

    if sets == 'none':
        sets = [dataset[0]['name'] for key, dataset in data.items()]



    groups = {}


    for key, dataset in data.items():
        print(key, dataset)
        filename = outfilebase + '.g' + key + '.data.csv'
        pfn = outfilebase + '.g' + key + '.toml'
        groups[key] = {'obs': filename, 'par': pfn}
        params['abc_params']['birthrate_coupling'] = 'all'
        params['abc_params']['birthrate_coupling_sets'] = []
        with open(filename, 'w') as out_csv:
            wtr = csv.DictWriter(out_csv, fieldnames=rdr.fieldnames)
            wtr.writeheader()
            for line in dataset:
                wtr.writerow(line)
        with open(pfn, 'w') as out_toml:
            new_params = deepcopy(params)
            if 'plot_params' not in params:
                params['plot_params'] = {}
                new_params['plot_params'] = {}
            if 'coupling_names' in params['plot_params']:
                new_params['plot_params']['coupling_names'] = params['plot_params']['coupling_names'][int(key)]
            else:
                new_params['plot_params']['coupling_names'] = ' '.join(sets[int(key)])
            toml.dump(new_params, out_toml)


@main.command()
@click.option('-i', 'infiles', type=click.Path(), multiple=True)
@click.option('-o', 'outfile', type=click.Path())
def merge(infiles, outfile):
    with open(infiles[0], 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        fieldnames = rdr.fieldnames[:]

    with open(outfile, 'w') as out_csv:
        wtr = csv.DictWriter(out_csv, fieldnames=fieldnames)
        wtr.writeheader()
        for infile in infiles:
            with open(infile, 'r') as in_csv:
                rdr = csv.DictReader(in_csv)
                # print(rdr.fieldnames, fieldnames)
                assert rdr.fieldnames == fieldnames
                for row in rdr:
                    wtr.writerow(row)


if __name__ == '__main__':
    main()
