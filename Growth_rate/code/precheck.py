"""
rudimentary test of input files to verify that they conform to expected content/format
"""


import toml
import csv
import click


@click.group()
def main():
    """
    input file testing
    """
    pass


@main.command()
@click.option('-p', 'paramfile', type=click.Path())
@click.option('-o', 'obsfile', type=click.Path())
def test(paramfile, obsfile, minimal=True):
    """
    input file testing
    """

    print('Testing parameter file:', paramfile)
    params = toml.load(paramfile)
    if not minimal:
        assert isinstance(params['simulation_params']['starting_cell_count'], (int, float)) or \
           params['simulation_params']['starting_cell_count'] == 'calculate'
        assert isinstance(params['simulation_params']['end_time'], (int, float)) or \
           params['simulation_params']['end_time'] == 'max_observed'

    if 'deathrate_interaction' in params['simulation_params']:
        assert isinstance(params['simulation_params']['deathrate_interaction'], (int, float))
        print('carrying capacity is:', 1/params['simulation_params']['deathrate_interaction'])
    else:
        assert params['simulation_params']['carrying_capacity'] > 0
        print('carrying capacity is:', params['simulation_params']['carrying_capacity'])

    if not minimal:
        assert isinstance(params['abc_params']['starting_population_size'], int)
        assert params['abc_params']['starting_population_size'] > 0
        assert isinstance(params['abc_params']['min_epsilon'], (int, float))
        assert params['abc_params']['min_epsilon'] > 0
        assert isinstance(params['abc_params']['max_populations'], (int))
        assert params['abc_params']['max_populations'] > 0
        assert isinstance(params['abc_params']['min_acceptance'], (int, float))
        assert params['abc_params']['min_acceptance'] >= 0

        if params['abc_params']['max_populations'] > 20:
            print('high ABC generations:',
                params['abc_params']['max_populations'],
                'simulations may take a long time')

        if params['abc_params']['starting_population_size'] > 1000:
            print('high ABC particle count:',
                params['abc_params']['starting_population_size'],
                'simulations may take a long time')

    assert params['abc_params']['simulator'] in ['rar-engine', 'bernoulli']

    assert len(params['abc_params']['rate_limits']) == 2
    for rate_limit in params['abc_params']['rate_limits']:
        if params['abc_params']['simulator'] == 'rar-engine':
            assert rate_limit > 0

    assert len(params['abc_params']['resolution_limits']) == 2
    for resolution_limit in params['abc_params']['resolution_limits']:
        assert resolution_limit > 0

    assert isinstance(params['abc_params']['parallel_simulations'], int)
    assert params['abc_params']['parallel_simulations'] > 0

    if not minimal:
        if params['abc_params']['birthrate_coupling_sets'] not in ['all', 'none']:
            assert len(params['plot_params']['coupling_names']) == len(params['abc_params']['birthrate_coupling_sets'])

    for filt in params['filters']:
        assert filt['name'] in ['copy', 'perfect', 'poisson', 'gauss-multiplicative', 'gauss-additive']
        # TODO test for presence of correct filter parameters

    with open(obsfile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for column in ['name', 'time', 'count']:
            assert column in rdr.fieldnames
        if 'sample' in rdr.fieldnames:
            print('column \'sample\' is present, but program looks for numbered samplings starting from 1')
            print('for instance: sample1, sample2, etc...')
        if 'dilution' in rdr.fieldnames:
            print('column \'dilution\' is present, but program looks for numbered dilutions starting from 1')
            print('for instance: dilution1, dilution2, etc...')


if __name__ == '__main__':
    main()
