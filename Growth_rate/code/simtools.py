"""
Shared tools for abc simulation and analysis
"""

import copy
import csv
# import statistics
import subprocess
import sys
from io import StringIO

import numpy as np
import toml


FORWARD_SAMPLING = 'RV'
BACKWARD_SAMPLING = 'MLE'


# change to 0 for less output, and 2 or 3 for more output
VERBOSITY = 1


# simulate a lb-process using the given parameters with external software
# n - starting number of cells
# t - series of time points when population will be measured (have to include 0)
# b - birth rate, can be a number or a list. In case of list, it is evenly spread over the timeline
#     with interpolation
# q - interaction death rate. Complementary part of quadratic term that just works by increasing
#     death rate.
# Returns a simulated timeline, calculated using c++ software.
# (Stochastic simulation using next reaction method)
# returns 3 vectors, time, size, rate
# time - time point for this datapoint
# size - size at that timepoint
# rate - growth rate at that timepoint (nice for visualizing interpolation)
def simulate_timeline(starting_population,
                      times,
                      birthrates,
                      deathrate_interaction,
                      simulator,
                      verbosity=VERBOSITY):
    """
    Simulate a lb-process using external software
    """
    # sanity checking
    assert starting_population > 0
    if simulator != 'bernoulli':
        # bernoulli simulator handles negative rates without issues
        # NOTE remember to update this if if neccessary
        for birthrate in birthrates:
            assert birthrate >= 0
    assert deathrate_interaction >= 0
    times = sorted(times)
    for t in times:
        assert t >= 0
    # run external

    cmd = 'code/bin/' + simulator + \
          ' -n ' + str(starting_population) + \
          ' -t \'' + str(['{:f}'.format(x) for x in times]) + '\'' \
          ' -b \'' + str(['{:f}'.format(x) for x in birthrates]) + '\'' \
          ' -q ' + str(deathrate_interaction)
          # ' -t \'' + str(times) + '\'' \
          # ' -b \'' + str(birthrates) + '\'' \
    if verbosity > 0:
        print(cmd)
    output = subprocess.getoutput(cmd)
    if verbosity > 1:
        print(output)
    # output from rar-engine is actually a .tsv file (printed in stdout)
    # use some trickery to parse it with csv
    buff = StringIO(output)
    time = []
    size = []
    rate = []

    rdr = csv.DictReader(buff, dialect='excel-tab')
    try:
        for line in rdr:
            time.append(float(line['time']))
            size.append(int(float(line['size'])))
            # Casting like this can maybe lose precision.
            # But python3 doesn't want to construct integers from scientific notation,
            # whereas the float-constructor handles anything TODO fix?
            rate.append(float(line['rate']))
    except ValueError:
        print('Timeline simulation does not conform to standard', file=sys.stderr)
        print(output, file=sys.stderr)
        exit(1)

    if verbosity > 2:
        print([x for x in zip(time, size, rate)])

    return np.array(time), np.array(size), np.array(rate)


def apply_noise(size, filters):
    """
    Apply list of noise filters in order.
    Implemented types:
      copy: filter does nothing
      perfect: simulates perfect sampling
      poisson: simulates random sampling (in any number of steps)
      gauss-multiplicative: gaussian noise with constant COV
      gauss-additive: gaussian noise with constant stdev
    """
    for filt in filters:
        if filt['name'] == ['copy']:
            continue
        elif filt['name'] == 'perfect':
            size *= filt['sample']
            # size = [x * y for x, y in zip(size, filt['sample'])]
        elif filt['name'] == 'poisson':
            size = np.random.poisson(size * filt['sample'])
            # size = [np.random.poisson(x * y) for x, y in zip(size, filt['sample'])]
        elif filt['name'] == 'gauss-multiplicative':
            size = np.round(size*np.random.normal(filt['mean'], filt['sigma'], size.size), 0).astype(int)
            # size = [np.random.normal(filt['mean'], filt['sigma']) * x for x in size]
        elif filt['name'] == 'gauss-additive':
            size += np.random.normal(filt['mean'], filt['sigma'], size.size)
            # size = [np.random.normal(filt['mean'], filt['sigma']) + x for x in size]

    return np.round(size)


PARAMS = {}

OBSERVED = {}


def get_samplings_dilutions(observed):
    """
    Get the list of all samplings and dilutions done to particular observation
    """
    samplings = [[] for __ in observed['time']]
    dilutions = [[] for __ in observed['time']]
    if VERBOSITY > 2:
        print('obs', observed)
    for j, __ in enumerate(observed['time']):
        i = 1
        while True:
            if 'sample' + str(i) in observed.keys():
                samplings[j].append(observed['sample' + str(i)][j])
            else:
                break
            i += 1
        i = 1
        while True:
            if 'dilute' + str(i) in observed.keys():
                dilutions[j].append(observed['dilute' + str(i)][j])
            else:
                break
            i += 1

    if VERBOSITY > 2:
        print(np.array(samplings), np.array(dilutions))
    return np.array(samplings), np.array(dilutions)
    # return np.array(zip(*samplings)), np.array(zip(*dilutions))


def apply_sampling(size, samplings, dilutions):
    """
    apply sampling methods to simulated data to make it comparable to observations
    """

    if BACKWARD_SAMPLING == 'MLE':
        for dilution in dilutions.transpose():
            # dilution = np.array(dilution)
            size = size/dilution
            # size = [x / y for x, y in zip(size, dilution)]
    else:
        sys.exit("Unsupported backward sampling method")

    if FORWARD_SAMPLING == 'MLE':
        for sample in samplings.transpose():
            # sample = np.array(sample)
            size = size*sample
            # size = [x * y for x, y in zip(size, sample)]
    elif FORWARD_SAMPLING == 'RV':
        for sample in samplings.transpose():
            # sample = np.array(sample)
            size = np.random.poisson(size*sample)
            # size = [np.random.poisson(x * y) for x, y in zip(size, sample)]
    else:
        sys.exit("Unsupported forward sampling method")

    return np.round(size)


def parse_observations(infile):
    """
    Parse csv of observations
    Creates a dictionary of dictionaries
    The inner ones hold single wells/colonies/whatever
    the outer ones holds that data coupled with their names
    """

    observed = {}

    with open(infile, 'r') as in_csv:
        rdr = csv.DictReader(in_csv)
        for obs in rdr:
            # id_string = '.'.join([obs[x] for x in ['name', 'generation', 'well']])
            id_string = obs['name']
            if id_string not in observed:
                observed[id_string] = {str(x): [] for x in rdr.fieldnames}
            for k in rdr.fieldnames:
                entry = obs[k]
                if entry:
                    if k in ['name', 'well']:
                        observed[id_string][k].append(str(entry))
                    elif k in ['generation', 'birthrate_group', 'deathrate_group']:
                        observed[id_string][k].append(int(entry))
                    else:
                        observed[id_string][k].append(float(entry))
                else:
                    if k in ['count', 'dead']:
                        observed[id_string][k].append(None)
                    else:
                        observed[id_string][k].append(1.0)

    #         for k, v in observed[id_string].items():
    #             if type(v[0]) in [int, float]:
    #                 observed[k] = np.array(v)

    # print(observed)

    global OBSERVED
    OBSERVED = copy.deepcopy(observed)
    return observed


def parse_params(paramfile, observed=None):
    """
    parse toml parameter file and observed data for lb-process parameters that are not the birthrate
    """
    # set model parameters
    global PARAMS
    PARAMS = toml.load(paramfile)

    # set defaults for variables where parameters are not mandatory
    if 'starting_cell_count' not in PARAMS['simulation_params']:
        PARAMS['simulation_params']['starting_cell_count'] = 'calculate'
    if 'end_time' not in PARAMS['simulation_params']:
        PARAMS['simulation_params']['end_time'] = 'max_observed'

    if 'min_starting_cell_count' not in PARAMS['simulation_params']:
        PARAMS['simulation_params']['min_starting_cell_count'] = 0
    if 'starting_population_size' not in PARAMS['abc_params']:
        PARAMS['abc_params']['starting_population_size'] = 100
    if 'min_epsilon' not in PARAMS['abc_params']:
        PARAMS['abc_params']['min_epsilon'] = 0.1
    if 'max_populations' not in PARAMS['abc_params']:
        PARAMS['abc_params']['max_populations'] = 10
    if 'min_acceptance' not in PARAMS['abc_params']:
        PARAMS['abc_params']['min_acceptance'] = 0.0
    if 'plot_params' not in PARAMS:
        PARAMS['plot_params'] = {}
        PARAMS['plot_params'][['population_measure']] = 'Cells'

    if 'distance_function' not in PARAMS['abc_params']:
        PARAMS['abc_params']['distance_function'] = 'linear'

    # if we specified carrying capacity
    if 'deathrate_interaction' not in PARAMS['simulation_params']:
        PARAMS['simulation_params']['deathrate_interaction'] = 1.0/PARAMS['simulation_params']['carrying_capacity']

    # finalize parsing
    if PARAMS['simulation_params']['starting_cell_count'] == 'calculate':
        if observed is None:
            sys.exit("Cannot compute starting cell count without observations")
        PARAMS['starting_population'] = {}
        for id_string, obs in observed.items():
            samplings, dilutions = get_samplings_dilutions(obs)
            samplings = samplings[0]
            dilutions = dilutions[0]
            pop = max(obs['count'][0], PARAMS['simulation_params']['min_starting_cell_count'])
            if FORWARD_SAMPLING == 'MLE' and BACKWARD_SAMPLING == 'MLE':
                for sample in samplings:
                    pop /= sample
                for dilution in dilutions:
                    pop *= dilution
                PARAMS['starting_population'][id_string] = lambda x=int(pop): x
            if FORWARD_SAMPLING == 'RV' and BACKWARD_SAMPLING == 'MLE':
                def f(x=pop, y=copy.deepcopy(samplings), z=copy.deepcopy(dilutions)):
                    # print(x, list(y), list(z))
                    for sample in y:
                        x /= sample
                    for dilution in z:
                        x = np.random.poisson(x * dilution)
                    return int(x)
                PARAMS['starting_population'][id_string] = f

    else:
        PARAMS['starting_population'] = lambda x=int(PARAMS['simulation_params']['starting_cell_count']): x
    # no need to simulate longer than observed segment
    PARAMS['end_time'] = {}
    if PARAMS['simulation_params']['end_time'] == 'max_observed':
        for id_string, obs in observed.items():
            if observed is None:
                sys.exit("Cannot compute end_time: 'max_observed' without observations")
            PARAMS['end_time'][id_string] = lambda x=max(obs['time']): x
    else:
        PARAMS['end_time'][id_string] = lambda x=float(PARAMS['end_time']): x

