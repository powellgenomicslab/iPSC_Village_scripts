
"""
Functions for reconstructing a reproduction rate timeline
using approximate bayesian computing (ABC)
and a logistic branching process (lb-process).
"""

import copy
import re
import sys

import click
import numpy as np
from pyabc import (ABCSMC, Distribution, RV)
from pyabc.sampler import MulticoreEvalParallelSampler
from pyabc.populationstrategy import AdaptivePopulationSize
from pyabc.populationstrategy import ConstantPopulationSize
from pyabc import History
import toml

import simtools


@click.group()
def main():
    pass


def flatten_observed(observed):
    flat = {}
    for id_string, obs in observed.items():
        for k, v in obs.items():
            flat[id_string + '.' + k] = v
    return flat


def rmsd(v1, v2):
    """
    simple rmsd between two vectors of equal length
    """
    assert v1.shape == v2.shape
    return np.sqrt(np.sum((v1 - v2)**2.0))
    # return math.sqrt(sum([(x1 - x2)**2.0 for x1, x2 in zip(v1, v2)]))


def select_time_match(time, *obs):
    """
    Find timepoints from simulation that best match timepoints in observations
    """
    # flatten observations into single sorted list
    obs_time = np.array(sorted([x for y in obs for x in y]))
    # set error limit (grows with number of observations)
    error = 1.0 * len(obs_time)
    selected = []
    matching = []
    iprev = 0
    for j, t_obs in enumerate(obs_time):
        diffmin = 99999
        i = 0
        for i, t_sim in enumerate(time[iprev:]):
            if abs(t_sim - t_obs) < diffmin:
                diffmin = abs(t_sim - t_obs)
                imin = iprev + i
            else:
                break
        selected.append(imin)
        matching.append(j)
        iprev += i - 1
    # find rmsd if x fit is good enough
    assert rmsd(obs_time, np.array([time[x] for x in selected])) < error
    # if it's good enough, we can use it
    return np.array(selected), np.array(matching)


def distance(simulation, observation):
    """
    rmsd between a simulated growth curve and a set of experimental datapoints
    """
    distances = []
    # print(simulation, observation)
    for id_string in simtools.OBSERVED:
        sim = {str(k): simulation[id_string + '.' + str(k)] for k in ['time', 'size', 'rate']}
        obs = {str(k): observation[id_string + '.' + str(k)] for k in simtools.OBSERVED[id_string]}
        # print(sim, obs)
        # selected_size = [sim['size'][x] for x in sim_selected]
        selected_size = sim['size']
        filters = copy.deepcopy(simtools.PARAMS['filters'])
        samplings, dilutions = simtools.get_samplings_dilutions(obs)
        # print(samplings, dilutions)
        # apply noise filters
        selected_count = simtools.apply_sampling(selected_size, samplings, dilutions)
        # print(type(selected_count))
        selected_count = simtools.apply_noise(selected_count, filters)
        # print(selected_count)
        if simtools.PARAMS['abc_params']['distance_function'] == 'linear':
            distances.append(np.sum(np.abs(np.array(obs['count']) - selected_count)))
        elif simtools.PARAMS['abc_params']['distance_function'] == 'rmsd':
            distances.append(rmsd(np.array(obs['count']), selected_count))
    return sum(distances)


def abc_model(params):
    """
    model for abc computation
    run one timeline for each observation
    """
    data = {}
    for id_string, obs in simtools.OBSERVED.items():
        re_birthrates = re.compile(r'r([0-9])')
        kvs = sorted([(k, v) for k, v in params.items() if re_birthrates.search(k)],
                     key=lambda x: int(re_birthrates.search(x[0]).group(1)))
        birthrate = [x[1] for x in kvs]
        deathrate_interaction = simtools.PARAMS['simulation_params']['deathrate_interaction']
        # print('obs', obs)
        time, size, rate = simtools.simulate_timeline(
            simtools.PARAMS['starting_population'][id_string](),
            obs['time'],
            birthrate,
            deathrate_interaction,
            simtools.PARAMS['abc_params']['simulator'],
            # verbosity=1
        )
        data[id_string] = {
            'time': time,
            'size': size,
            'rate': rate,
        }
    data = flatten_observed(data)
    data['simulation'] = True # tag as simulation data for distance calculation
    return data


def abc_distance(a, b):
    """
    Distance for abc computation. Simply runs distance using abc model dicts
    """
    # print(a, b)
    if 'simulation' in a:
        return distance(a, b)
    else:
        return distance(b, a)


def abc_setup(birthrate_groups):
    """
    create abc model
    parameters are stored in the global simtools.PARAMS dict
    """
    for curve_resolution in simtools.PARAMS['abc_params']['resolution_limits']:
        assert curve_resolution > 0 and curve_resolution <= 9
    abc_priors = []
    for resolution_limit in range(
            simtools.PARAMS['abc_params']['resolution_limits'][0],
            simtools.PARAMS['abc_params']['resolution_limits'][1] + 1):
        abc_prior_dict = {}
        for i in range(resolution_limit):
            abc_prior_dict['r' + str(i)] = \
                RV("uniform", simtools.PARAMS['abc_params']['rate_limits'][0],
                abs(simtools.PARAMS['abc_params']['rate_limits'][1] - \
                simtools.PARAMS['abc_params']['rate_limits'][0]))
        abc_priors.append(Distribution(birthrate=copy.deepcopy(abc_prior_dict)))
    print('priors', abc_priors)
    #abc = ABCSMC([abc_model for __ in abc_priors], abc_priors, abc_distance,
    #             population_size=AdaptivePopulationSize(
    #                 int(simtools.PARAMS['abc_params']['starting_population_size']),
    #                 0.15,
    #                 max_population_size=int(simtools.PARAMS['abc_params']['max_population_size']),
    #                 min_population_size=int(simtools.PARAMS['abc_params']['min_population_size'])),
    #             sampler=MulticoreEvalParallelSampler(
    #                 simtools.PARAMS['abc_params']['parallel_simulations']))
    abc = ABCSMC([abc_model for __ in abc_priors], abc_priors, abc_distance,
                 population_size=ConstantPopulationSize(
                     int(simtools.PARAMS['abc_params']['starting_population_size'])),
                 sampler=MulticoreEvalParallelSampler(
                     simtools.PARAMS['abc_params']['parallel_simulations']))
    return abc



@main.command()
@click.option('-p', '--paramfile', type=click.Path())
@click.option('-o', '--obsfile', type=click.Path())
@click.option('-d', '--dbfile', type=click.Path())
def reconstruct(paramfile, obsfile, dbfile):
    """
    Reconstruct a likely reproduction rate function given a set of experimental observations
    """

    # Generate observed dictionary from input observations
    observed = simtools.parse_observations(obsfile)
    print('Observed data:', observed)

    # Set simulation parameters
    simtools.parse_params(paramfile, observed)
    # print(simtools.PARAMS)
    print('Starting populations (poisson distributed)')
    for k, v in simtools.PARAMS['starting_population'].items():
        print(k, v())
    print('Simulation end times')
    for k, v in simtools.PARAMS['end_time'].items():
        print(k, v())
    print('Using simulator:', simtools.PARAMS['abc_params']['simulator'])

    observed = flatten_observed(observed)
    print('Observed data (flat):', observed)

    # generate abc model
    abc = abc_setup(len({v[0] for k, v in observed.items() if 'birthrate_group' in k}))

    db_path = 'sqlite:///' + dbfile
    print('Saving database in:', db_path, file=sys.stderr)

    # run abc
    print('Constructing ABC', file=sys.stderr)
    abc.new(db_path, observed)
    print('Running ABC', file=sys.stderr)
    abc.run(minimum_epsilon=simtools.PARAMS['abc_params']['min_epsilon'],
            max_nr_populations=simtools.PARAMS['abc_params']['max_populations'],
            min_acceptance_rate=simtools.PARAMS['abc_params']['min_acceptance'])


if __name__ == '__main__':
    main()
