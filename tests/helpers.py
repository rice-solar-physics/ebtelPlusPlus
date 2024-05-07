"""
Helper functions for tests
"""
import os
import pathlib
import sys

import astropy.units as u
import h5py
import numpy as np

TEST_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = TEST_DIR / 'data'
REPO_DIR = TEST_DIR.parent
sys.path.append(str(REPO_DIR / 'examples'))
from util import run_ebtel


def run_ebtelplusplus(config):
    res = run_ebtel(config, REPO_DIR)
    return res


def generate_idl_test_data(ebtel_idl_path, config):
    # Import here to avoid making this a hard dependency for this whole
    # script.
    import hissw
    flags = []
    if 'dem' not in config or not config['dem']['use_new_method']:
        flags += ['dem_old']
    if not config['use_flux_limiting']:
        flags += ['classical']
    # Set expansions keywords
    # NOTE: It is assumed that A_TR is always 1 here which we can do wlog
    keywords = {
        'a_tr': 1.0,
        'a_c': 1/config['area_ratio_tr_corona'],
        'a_0': config['area_ratio_0_corona']/config['area_ratio_tr_corona'],
        'l_fact': 1.0 - config['loop_length_ratio_tr_total'],
    }
    # keywords = [f'{k}={v:.6f}' for k,v in keywords.items()]
    time = np.arange(0, config['total_time']-config['tau'], config['tau'])
    heat = np.ones(time.shape) * config['heating']['background']
    for _e in config['heating']['events']:
        e = _e['event']
        # Rise
        i = np.where(np.logical_and(time >= e['rise_start'], time < e['rise_end']))
        heat[i] += e['magnitude'] * (time[i] - e['rise_start']) / (e['rise_end'] - e['rise_start'])
        # Plateau
        i = np.where(np.logical_and(time >= e['rise_end'], time < e['decay_start']))
        heat[i] += e['magnitude']
        # Decay
        i = np.where(np.logical_and(time >= e['decay_start'], time <= e['decay_end']))
        heat[i] += e['magnitude'] * (e['decay_end'] - time[i])/(e['decay_end'] - e['decay_start'])

    args = {
        'time': time.tolist(),
        'loop_length': config['loop_length'],
        'heat': heat.tolist(),
        'flags': flags,
        **keywords,
    }
    idl = hissw.Environment(extra_paths=[ebtel_idl_path])
    script = """time = {{ time | force_double_precision }}
heat = {{ heat | force_double_precision }}
loop_length = {{ loop_length | force_double_precision }}
ebtel2,time,heat,loop_length,temperature,density,pressure,velocity,$
    a_tr={{ a_tr | force_double_precision }},$
    a_c={{ a_c | force_double_precision }},$
    a_0={{ a_0 | force_double_precision }},$
    l_fact={{ l_fact | force_double_precision }}{% if flags %}, /{{ flags | join(', /') }}{% endif %}
    """
    return idl.run(script, args=args)


def read_idl_test_data(data_filename, ebtel_idl_path, config):
    varnames = ['time', 'temperature', 'density', 'pressure', 'velocity']
    varunits = ['s', 'K', 'cm-3', 'dyne cm-2', 'cm s-1']
    # Generate and save if it does not exist
    if ebtel_idl_path is not None: #not os.path.isfile(data_filename):
        data = generate_idl_test_data(ebtel_idl_path, config)
        data_array = np.zeros(data['time'].shape+(len(data),))
        for i, v in enumerate(varnames):
            data_array[:, i] = data[v]
        np.savetxt(data_filename, data_array)
    # Load data into a dictionary
    data = np.loadtxt(data_filename)
    return {v: u.Quantity(data[:, i], vu) for i, (v,vu) in enumerate(zip(varnames, varunits))}


def read_hydrad_test_data(data_filename, tau, heating):
    data = {}
    with h5py.File(data_filename, 'r') as hf:
        grp = hf[f'/{heating}/tau{tau:.0f}']
        data['time'] = u.Quantity(hf['time'], hf['time'].attrs['unit'])
        data['electron_temperature'] = u.Quantity(grp['electron_temperature'],
                                                  grp['electron_temperature'].attrs['unit'])
        data['ion_temperature'] = u.Quantity(grp['ion_temperature'],
                                             grp['ion_temperature'].attrs['unit'])
        data['density'] = u.Quantity(grp['density'],
                                     grp['density'].attrs['unit'])
    return data


def plot_comparison(r_cpp, r_idl):
    import matplotlib.pyplot as plt
    plt.subplot(121)
    plt.plot(r_cpp['time'], r_cpp['electron_temperature'], label='ebtel++')
    plt.plot(r_idl['time'], r_idl['temperature'], label='IDL')
    plt.legend()
    plt.subplot(122)
    plt.plot(r_cpp['time'], r_cpp['density'])
    plt.plot(r_idl['time'], r_idl['density'])
    plt.show()
