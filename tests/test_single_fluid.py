"""
Test whether single-fluid option keeps electron and ion temperatures equal
"""
import os
from collections import OrderedDict

import pytest
import numpy as np

from .helpers import run_ebtelplusplus


@pytest.fixture
def base_config():
    base_config = {
        'total_time': 5e3,
        'tau': 1.0,
        'tau_max': 10.0,
        'loop_length': 4e9,
        'loop_length_ratio_tr_total': 0.0,
        'area_ratio_tr_corona': 1.0,
        'area_ratio_0_corona': 1.0,
        'saturation_limit': 1/6,
        'force_single_fluid': True,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'calculate_dem': False,
        'save_terms': False,
        'use_adaptive_solver': True,
        'adaptive_solver_error': 1e-6,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 2.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.075,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 0.5,
            'background': 1e-6,
        }),
    }
    return base_config


def test_single_fluid_gentle(base_config, tmpdir):
    base_config['heating']['events'] = [
        {'event': {'rise_start': 0.0, 'rise_end': 250.0, 'decay_start': 1000.0,
                   'decay_end': 2000.0, 'magnitude': 0.005}}]
    results = run_ebtelplusplus(base_config)
    # Temperature
    assert np.allclose(results['electron_temperature'],
                       results['ion_temperature'], atol=0., rtol=1e-10)
    # Pressure
    assert np.allclose(results['electron_pressure'], results['ion_pressure'],
                       atol=0., rtol=1e-10)


def test_single_fluid_impulsive(base_config, tmpdir):
    base_config['heating']['events'] = [
        {'event': {'rise_start': 0.0, 'rise_end': 100.0, 'decay_start': 100.0,
                   'decay_end': 200.0, 'magnitude': 0.1}}]
    results = run_ebtelplusplus(base_config)
    # Temperature
    assert np.allclose(results['electron_temperature'],
                       results['ion_temperature'], atol=0., rtol=1e-10)
    # Pressure
    assert np.allclose(results['electron_pressure'], results['ion_pressure'],
                       atol=0., rtol=1e-10)
