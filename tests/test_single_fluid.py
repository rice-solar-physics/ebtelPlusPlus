"""
Test whether single-fluid option keeps electron and ion temperatures equal
"""
import os
from collections import OrderedDict

import astropy.units as u
import pytest

from .helpers import run_ebtelplusplus


@pytest.fixture
def base_config():
    base_config = {
        'total_time': 5e3,
        'tau': 1.0,
        'tau_max': 10.0,
        'loop_length': 4e9,
        'saturation_limit': 1/6,
        'force_single_fluid': True,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'calculate_dem': False,
        'save_terms': False,
        'use_adaptive_solver': True,
        'radiation': 'power_law',
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


@pytest.mark.parametrize(['rise_start', 'rise_end', 'decay_start', 'decay_end', 'magnitude'], [
    (0, 250, 1000, 2000, 0.005),  # gentle case
    (0, 100, 100, 200, 0.1),  # more impulsive case
])
def test_single_fluid_gentle(base_config, rise_start, rise_end, decay_start, decay_end, magnitude):
    base_config['heating']['events'] = [
        {'event': {
            'rise_start': rise_start,
            'rise_end': rise_end,
            'decay_start': decay_start,
            'decay_end': decay_end,
            'magnitude': magnitude,
        }}
    ]
    results = run_ebtelplusplus(base_config)
    assert u.allclose(results['electron_temperature'], results['ion_temperature'], rtol=1e-10)
    assert u.allclose(results['electron_pressure'], results['ion_pressure'], rtol=1e-10)
