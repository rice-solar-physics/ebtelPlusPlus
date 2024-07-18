"""
Compare results of adaptive and static solvers
"""
from collections import OrderedDict

import pytest
import numpy as np

from .helpers import run_ebtelplusplus
from util import EbtelPlusPlusError

@pytest.fixture(scope='module')
def base_config():
    base_config = {
        'total_time': 5e3,
        'tau': 1.0,
        'tau_max': 10.0,
        'loop_length': 4e9,
        'saturation_limit': 1/6,
        'force_single_fluid': False,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'calculate_dem': True,
        'save_terms': False,
        'radiation': 'power_law',
        'adaptive_solver_error': 1e-6,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 2.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.075,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 1.0,
            'background': 1e-6,
            'events': [
                {'event': {'rise_start': 0.0, 'rise_end': 100.0, 'decay_start': 100.0,
                           'decay_end': 200.0, 'magnitude': 0.1}}
            ]
        }),
        'dem': OrderedDict({
            'use_new_method': True,
            'temperature': {'bins': 451, 'log_min': 4, 'log_max': 8.5},
        }),
    }
    return base_config


@pytest.fixture(scope='module')
def adaptive_results(base_config):
    config = base_config.copy()
    config['use_adaptive_solver'] = True
    return run_ebtelplusplus(config)    


@pytest.fixture(scope='module')
def static_results(base_config):
    config = base_config.copy()
    config['use_adaptive_solver'] = False
    return run_ebtelplusplus(config)
    
@pytest.mark.parametrize(['name', 'atol'], [
    ('electron_temperature', 1e4),
    ('ion_temperature', 1e4),
    ('density', 1e7),
    ('electron_pressure', 1e-4),
    ('ion_pressure', 1e-4),
    ('velocity', 1e4),
])
def test_quantities_equal_adaptive_static(adaptive_results, static_results, name, atol):
    t_static = static_results['time'].to_value('s')
    t_adaptive = adaptive_results['time'].to_value('s')
    adapt_interp = np.interp(t_static, t_adaptive, adaptive_results[name].value)
    # NOTE: Skip the first 5 steps b/c there is always one anomalous point that gives
    # an error > 10%; due to static case not rising fast enough
    assert np.allclose(adapt_interp[5:], static_results[name][5:].to_value(adaptive_results[name].unit),
                       rtol=1e-2, atol=atol)

@pytest.mark.parametrize('value', [-1e-5, 0, 1e-15])
def test_insufficient_heating(base_config, value):
    config = base_config.copy()
    config['use_adaptive_solver'] = False
    config['heating']['background'] = value
    with pytest.raises(EbtelPlusPlusError):
        run_ebtelplusplus(config)

@pytest.mark.parametrize('use_adaptive_solver', [True, False])
def test_NaNs_in_solver(base_config, use_adaptive_solver):
    config = base_config.copy()
    config['use_adaptive_solver'] = use_adaptive_solver
    config['heating']['events'] = [
                {'event': {'rise_start': 0.0, 'rise_end': 100.0, 'decay_start': 100.0,
                           'decay_end': 200.0, 'magnitude': -10.0}}
            ]
    with pytest.raises(EbtelPlusPlusError):
        run_ebtelplusplus(config)

