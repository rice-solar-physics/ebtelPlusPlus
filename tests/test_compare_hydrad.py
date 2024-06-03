"""
Compare output of EBTEL IDL and ebtel++
"""
import pytest
from collections import OrderedDict

import astropy.units as u
from scipy.interpolate import interp1d

from .helpers import DATA_DIR, read_hydrad_test_data, run_ebtelplusplus


@pytest.fixture
def base_config():
    base_config = {
        'total_time': 5e3,
        'tau': 0.1,
        'tau_max': 10,
        'loop_length': 4e9,
        'saturation_limit': 1,
        'force_single_fluid': False,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'use_power_law_radiative_losses': True,
        'calculate_dem': False,
        'save_terms': False,
        'use_adaptive_solver': True,
        'radiation': 'power_law',
        'adaptive_solver_error': 1e-8,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 6.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.075,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 1.,
            'background': 3.5e-5,}),
    }
    return base_config


@pytest.mark.parametrize('tau', [200.0, 500.0])
@pytest.mark.parametrize('heating_type', ['single', 'electron', 'ion'])
def test_compare_hydrad_single_event_peak_values(base_config, tau, heating_type):
    config = base_config.copy()
    if heating_type == 'single':
        config['force_single_fluid'] = True
        config['heating']['partition'] = 0.5
    elif heating_type == 'electron':
        config['force_single_fluid'] = False
        config['heating']['partition'] = 1.0
    elif heating_type == 'ion':
        config['force_single_fluid'] = False
        config['heating']['partition'] = 0.0
    else:
        raise ValueError('Unrecognized heating type')
    total_energy = 10.0
    config['heating']['events'] = [
        {'event': {
            'magnitude': 2 * total_energy/tau,
            'rise_start': 0.,
            'rise_end': tau/2.,
            'decay_start': tau/2.,
            'decay_end': tau,
        }}
    ]
    r_ebtel = run_ebtelplusplus(config)
    r_hydrad = read_hydrad_test_data(DATA_DIR / 'hydrad_results.h5', tau, heating_type)
    # Require all quantities at the peak to be <20% in accordance with the comparisons
    # done in Barnes et al. (2016a)
    for name in ['electron_temperature', 'ion_temperature', 'density']:
        f_interp = interp1d(r_ebtel['time'].to_value('s'), r_ebtel[name].to_value(r_hydrad[name].unit), 
                            fill_value='extrapolate')
        x_ebtel_interp = u.Quantity(f_interp(r_hydrad['time'].to_value('s')), r_hydrad[name].unit)
        i_peak = r_hydrad[name].argmax()
        assert u.allclose(x_ebtel_interp[i_peak], r_hydrad[name][i_peak], rtol=0.20)
