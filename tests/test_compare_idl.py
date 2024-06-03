"""
Compare output of EBTEL IDL and ebtel++
"""
import pytest
from collections import OrderedDict

import astropy.units as u

from .helpers import DATA_DIR, read_idl_test_data, run_ebtelplusplus


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
        'use_flux_limiting': False,
        'use_adaptive_solver': False,
        'calculate_dem': False,
        'radiation': 'power_law',
        'save_terms': False,
        'adaptive_solver_error': 1e-6,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 2.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.0,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 0.5,
            'background': 1e-6,
        }),
    }
    return base_config


@pytest.mark.xfail
def test_compare_idl_single_event(base_config, tmpdir, ebtel_idl_path):
    config = base_config.copy()
    config['heating']['events'] = [
        {'event': {
            'rise_start': 0.0, 'rise_end': 100.0, 'decay_start': 100.0,
            'decay_end': 200.0, 'magnitude': 0.1}}]
    r_cpp = run_ebtelplusplus(config)
    r_idl = read_idl_test_data(DATA_DIR / 'data' / 'idl_single_event.json',
                               ebtel_idl_path, config)
    assert u.allclose(r_cpp['electron_temperature'], r_idl['temperature'], rtol=1e-2)
    assert u.allclose(r_cpp['ion_temperature'], r_idl['temperature'], rtol=1e-2)
    assert u.allclose(r_cpp['density'], r_idl['density'], rtol=1e-2)
    assert u.allclose(r_cpp['electron_pressure']+r_cpp['ion_pressure'], r_idl['pressure'], 
                      rtol=1e-2)
    assert u.allclose(r_cpp['velocity'], r_idl['velocity'], rtol=1e-2)
