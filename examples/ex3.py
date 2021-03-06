"""
Third example: single-fluid, single asymmetric pulse
"""

import os

import numpy as np

from plot import make_figure
from util import run_ebtel, read_xml


top_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


if __name__ == '__main__':
    # Configure and run ebtel++
    base = read_xml(
        (os.path.join(top_dir, 'config', 'ebtel.example.cfg.xml')))
    base['calculate_dem'] = True
    base['heating']['partition'] = 0.5
    base['force_single_fluid'] = True
    base['heating']['events'] = [{'event': {
        'rise_start': 0.0,
        'rise_end': 250.0,
        'decay_start': 1000.0,
        'decay_end': 2000.0,
        'magnitude': 0.005
    }}]
    results = run_ebtel(base, top_dir)
    # Time average DEM
    results['dem_total'] = np.average(
        results['dem_tr']+results['dem_corona'],
        axis=0,
        weights=np.gradient(results['time']))
    results['dem_tr'] = np.average(
        results['dem_tr'],
        axis=0,
        weights=np.gradient(results['time']))
    results['dem_corona'] = np.average(
        results['dem_corona'],
        axis=0,
        weights=np.gradient(results['time']))
    # Plot
    figure_filename = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'ex3.png')
    make_figure(results, figure_filename)
