An ebtel++ run is configured by a single XML configuration file. The table below lists all of the input parameters needed for a run. An [example configuration file](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/config/ebtel.example.cfg.xml) is included in the GitHub repository. The configuration file can be written manually or configured via Python as is done in the [included examples](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples).

| Parameter | Type | Description |
|:---------:|:----:|:------------|
| **total_time** | `float` | duration of the simulation (in s)|
| **tau** | `float` | timestep (in s); if using adaptive solver, the initial timestep |
| **tau_max** | `float` | maximum allowed timestep (in s) when using adaptive solver |
| **loop_length** | `float` | Loop half-length (in cm) |
| **saturation_limit** | `float` | Flux limiter, $f$ in section 2.1 of [Barnes et al. (2016)][barnes_2016] |
| **force_single_fluid** | `bool` | if True, electron and ion populations forced into equilibrium |
| **use_c1_loss_correction** | `bool` | use correction in Eq. 16 of [Cargill et al. (2012a)][cargill_2012a] |
| **use_c1_grav_correction** | `bool` | use correction in Eq. 12 of [Cargill et al. (2012a)][cargill_2012a] |
| **use_flux_limiting** | `bool` | impose a flux limiter according to Eq. 22 of [Klimchuk et al. (2008)][klimchuk_2008] |
| **calculate_dem** | `bool` | if True, do the TR and coronal DEM calculation; increases compute time significantly |
| **save_terms** | `bool` | if True, save heat flux, $c_1$ parameter, and radiative loss to a separate file `<output_filename>.terms` |
| **use_adaptive_solver** | `bool` | if True, use adaptive timestep; significantly smaller compute times. In both cases, a Runge-Kutta Cash-Karp integration method is used (see section 16.2 of [Press et al. (1992)][press_num_recipes])  |
| **output_filename** | `string` | path to output file |
| **adaptive_solver_error** | `float` | Allowed truncation error in adaptive timestep routine |
| **adaptive_solver_safety** | `float` | Refinement factor, between 0 and 1, used if timestep becomes too large and solution contains NaNs. Especially important for short, infrequently heated loops. Also controls decreases in timestep due to thermal conduction timestep. Suggested value is 0.5 |
| **c1_cond0** | `float` | Nominal value of $c_1$ during the conduction phase; see Appendix A of [Barnes et al. (2016)][barnes_2016] |
| **c1_rad0** | `float` | Nominal value of $c_1$ during radiative phase; see Eq. 16 of [Cargill et al. (2012a)][cargill_2012a] |
| **helium_to_hydrogen_ratio** | `float` | Ratio of helium to hydrogen abundance; used in correction to ion mass, ion equation of state |
| **surface_gravity** | `float` | Surface gravity in units of solar surface gravity; should be set to 1.0 unless using for extra-solar cases |
| **loop_length_ratio_tr_total** | `float` | Ratio between the length of the TR and the total loop length. Typically, a value of 0.15 is used. |
| **area_ratio_tr_corona** | `float` | Ratio between the cross-sectional area averaged over the transition region and averaged over the corona |
| **area_ratio_0_corona** | `float` | Ratio between the cross-sectional area at the TR-corona boundary and the cross-sectional area averaged over the corona |

### Heating
The time dependent heating is configured in a separate node. It includes the following parameters,

| Parameter | Type | Description |
|:-------:|:------:|:-----------|
| **partition** | `float` | partition of heating between electrons and ions, between 0 and 1; 1 is pure electron heating, 0 pure ion heating |
| **background** | `float` | constant background heating (in ergs cm$^{-3}$ s$^{-1}$) |

The heating function is constructed by a list of discrete events and should be specified in the following way,

```XML
<events>
  <event magnitude="0.1" rise_start="0.0" rise_end="50.0" decay_start="50.0" decay_end="100.0"/>
  <event magnitude="0.05" rise_start="1000.0" rise_end="1250.0" decay_start="1350.0" decay_end="1450.0"/>
</events>
```

Here, we've configured two separate heating events. The first starts at 0 seconds, rises linearly to a maximum heating rate of 0.1 erg cm$^{-3}$ s$^{-1}$ in 50 seconds, and then immediately falls off with the event concluding at 100 seconds, i.e. a triangular heating profile. The second starts at 1000 seconds, rises to a maximum heating rate of 0.05 erg cm$^{-3}$ s$^{-1}$ in 250 seconds, is sustained at 0.05 for 100 seconds and then the event concludes at 1450 seconds.

Using this format, it is easy to specify either symmetric or asymmetric events of many different shapes. For more examples, see the [example configuration file](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/config/ebtel.example.cfg.xml) or the included [examples](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples).

### Differential Emission Measure
Optionally, ebtel++ can can also calculate the differential emission measure (DEM) in both the transition region and the corona. See sections 2.2 and 3 of [Klimchuk et al. (2008)][klimchuk_2008] for the details of this calculation. To enable this calculation, set `calculate_dem` to `True` in the configuration file (as described above). Note that this will result in much longer computation times.

If `calculate_dem` is set to `True`, ebtel++ will read in the `dem` node from the configuration file which should be structured as follows,
```XML
<dem>
  <use_new_method>True</use_new_method>
  <temperature bins="451" log_min="4" log_max="8.5"/>
</dem>
```

If `use_new_method` is set to True (False), the transition region DEM is calculated using the method outlined in section 3 (the appendix) of [Klimchuk et al. (2008)][klimchuk_2008]. The `temperature` node configures the range and number of bins used when calculating the DEM. Here, for example, there are 450 bins of equal width between $10^4$ and $10^{8.5}$ K.

If you do not need to calculate the DEM, set the `calculate_dem` parameter to False and this section of the configuration file need not be included.

## Output
Once the EBTEL run has finished, the results are printed to the file specified in `output_filename` in the configuration file (as described above). Several examples of how to parse the results in Python can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples). In general, the results file follows the structure,

| | | | | | | | |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| $t_0$ | $T_e(t_0)$ | $T_i(t_0)$ | $n(t_0)$ | $p_e(t_0)$ | $p_i(t_0)$ | $v(t_0)$ | $h(t_0)$ |
| ... | ... | ... | ... | ... | ... | ... | ... |
| $t_i$ | $T_e(t_i)$ | $T_i(t_i)$ | $n(t_i)$ | $p_e(t_i)$ | $p_i(t_i)$ | $v(t_i)$ | $h(t_i)$ |
| ... | ... | ... | ... | ... | ... | ... | ... |
| $t_{N-1}$ | $T_e(t_{N-1})$ | $T_i(t_{N-1})$ | $n(t_{N-1})$ | $p_e(t_{N-1})$ | $p_i(t_{N-1})$ | $v(t_{N-1})$ | $h(t_{N-1})$ |


Here $t$ is the time, $T_e$ is the electron temperature, $T_i$ is the ion temperature, $p_e$ is the electron pressure, $p_i$ is the ion pressure, $n$ is the density, $v$ is the velocity, and $h$ is the heating rate.

If `calculate_dem` is set to True, the TR and coronal DEM results are printed to `<output_filename>.dem_tr` and `<output_filename>.dem_corona`, respectively. These output files are structured in the following way,

| | | | | |
|:----:|:----:|:----:|:----:|:----:|
| $T_0$ | ... | $T_j$ | ... | $T_{M-1}$ |
| DEM$(t_0,T_0)$ | ... | DEM$(t_0,T_j)$ | ... | DEM$(t_0,T_{M-1})$ |
| ... | ... | ... | ... | ... |
| DEM$(t_i,T_0)$ | ... | DEM$(t_i,T_j)$ | ... | DEM$(t_i,T_{M-1})$ |
| ... | ... | ... | ... | ... |
| DEM$(t_{N-1},T_0)$ | ... | DEM$(t_{N-1},T_j)$ | ... | DEM$(t_{N-1},T_{M-1})$ |

where $M$ is the number of temperature bins and $N$ is again the number of timesteps.

[klimchuk_2008]: http://adsabs.harvard.edu/abs/2008ApJ...682.1351K "Klimchuk et al. (2008)"
[cargill_2012a]: http://adsabs.harvard.edu/abs/2012ApJ...752..161C "Cargill et al. (2012a)"
[cargill_2012b]: http://adsabs.harvard.edu/abs/2012ApJ...758....5C "Cargill et al. (2012b)"
[barnes_2016]: http://adsabs.harvard.edu/abs/2016ApJ...829...31B "Barnes et al. (2016)"
[press_num_recipes]: http://dl.acm.org/citation.cfm?id=148286 "Press et al. (1992)"
