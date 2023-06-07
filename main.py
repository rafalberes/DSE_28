import numpy as np
import constants as const
import satellite
import tudat
import determine_orbit
import eclipse

## Set up the tudat environment
tudat_env = tudat.set_up_tudat_environment()

## Determine the optimal orbits
N_orbits, N_sat, j, k, DeltaL, a, e, i, r_p, r_a, T, OMEGA, omega, nu = determine_orbit.determine_orbit()

## Create the satellite objects with the tudat object as well
SATS = satellite.create_sats(N_orbits, N_sat, j, k, DeltaL, a, e, i, r_p, r_a, T, OMEGA, omega, nu)

## eclipse calculations
eclipse.calc_eclipse_times(SATS)

## TODO: perturbation calculations


## Save satellite objects
for sat in SATS:
	sat.save_sat(sat.name)
