import satellite
import tudat
import determine_orbit
import eclipse


## Determine the optimal orbits
N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu = determine_orbit.determine_orbit()

## Create the satellite objects with the tudat object as well
SATS = satellite.create_sats(N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu)
# tudat_env = tudat.set_up_and_run_tudat(SATS)

## eclipse calculations
eclipse.calc_eclipse_times(SATS)

## TODO: perturbation calculations


## Save satellite objects
for sat in SATS:
	sat.save_sat(sat.name)
