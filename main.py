import satellite
import determine_orbit
import eclipse
import propulsion


## Determine the optimal orbits
Orbital_parameters = determine_orbit.determine_orbit()

## Create the satellite objects with the tudat object as well
SATS = satellite.create_sats(Orbital_parameters)

## eclipse calculations
eclipse.calc_eclipse_times(SATS)

## TODO: perturbation calculations


## Propulsion calculations
# Preliminary values required
SATS[0].DeltaH_desaturation = 3 * 31  # [Nms]
SATS[0].DeltaV_maintenance = 3.75  # [m/s]
SATS[0].freq_desaturation = 1 / 2  # [/day]
SATS[0].freq_maintenance = 1 / 84  # [/day]
propulsion.determine_propellant_mass(SATS[0].dry_mass, SATS[0].lifetime, SATS[0])

## Save satellite objects
for sat in SATS:
	sat.save_sat(sat.name)
