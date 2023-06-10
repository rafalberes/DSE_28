from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import constants as const
import satellite
import tudat


def determine_orbit():
	# We have selected a three type orbit. From this determine the maximum available altitude to achieve this.
	# Knowing the j is 14, from trial and error
	
	N_orbits = 3
	
	j = 14
	k = 1
	N_sat = ceil(k * N_orbits / const.Temporal_Res)
	DeltaL = 2 * np.pi * k / j
	
	lat_measured = np.rad2deg(DeltaL / N_orbits) / const.m_to_deg
	
	accuracy = 25_000
	r_p_alt = (.5 * lat_measured / np.tan(.5 * const.FoV_PL)) // accuracy * accuracy + accuracy
	
	r_p = r_p_alt + const.R_E
	
	T = const.T_E * k / j
	a = ((T / (2 * np.pi)) ** 2 * const.mu_E) ** (1 / 3)
	
	Part1 = ((-1) * k / j * 2 * np.pi + 2 * np.pi * 2 * np.pi * np.sqrt(a ** 3 / const.mu_E) / const.T_E) * \
	        (-a ** 2 * (1 - (1 - r_p / a) ** 2) ** 2 / (3 * np.pi * const.J_2 * const.R_E ** 2))
	
	if np.abs(Part1) < 1e-10:
		i = np.pi / 2
	else:
		i = np.arccos(Part1)
	
	OMEGA = np.arange(0, DeltaL, DeltaL / N_orbits)
	
	## Argument of peri-centre
	lat = np.deg2rad(90)
	# Sketchy source: https://physics.stackexchange.com/questions/47094/convert-latitude-of-lowest-altitude-to-argument-of-perigee
	omega = np.arcsin((np.sin(lat)) / (np.sin(i)))
	
	e = 1 - r_p / a
	r_a = a * (1 + e)
	
	nu = np.array([0, 0, 0])
	
	DeltaL1 = -4 * np.pi ** 2 * np.sqrt(a ** 3 / const.mu_E) / const.T_E
	DeltaL2 = - (3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(i)) / (a ** 2 * (1 - e ** 2) ** 2)
	
	if __name__ == "__main__":
		print_orbital_parameters(N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu)
	
	DeltaV_insertion = calc_DeltaV_insertion(e, a, np.pi, DeltaL/N_orbits)
	DeltaV_deorbit = 179.108
	
	Orbital_parameters = {
		"N_orbits": N_orbits,
		"N_sat": N_sat,
		"j": j,
		"k": k,
		"DeltaL": DeltaL,
		"DeltaL1": DeltaL1,
		"DeltaL2": DeltaL2,
		"a": a,
		"e": e,
		"i": i,
		"r_p": r_p,
		"r_a": r_a,
		"T": T,
		"OMEGA": OMEGA,
		"omega": omega,
		"nu": nu,
		"DeltaV_insertion": DeltaV_insertion,
		"DeltaV_deorbit":DeltaV_deorbit,
	}
	
	return Orbital_parameters


def calc_DeltaV_insertion(e: float, a: float, nu: float = np.deg2rad(180), Delta_OMEGA: float = 8.57):
	v = np.sqrt(const.mu_E * ((2 + 2 * e * np.cos(nu)) / (a * (1 - e ** 2)) - 1 / a))
	DeltaV_simpl = 2 * v * np.sin(Delta_OMEGA / 2)

	if __name__ == "__main__":
		print(
			f"###### DeltaV required for insertion ######\n"
			f"Change of inclination       : {np.rad2deg(Delta_OMEGA)} [deg]\n"
			f"Velocity at apogee:         : {v} [m/s]\n"
			f"Simple calculations Delta V : {DeltaV_simpl} [m/s]\n"
		)
	
	return DeltaV_simpl


def print_orbital_parameters(N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu):
	r_40_lat = (a * (1 - e ** 2)) / (1 + e * np.cos(np.deg2rad(50)))
	r = 800e3 + const.R_E
	theta_r = np.arccos((a * (1 - e ** 2) / r - 1) / e)
	
	print(
		f"###### Orbital parameters ######\n"
		f"{N_orbits} types of orbits are used that will repeat in {k} days with {N_sat} satellites.\n"
		f"In the {k} days, {j} orbits are completed.\n"
		f"These orbits have:\n"
		f"DeltaL                          : {np.rad2deg(DeltaL)} [Deg]\n"
		f"DeltaL1                         : {np.rad2deg(DeltaL1)} [Deg]\n"
		f"DeltaL2                         : {np.rad2deg(DeltaL2)} [Deg]\n"
		f"semi-major axis, a              : {a} [m]\n"
		f"Eccentricity, e                 : {e} [-]\n"
		f"Inclination, i                  : {np.rad2deg(i)} [Deg]\n"
		f"Peri-centre radius, r_p         : {r_p} [m]\n"
		f"                                  {r_p - const.R_E} [m]\n"
		f"Apo-centre radius, r_a          : {r_a} [m]\n"
		f"                                  {r_a - const.R_E} [m]\n"
		f"Orbital period, T               : {T} [s]\n"
		f"                                  {T / 60} [min]\n"
		f"Ascending nodes, OMEGA          : {np.rad2deg(OMEGA)} [Deg]\n"
		f"Argument of peri-centre, omega  : {np.rad2deg(omega)} [Deg]\n"
		f"True anomaly, nu                : {np.rad2deg(nu)} [Deg]\n"
		f"\n"
		f"Altitude at 40 deg latitude     : {r_40_lat - const.R_E} [m]\n"
		f"True anomaly where alt is 800 km: {np.rad2deg(theta_r)} [deg]\n"
		f"Latitude where alt is 800 km    : {90 - np.rad2deg(theta_r)} [deg]\n"
	)


def coverage(tudat_env, sats):
	subset = tudat_env.create_plot(hours=12)  # sats[0].T*7/3600)
	for sat in sats:
		print(sat)
		
		if subset is None:
			latitude = np.rad2deg(tudat_env.dep_vars_array[:, 10])[::10]
			longitude = np.rad2deg(tudat_env.dep_vars_array[:, 11])[::10] + np.rad2deg(sat.OMEGA)
		else:
			latitude = np.rad2deg(tudat_env.dep_vars_array[:, 10])[:subset:10]
			longitude = np.rad2deg(tudat_env.dep_vars_array[:, 11])[:subset:10] + np.rad2deg(sat.OMEGA)
		
		FoV_lat = 8.68 / 2
		FoV_lon = 8.68 / 2
		print(len(latitude))
		
		coords = np.array([[0, 0]])
		
		for lat, lon in zip(latitude, longitude):
			lat = round(lat, 1)
			lon = round(lon, 1)
			latitude_range = np.linspace(lat - FoV_lat, lat + FoV_lat, 25)
			longitude_range = np.linspace(lon - FoV_lon, lon + FoV_lon, 25)
			
			latitude_range = np.atleast_2d(latitude_range).T
			longitude_range = np.atleast_2d(longitude_range).T
			
			arr = np.append(longitude_range, latitude_range, axis=1)
			coords = np.append(coords, arr, axis=0)
		
		plt.scatter(coords[:, 0], coords[:, 1], s=.5)
		plt.scatter(longitude, latitude, s=1)
	
	plt.show()


if __name__ == "__main__":
	Orbital_parameters = determine_orbit()
	SATS = satellite.create_sats(Orbital_parameters)
	
	for sat in SATS:
		satref = 25
		sat.reference_area = 6.59 + satref
		sat.drag_coefficient = 1.17
		sat.radiation_reference_area = 19.76 + satref
		sat.solar_pressure_coefficient = 1.2
		sat.rho = 1.5645e-13
		sat.dry_mass = 1254
	
	# tudat_env = tudat.set_up_and_run_tudat(SATS)
	# coverage(tudat_env, SATS)
	
	
