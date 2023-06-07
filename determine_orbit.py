from math import ceil

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import constants as const


def determine_orbit():
	# We have selected a three type orbit. From this determine the maximum available altitude to achieve this.
	# Knowing the j is 14, from trial and error
	
	N_orbits = 3
	
	j = 14
	k = 1
	N_sat = ceil(k * N_orbits / const.Temporal_Res)
	DeltaL = 2*np.pi * k / j

	lat_measured = np.rad2deg(DeltaL / N_orbits) / const.m_to_deg
	
	accuracy = 25_000
	r_p_alt = (.5 * lat_measured / np.tan(.5 * const.FoV_PL)) // accuracy * accuracy + accuracy
	
	r_p = r_p_alt + const.R_E
	
	arange = np.arange(r_p, r_p + 2000e3, 10000)
	
	results = np.array([[0, 0, 0, 0]])
	
	for a in arange:
		Part1 = ((-1) * k / j * 2 * np.pi + 2 * np.pi * 2 * np.pi * np.sqrt(a ** 3 / const.mu_E) / const.T_E)
		Part2 = (-a ** 2 * (1 - (1 - r_p / a) ** 2) ** 2 / (3 * np.pi * const.J_2 * const.R_E ** 2))
		i = np.arccos(Part1 * Part2)
		results = np.append(results, [[k, j, a, i]], axis=0)
	
	k, j, a, i = results[np.where(results[:, 3] > np.deg2rad(80))[0][0]]
	
	OMEGA = np.arange(0, DeltaL, DeltaL / N_orbits)
	
	## Argument of peri-centre
	lat = np.deg2rad(70)
	# Sketchy source: https://physics.stackexchange.com/questions/47094/convert-latitude-of-lowest-altitude-to-argument-of-perigee
	omega = np.arcsin((np.sin(lat))/(np.sin(i)))
	
	e = 1 - r_p / a
	r_a = a / (1 + e)
	
	T = 2 * np.pi * np.sqrt(a ** 3 / const.mu_E)
	
	nu = np.array([0, 0, 0])
	
	DeltaL1 = -4*np.pi ** 2 * np.sqrt(a ** 3 / const.mu_E) / const.T_E
	DeltaL2 = - (3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(i)) / (a ** 2 * (1 - e ** 2) ** 2)

	return N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu


def store_orbital_parameters(N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu):
	Sats[0].N_orbits = N_orbits


def print_orbital_parameters():
	N_orbits, N_sat, j, k, DeltaL, DeltaL1, DeltaL2, a, e, i, r_p, r_a, T, OMEGA, omega, nu = determine_orbit()
	print(f"{N_orbits} types of orbits are used that will repeat in {k} days with {N_sat} satellites.\n"
	      f"In the {k} days, {j} orbits are completed.\n"
	      f"These orbits have:\n"
	      f"DeltaL                        : {np.rad2deg(DeltaL)} [Deg]\n"
	      f"semi-major axis, a            : {a} [m]\n"
	      f"Eccentricity, e               : {e} [-]\n"
	      f"Inclination, i                : {np.rad2deg(i)} [Deg]\n"
	      f"Peri-centre radius, r_p       : {r_p} [m]\n"
	      f"                                {r_p-const.R_E} [m]\n"
	      f"Apo-centre radius, r_a        : {r_a} [m]\n"
	      f"                                {r_a-const.R_E} [m]\n"
	      f"Orbital period, T             : {T} [s]\n"
	      f"                                {T/60} [min]\n"
	      f"Ascending nodes, OMEGA        : {np.rad2deg(OMEGA)} [Deg]\n"
	      f"Argument of peri-centre, omega: {np.rad2deg(omega)} [Deg]\n"
	      f"True anomaly, nu              : {np.rad2deg(nu)} [Deg]"
	      )
		

if __name__ == "__main__":
	print_orbital_parameters()
	