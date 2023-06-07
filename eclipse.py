import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

import constants as const


def calc_eclipse_times(sats):
	eclipse_time = calc_eclipse_time(sats[0])
	for sat in sats:
		sat.eclipse_time = eclipse_time


def calc_eclipse_time(sat, show: bool = False):
	OMEGA = sat.OMEGA
	omega = sat.omega
	i = sat.i
	a = sat.a
	e = sat.e
	
	n_array = np.arange(0, 366, 1)
	theta = np.linspace(0, 2 * np.pi, 360)
	
	OMEGA_array = np.linspace(0, 2 * np.pi, int(sat.j))
	max_duration_array = np.empty(np.shape(OMEGA_array))
	
	for idx, OMEGA in enumerate(OMEGA_array):
		l1 = np.cos(OMEGA) * np.cos(omega) - np.sin(OMEGA) * np.sin(omega) * np.cos(i)
		m1 = np.sin(OMEGA) * np.cos(omega) + np.cos(OMEGA) * np.sin(omega) * np.cos(i)
		n1 = np.sin(omega) * np.sin(i)
		
		l2 = -np.cos(OMEGA) * np.sin(omega) - np.sin(OMEGA) * np.cos(omega) * np.cos(i)
		m2 = -np.sin(OMEGA) * np.sin(omega) + np.cos(OMEGA) * np.cos(omega) * np.cos(i)
		n2 = np.cos(omega) * np.sin(i)
		
		X_sun, Y_sun, Z_sun = calc_coord_sun(n_array)
	
		## Get the time, matching with the theta's
		time, E = calc_time(a, e, theta)
		dt = np.diff(time)
		avg_dt = np.average(dt)
		dt = np.append(dt, avg_dt)
		duration_array = np.zeros(np.shape(n_array))
	
		for n in n_array:
			eclipse_array = (n+1) * np.ones(np.shape(theta))
			a_bar = calc_a_bar(l1, m1, n1, X_sun[n], Y_sun[n], Z_sun[n])
			b_bar = calc_b_bar(l2, m2, n2, X_sun[n], Y_sun[n], Z_sun[n])
			
			psi = calc_psi(theta, a_bar, b_bar)
			S = calc_shadow_function(a, e, theta, a_bar, b_bar)
			
			eclipse_array[(S > 0)*(psi > .5*np.pi)] = 0
			duration_array[n] = np.sum(dt[eclipse_array == 0])
			
		max_duration_array[idx] = np.max(duration_array)
	
	return np.max(max_duration_array)
	

def calc_time(a, e, theta_array):
	E_array = np.empty(np.shape(theta_array))

	for idx, theta in enumerate(theta_array):
		E_fun = lambda E: np.sqrt((1+e)/(1-e)) * np.tan(E / 2) - np.tan(theta / 2)
		E = fsolve(E_fun, 3.14)
		if E < 0:
			E += 2*np.pi
		E_array[idx] = E
		
	t = 1 / np.sqrt(const.mu_E / (a**3)) * (E_array - e * np.sin(E_array))
	
	return t, E_array


def calc_psi(theta: np.ndarray, a_bar: float, b_bar: float):
	cos_psi = a_bar * np.cos(theta) + b_bar * np.sin(theta)
	psi = np.arccos(cos_psi)
	return psi


def calc_a_bar(l1: float, m1: float, n1: float, X_sun: float, Y_sun: float, Z_sun: float):
	return (l1 * X_sun + m1 * Y_sun + n1 * Z_sun) / const.AU


def calc_b_bar(l2: float, m2: float, n2: float, X_sun: float, Y_sun: float, Z_sun: float):
	return (l2 * X_sun + m2 * Y_sun + n2 * Z_sun) / const.AU


def calc_coord_sun(n: np.ndarray):
	g = np.deg2rad(357.528) + np.deg2rad(0.9856003) * n
	epsilon = np.deg2rad(23.439) - np.deg2rad(4e-7) * n
	L = np.deg2rad(280.46) + np.deg2rad(0.9856474) * n
	Lambda = L + np.deg2rad(1.915) * np.sin(g) + np.deg2rad(0.02) * np.sin(2*g)
	R = const.AU * (1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2*g))
	
	X_sun = R * np.cos(Lambda)
	Y_sun = R * np.cos(epsilon) * np.sin(Lambda)
	Z_sun = R * np.sin(epsilon) * np.sin(Lambda)
	
	return X_sun, Y_sun, Z_sun


def calc_shadow_function(a: float, e: float, theta: np.ndarray, a_bar: float, b_bar: float):
	p = a * (1 - e ** 2)
	return const.R_E ** 2 * (1 + e * np.cos(theta)) ** 2 + p ** 2 * (a_bar * np.cos(theta) + b_bar * np.sin(theta)) ** 2 - p ** 2


if __name__ == "__main__":
	import satellite
	SATS = satellite.load_sat([1, 2, 3])
	print(SATS[0].__dict__)
	duration_eclipse_max = calc_eclipse_time(SATS[0])
	print("Max duration eclipse: ", duration_eclipse_max, "[s]")
	print("Max duration eclipse: ", duration_eclipse_max/60, "[min]")
	print("Max duration eclipse: ", duration_eclipse_max / (2*np.pi*np.sqrt(SATS[0].a**3/const.mu_E)), "[%]")
	