import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import constants as const
import tudat
import os
import pickle


class Satellite:
	def __init__(self):
		
		## Orbital parameters
		self.N_orbits: int = None
		""""Number of distinct types of orbits, changes in true anomaly and ascending node [-]"""
		self.N_sat: int = None
		"""Number of satellites in the constellation [-]"""
		self.j: int = None
		"""Number of orbits before the ground track is repeated [-]"""
		self.k: int = None
		"""Number of days before the ground track is repeated [Day]"""
		self.T: float = None
		"""Orbital period [s]"""
		self.r_p: float = None
		"""Peri-centre radius [m]"""
		self.e: float = None
		"""Eccentricity [-]"""
		self.a: float = None
		"""Semi-major axis [m]"""
		self.r_a: float = None
		"""Apo-centre radius [m]"""
		self.i: float = None
		"""Inclination angle [Rad]"""
		self.DeltaL: float = None
		"""Total shift in longitude per orbit [rad]"""
		self.omega: float = None
		"""Argument of the peri-centre [Rad] NEEDS REVISION"""
		self.OMEGA: float = None
		"""Longitude of the ascending node [Rad]"""
		self.nu: float = None
		"""True anomaly, angle from peri-centre of where the sat is now (starting point) [Rad]"""
		self.eclipse_time: float = None
		"""Maximum eclipse time during orbit [s]"""
		## Payload characteristics
		self.FoV_PL = None  # FoV of the payload [Rad]
		
		## Mission Requirements
		self.Temporal_Res = None  # Temporal resolution of the mission [days]
		
		self.mass = None  # Satellite mass [kg]
		self.frontal_area = None  # Frontal area [m^2]
		
		## Satellite Characteristics
		self.name = None
		self.reference_area = None
		self.drag_coefficient = None
		self.radiation_reference_area = None
		self.solar_pressure_coefficient = None
			
	def save_sat(self, sat_name: str, verbose=True) -> None:
		"""
		:param n_sat: Satellite number, int: 1 to 3
		:param verbose: Bool for printing confirmation/errors
		"""

		# # check satellite number/s validity
		# if n_sat in [1, 2, 3, 4, 5, 6]:
		# 	pass
		# else:
		# 	print(f"Invalid satellite number: {n_sat}, save failed")
		# 	return

		# set path to save folder
		cwd = os.getcwd().replace('\\', '/')			# unix compatibility
		dir_list = cwd.split('/')						# split directories
		try:
			dir_depth = dir_list[::-1].index("DSE_28")	# subdirectories away from "DSE_28"
		except ValueError:
			print('unable to find DSE_28 parent directory')
			return

		filepath = dir_depth * '../' + 'satellites/' + sat_name

		# save satellites
		with open(filepath, 'wb') as file:
			pickle.dump(self, file)

		print(f"Saved satellite: {sat_name}")
		return


def load_sat(n_sats: int | list[int] | str, verbose=True) -> list | None:
	"""
	:param n_sats: Satellite number, int: 1 to 6 OR str: "all"
	:param verbose: Bool for printing confirmation/errors
	:return:
	"""

	# check satellite number/s validity
	if type(n_sats) == int:
		n_sats = [n_sats]
	elif type(n_sats) == list:
		for n in n_sats:
			if n in [1, 2, 3, 4, 5, 6]:
				pass
			else:
				print(f"Invalid satellite number: {n}, save failed")
				return
	elif type(n_sats) == str:
		if n_sats.upper() == "ALL":
			n_sats = [1, 2, 3, 4, 5, 6]
	else:
		print("Unknown n_sat input, save failed")
		return

	# set path to save folder
	cwd = os.getcwd().replace('\\', '/')  # unix compatibility
	dir_list = cwd.split('/')  # split directories
	try:
		dir_depth = dir_list[::-1].index("DSE_28")  # subdirectories away from "DSE_28"
	except ValueError:
		print('unable to find DSE_28 parent directory')
		return

	filepath = dir_depth * '../' + 'satellites/sat'

	# load satellites
	sats = []
	for n_sat in n_sats:
		with open(filepath + str(n_sat), 'rb') as file:
			try:
				sats.append(pickle.load(file))
			except EOFError:
				print(f"Could not load sat {n_sat}")
				n_sats.remove(n_sat)
	print(f"Loaded satellites: {', '.join(map(str, n_sats))}")

	return sats


def create_sats(N_orbits: int, N_sat: int, j: int, k: int, DeltaL: float,
                a: float, e: float, i: float, r_p: float, r_a: float,
                T: float, OMEGA: np.ndarray, omega: float, nu: np.ndarray):
	SATS = np.empty(3, dtype="object")
	for n in range(N_sat):
		sat = Satellite()
		sat.N_orbits = N_orbits
		sat.N_sat = N_sat
		sat.j = j
		sat.k = k
		sat.DeltaL = DeltaL
		sat.a = a
		sat.e = e
		sat.i = i
		sat.r_p = r_p
		sat.r_a = r_a
		sat.T = T
		sat.OMEGA = OMEGA[n]
		sat.omega = omega
		sat.nu = nu[n]
		sat.name = f"Sat{n+1}"
		SATS[n] = sat
	
	return SATS


if __name__ == "__main__":
	pass
	## tudat
	# Sat1.tudat_env = tudat.tudat_environment()
	# Sat1.tudat_env.create_sat(Sat1, "Sat1")
	# Sat1.tudat_env.set_up_aerodynamics("Sat1")
	# Sat1.tudat_env.set_up_solar_pressure("Sat1")
	# Sat1.tudat_env.propagation_setup(["Sat1"])
	# Sat1.tudat_env.set_up_acceleration(["Sat1"])
	# Sat1.tudat_env.set_up_initial_states(["Sat1"])
	# Sat1.tudat_env.set_initial_state(Sat1, "Sat1", ["Sat1"])
	# Sat1.tudat_env.finalise_initial_states()
	# Sat1.tudat_env.define_vars_to_store("Sat1")
	# Sat1.tudat_env.define_propagator_settings()
	# Sat1.tudat_env.simulate()
	# Sat1.tudat_env.plot_ground_track()
