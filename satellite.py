import numpy as np
import os
import pickle
import csv


class Satellite:
	def __init__(self):
		
		## Orbital parameters
		self.N_orbits: int | None = None
		""""Number of distinct types of orbits, changes in true anomaly and ascending node [-]"""
		self.N_sat: int | None = None
		"""Number of satellites in the constellation [-]"""
		self.j: int | None = None
		"""Number of orbits before the ground track is repeated [-]"""
		self.k: int | None = None
		"""Number of days before the ground track is repeated [Day]"""
		self.T: float | None = None
		"""Orbital period [s]"""
		self.r_p: float | None = None
		"""Peri-centre radius [m]"""
		self.e: float | None = None
		"""Eccentricity [-]"""
		self.a: float | None = None
		"""Semi-major axis [m]"""
		self.r_a: float | None = None
		"""Apo-centre radius [m]"""
		self.i: float | None = None
		"""Inclination angle [Rad]"""
		self.DeltaL: float | None = None
		"""Total shift in longitude per orbit [rad]"""
		self.DeltaL1: float | None = None
		"""Shift in longitude per orbit due to rotation of the Earth [Rad]"""
		self.DeltaL2: float | None = None
		"""Shift in longitude per orbit due to J2 effect [Rad]"""
		self.omega: float | None = None
		"""Argument of the peri-centre [Rad] NEEDS REVISION"""
		self.OMEGA: float | None = None
		"""Longitude of the ascending node [Rad]"""
		self.nu: float | None = None
		"""True anomaly, angle from peri-centre of where the sat is now (starting point) [Rad]"""
		self.eclipse_time: float | None = None
		"""Maximum eclipse time during orbit [s]"""
		self.rho: float | None = None
		""""The mean density of the orbit [kg/m^3]"""
		
		### DeltaV requirements
		self.DeltaV_insertion: float = 0
		"""DeltaV required in the mission for orbit insertion [m/s]"""
		self.DeltaV_maintenance: float | None = None
		"""DeltaV required per orbit maintenance procedure [m/s]"""
		self.DeltaH_desaturation: float | None = None
		"""Delta angular momentum required per desaturation of reaction wheels [Nms]"""
		self.DeltaV_deorbit: float | None = None
		"""DeltaV required for the de-orbit manoeuvre [m/s]"""
		self.freq_desaturation: float | None = None
		"""Frequency of desaturation [/day]"""
		self.freq_maintenance: float | None = None
		"""Frequency of orbit maintenance [/day]"""
		
		### Propulsion
		self.I_sp_insertion: float = 285  # [s]
		"""Specific impulse of insertion engine [s]"""
		self.I_sp_desaturation: float = 285  # [s]
		"""Specific impulse of desaturation engine [s]"""
		self.I_sp_maintenance: float = 285  # [s]
		"""Specific impulse of maintenance engine [s]"""
		self.I_sp_deorbit: float = 285  # [s]
		"""Specific impulse of deorbit engine [s]"""
		
		## Payload characteristics
		self.FoV_PL = None  # FoV of the payload [Rad]
		
		## Mission Requirements
		self.Temporal_Res = None  # Temporal resolution of the mission [days]
		
		self.dry_mass = 1254  # Satellite dry_mass [kg]
		self.frontal_area = None  # Frontal area [m^2]
		self.Ld: float = None
		"""The total displacement allowed [m]"""
		## Satellite Characteristics
		self.name: str | None = None
		"""Name of the satellite"""
		self.reference_area: float | None = None
		"""Reference area of satellite's cross-section [m^2]"""
		self.drag_coefficient: float | None = None
		""""Drag coefficient of the satellite [-]"""
		self.radiation_reference_area: float | None = None
		""""Reference area of half the satellite's total area [m^2]"""
		self.solar_pressure_coefficient = None
		""""Solar pressure coefficient indicating to which scale radiation is absorbed or reflected [-]"""
		self.drag_coefficient: float | None = None
		"""Drag coefficient of the satellite [-]"""
		self.radiation_reference_area: float | None = None
		"""Reference area of half the satellite's total area [m^2]"""
		self.solar_pressure_coefficient: float | None = None
		"""Solar pressure coefficient indicating to which scale radiation is absorbed or reflected [-]"""
		self.lifetime: float | None = 7.5
		"""Lifetime of the satellite [y]"""
		self.reflectivity: float | None = None
		""""Reflectivity of the satellite"""

	def save_sat(self, name: str = '', verbose=False) -> None:
		"""
		:param name: Name of file to be saved, defaults to 'name' instance attribute of object
		:param verbose: Bool for printing confirmation/errors
		"""

		# set default name
		if name == '':
			name = self.name

		# set path to save folder
		cwd = os.getcwd().replace('\\', '/')
		dir_list = cwd.split('/')
		try:
			dir_depth = dir_list[::-1].index("DSE_28")
		except ValueError:
			if verbose:
				print('unable to find DSE_28 parent directory')
			return

		filepath = dir_depth * '../' + 'satellites'
		if os.path.isdir(filepath):
			pass
		else:
			os.mkdir(filepath)
		filepath += '/' + str(name)

		# save satellites
		with open(filepath + '.pkl', 'wb') as file:
			pickle.dump(self, file)

		with open(filepath + '.csv', 'w') as csv_file:
			writer = csv.writer(csv_file)
			for key, value in self.__dict__.items():
				writer.writerow([key, value])

		if verbose:
			print(f"Saved: {name}")
		return


def load_sat(sat_names: int | str | list[int | str], verbose: bool = False):
	"""
	:param sat_names: File names to be loaded (individual or list), integers automatically preceded with "Sat"
	:param verbose: Toggle for printing confirmation/errors
	:return: Satellite class object OR list of multiple
	"""

	# format name input
	if type(sat_names) == int:
		sat_names = ['Sat' + str(sat_names) + '.pkl']
	elif type(sat_names) == str:
		if sat_names[-4:] == '.pkl':
			sat_names = [sat_names]
		else:
			sat_names = [sat_names + '.pkl']
	elif type(sat_names) == list:
		for n, name in enumerate(sat_names):
			if type(name) == int:
				sat_names[n] = 'Sat' + str(name) + '.pkl'
	else:
		if verbose:
			print("Unknown sat_names input, load failed")
		return

	# set path to save folder
	cwd = os.getcwd().replace('\\', '/')  # unix compatibility
	dir_list = cwd.split('/')  # split directories
	try:
		dir_depth = dir_list[::-1].index("DSE_28")  # subdirectories away from "DSE_28"
	except ValueError:
		if verbose:
			print('unable to find DSE_28 parent directory')
		return

	filepath = dir_depth * '../' + 'satellites/'

	# load satellites
	sats = []
	for sat_name in sat_names:
		try:
			with open(filepath + str(sat_name), 'rb') as file:
				sats.append(pickle.load(file))
		except (FileNotFoundError, EOFError):
			if verbose:
				print(f"Could not load: {sat_name}")
			sats.append(None)
			sat_names.remove(sat_name)
	if verbose:
		print(f"Loaded: {', '.join(map(str, sat_names))}")
	if len(sats) == 1:
		sats = sats[0]

	return sats


def create_sats(Orbital_parameters: dict):
	SATS = np.empty(3, dtype="object")
	for n in range(Orbital_parameters["N_sat"]):
		sat = Satellite()
		sat.N_orbits = Orbital_parameters["N_orbits"]
		sat.N_sat = Orbital_parameters["N_sat"]
		sat.j = Orbital_parameters["j"]
		sat.k = Orbital_parameters["k"]
		sat.DeltaL = Orbital_parameters["DeltaL"]
		sat.DeltaL1 = Orbital_parameters["DeltaL1"]
		sat.DeltaL2 = Orbital_parameters["DeltaL2"]
		sat.a = Orbital_parameters["a"]
		sat.e = Orbital_parameters["e"]
		sat.i = Orbital_parameters["i"]
		sat.r_p = Orbital_parameters["r_p"]
		sat.r_a = Orbital_parameters["r_a"]
		sat.T = Orbital_parameters["T"]
		sat.OMEGA = Orbital_parameters["OMEGA"][n]
		sat.omega = Orbital_parameters["omega"]
		sat.nu = Orbital_parameters["nu"][n]
		sat.DeltaV_insertion = Orbital_parameters["DeltaV_insertion"]
		sat.DeltaV_deorbit = Orbital_parameters["DeltaV_deorbit"]
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
