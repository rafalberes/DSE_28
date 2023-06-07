import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import constants as const
import tudat
import os
import pickle

class Satellite:
	def __init__(self, name: str ="Taking Control", i: float = 0, e: float = 0,
	             FoV_hor: float = np.deg2rad(20) , FoV_ver: float =np.deg2rad(10), mass: float =1000, frontal_area: float =2,
	             Res_spac: float =10, Res_temp: float =1, reference_area: float =4, drag_coefficient: float =1.17,
				 radiation_reference_area: float =19.76, solar_pressure_coefficient: float =1.2):
		
		self.T = None  # Orbital period [s]
		self.r_p = None  # Peri-centre radius [m]
		self.e = e  # Eccentricity [-]
		self.a = None  # Semi-major axis [m]
		self.r_a = None  # Apo-centre radius [m]
		self.i = i  # Inclination angle [Rad]
		self.DeltaL1 = None  # Shift in longitude per obit due to rotation of the Earth
		self.DeltaL2 = None  # Shift in longitude per orbit due to the J2 effect [rad]
		self.DeltaL = None  # Total shift in longitude per orbit [rad]
		self.DeltaO = None  # Spacing in longitude needed of orbits within one earth-repeat orbit [rad]
		self.omega = np.deg2rad(70)  # Argument of the peri-centre [Rad] NEEDS REVISION
		self.OMEGA = np.deg2rad(20)  # Longitude of the ascending node [Rad]
		self.nu = np.deg2rad(45)  # True anomaly, angle from ascending node of where the sat is now [Rad]
		
		## Payload characteristics
		self.FoV_hor = FoV_hor  # Field of view in horizontal direction [Rad]
		self.FoV_ver = FoV_ver  # Field of view in vertical direction [Rad]
		self.PL_pix_hor = 2000  # Horizontal amount of pixels of the payload [-]
		self.PL_pix_ver = 1000  # Vertical amount of pixels of the payload [-]
		
		## Mission Requirements
		self.Res_spac = Res_spac  # Spacial resolution [m/pixel]
		self.Res_temp = Res_temp  # Temporal resolution [/hour]
		self.width_ver = None  # Latitude distance able to observe at a time [m]
		self.width_hor = None  # Longitude distance able to observe at a time [m]
		
		self.mass = mass  # Satellite mass [kg]
		self.frontal_area = frontal_area  # Frontal area [m^2]
		
		#self.calc_orbital_parameters()

		## Satellite Characteristics
		self.name = name
		self.reference_area = reference_area
		self.drag_coefficient = drag_coefficient
		self.radiation_reference_area = radiation_reference_area
		self.solar_pressure_coefficient =solar_pressure_coefficient
	
	def calc_orbital_parameters(self):
		self.calc_rp_and_DeltaO()
		
		self.a = self.r_p / (1 - self.e)
		self.r_a = self.a * (1 + self.e)
		
		self.T = 2 * np.pi * np.sqrt(self.a ** 3 / const.mu_E)
		
		self.DeltaL1 = -2 * np.pi * self.T / const.T_E
		self.DeltaL2 = (- 3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(self.i)) / (self.a ** 2 * (1 - self.e ** 2) ** 2)
		self.DeltaL = self.DeltaL1 + self.DeltaL2
		
	def calc_rp_and_DeltaO(self):
		## Horizontal
		self.width_hor = self.Res_spac * self.PL_pix_hor
		r_p_altitude_hor = 1 / 2 * self.width_hor / (np.tan(1 / 2 * self.FoV_hor))
		
		## Horizontal
		self.width_ver = self.Res_spac * self.PL_pix_ver
		r_p_altitude_ver = 1 / 2 * self.width_ver / (np.tan(1 / 2 * self.FoV_ver))
		
		r_p_altitude = min(r_p_altitude_hor, r_p_altitude_ver)
		print(f"Altitude r_p: {r_p_altitude * 10 ** -3} [km]")
		
		if r_p_altitude < 160e3:
			print("Equipment not sufficient")
			quit()
		else:
			self.r_p = r_p_altitude + const.R_E
			self.DeltaO = np.deg2rad(r_p_altitude * const.m_to_deg)
	
	def calc_DeltaL(self, AmountDeltaL,precission):
		## Determine range of DeltaLs:
		self.DeltaO = 1 * np.pi/180 ## Temporary 1 deg [rad]
		DeltaLRange = np.arange(1,AmountDeltaL,1)*self.DeltaO
		DeltaORounded = int(self.DeltaO*10**precission) ## Create integers by scaling OoM
		print(DeltaORounded)
		DeltaLRangeRounded = np.arange(1,AmountDeltaL,1)*DeltaORounded # DeltaL = integer * DeltaO
		#DeltaLRangeRounded = np.deg2rad(np.array([10,20,30,40,60,80]))
		## Find earth repeat orbit:
		PiRoundedInt = int(np.pi*2*10**precission) ## 62831 value
		#print(np.lcm([1,3,5,7],3)) ## As example of functions
		#DeltaLRangeRounded = int(30/180*np.pi*10**precission)
		lcm = np.lcm(DeltaLRangeRounded,PiRoundedInt, dtype='int64') * 10 ** (-precission)
		print(lcm)
		js = np.divide(lcm,DeltaLRangeRounded)
		ks = np.divide(lcm,PiRoundedInt)
		print("js:",js,"ks:",ks)
		#print(DeltaLRange)
		#print(js[0]*DeltaLRange[0], ks[0]*2*np.pi) ##Test consistency

	def calc_DeltaL2(self):
		self.DeltaO = 30*np.pi/180 #1/3 of the 90 degrees covered!
		self.e = 0








	
	def calc_e_from_DeltaL(self, DeltaL_required):
		def func_e_pos(e):
			func = -4 * np.pi ** 2 * np.sqrt((self.r_p ** 3 / ((1 - e) ** 3)) / const.mu_E) / const.T_E - \
			       (3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(self.i) /
			        (self.r_p ** 2 / ((1 - e) ** 2 * (1 - e ** 2) ** 2))) + DeltaL_required
			return func
		
		def func_e_neg(e):
			func = -4 * np.pi ** 2 * np.sqrt((self.r_p ** 3 / ((1 - e) ** 3)) / const.mu_E) / const.T_E - \
			       (3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(self.i) /
			        (self.r_p ** 2 / ((1 - e) ** 2 * (1 - e ** 2) ** 2))) - DeltaL_required
			return func
		
		self.DeltaL = DeltaL_required
		
		e = fsolve(func_e_pos, 0.5)[0]
		print(e)
		
		if e < 0:
			e = fsolve(func_e_neg, 0.5)[0]
			print(e)
		
		if e < 0 or e >= 1:
			print("Not realistic e:", e)
			quit()
		else:
			self.e = e
			self.calc_orbital_parameters()
			
	def save_sat(self, n_sat: int, verbose=True) -> None:
		"""
		:param n_sat: Satellite number, int: 1 to 6
		:param verbose: Bool for printing confirmation/errors
		"""

		# check satellite number/s validity
		if n_sat in [1, 2, 3, 4, 5, 6]:
			pass
		else:
			print(f"Invalid satellite number: {n_sat}, save failed")
			return

		# set path to save folder
		cwd = os.getcwd().replace('\\', '/')			# unix compatibility
		dir_list = cwd.split('/')						# split directories
		try:
			dir_depth = dir_list[::-1].index("DSE_28")	# subdirectories away from "DSE_28"
		except ValueError:
			print('unable to find DSE_28 parent directory')
			return

		filepath = dir_depth * '../' + 'satellites/sat' + str(n_sat)

		# save satellites
		with open(filepath, 'wb') as file:
			pickle.dump(self, file)

		print(f"Saved satellite: {n_sat}")
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

	return n_sats


if __name__ == "__main__":
	Sat1 = Satellite(np.deg2rad(80))
	#print(Sat1.__dict__)
	Sat1.calc_e_from_DeltaL(np.deg2rad(-30))
	#print(Sat1.__dict__)
	Sat1.calc_DeltaL(AmountDeltaL=200,precission=5)
	
	## tudat
	Sat1.tudat_env = tudat.tudat_environment()
	Sat1.tudat_env.create_sat(Sat1, "Sat1")
	Sat1.tudat_env.set_up_aerodynamics("Sat1")
	Sat1.tudat_env.set_up_solar_pressure("Sat1")
	Sat1.tudat_env.propagation_setup(["Sat1"])
	Sat1.tudat_env.set_up_acceleration(["Sat1"])
	Sat1.tudat_env.set_up_initial_states(["Sat1"])
	Sat1.tudat_env.set_initial_state(Sat1, "Sat1", ["Sat1"])
	Sat1.tudat_env.finalise_initial_states()
	Sat1.tudat_env.define_vars_to_store("Sat1")
	Sat1.tudat_env.define_propagator_settings()
	Sat1.tudat_env.simulate()
	Sat1.tudat_env.plot_ground_track()
	