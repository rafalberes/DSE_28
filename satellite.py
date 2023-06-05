import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import constants as const


class Satellite:
	def __init__(self, i=0, e=0,
	             FoV_hor=np.deg2rad(20), FoV_ver=np.deg2rad(10), mass=1000, frontal_area=2,
	             Res_spac=10, Res_temp=1):
		self.T = None  # Orbital period [s]
		self.r_p = None  # Peri-centre radius [m]
		self.e = e  # Eccentricity [-]
		self.a = None  # Semi-major axis [m]
		self.r_a = None  # Apo-centre radius [m]
		self.i = i  # Inclination angle [Rad]
		self.DeltaL1 = None
		self.DeltaL2 = None
		self.DeltaL = None
		
		## Payload characteristics
		self.FoV_hor = FoV_hor  # Field of view in horizontal direction [Rad]
		self.FoV_ver = FoV_ver  # Field of view in vertical direction [Rad]
		self.PL_pix_hor = 2000  # Horizontal amount of pixels of the payload [-]
		self.PL_pix_ver = 1000  # Vertical amount of pixels of the payload [-]
		
		## Mission Requirements
		self.Res_spac = Res_spac  # Spacial resolution [m/pixel]
		self.Res_temp = Res_temp  # Temporal resolution [/hour]
		
		self.mass = mass  # Satellite mass [kg]
		self.frontal_area = frontal_area  # Frontal area [m^2]
		
		self.calc_orbital_parameters()
	
	def calc_orbital_parameters(self):
		self.calc_rp_from_res_spac()
		
		self.a = self.r_p / (1 - self.e)
		self.r_a = self.a * (1 + self.e)
		
		self.T = 2 * np.pi * np.sqrt(self.a ** 3 / const.mu_E)
		
		self.DeltaL1 = -2 * np.pi * self.T / const.T_E
		self.DeltaL2 = (- 3 * np.pi * const.J_2 * const.R_E ** 2 * np.cos(self.i)) / (self.a ** 2 * (1 - self.e ** 2) ** 2)
		self.DeltaL = self.DeltaL1 + self.DeltaL2
		
	def calc_rp_from_res_spac(self):
		## Horizontal
		width_hor = self.Res_spac * self.PL_pix_hor
		r_p_altitude_hor = 1 / 2 * width_hor / (np.tan(1 / 2 * self.FoV_hor))
		
		## Horizontal
		width_ver = self.Res_spac * self.PL_pix_ver
		r_p_altitude_ver = 1 / 2 * width_ver / (np.tan(1 / 2 * self.FoV_ver))
		
		self.r_p = min(r_p_altitude_hor, r_p_altitude_ver) + const.R_E
	
	def calc_DeltaO(self):
		pass
	
	def calc_DeltaL(self):
		pass
	
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
		
		e = fsolve(func_e_pos, 0)[0]
		print(e)
		
		if e < 0:
			e = fsolve(func_e_neg, 0)[0]
			print(e)
		
		if e < 0 or e >= 1:
			print("Not realistic e:", e)
		else:
			self.e = e
			self.calc_orbital_parameters()
			

if __name__ == "__main__":
	Sat1 = Satellite(np.deg2rad(70))
	print(Sat1.__dict__)
	Sat1.calc_e_from_DeltaL(np.deg2rad(-23))
	print(Sat1.__dict__)
	