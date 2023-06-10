import numpy as np
from matplotlib import pyplot as plt

import constants as const
import satellite


def determine_propellant_mass(
		dry_mass: float, mission_duration: float, sat_object
):
	"""
	:param dry_mass: Dry dry_mass of the satellite [kg]
	:param mission_duration: Lifetime of the satellite [year]
	:param sat_object: Satellite class object
	:return:
	"""
	
	def create_day_array(manoeuvre_type: str, DeltaV_H: float, I_sp: float, freq: float | None = None):
		if manoeuvre_type == "insertion":
			return np.array([[0, manoeuvre_type, DeltaV_H, I_sp]], dtype="object")
		elif manoeuvre_type == "deorbit":
			return np.array([[mission_days, manoeuvre_type, DeltaV_H, I_sp]], dtype="object")
		else:
			days = np.arange(int(1 / freq), mission_days + 1, int(1 / freq))
			days_array = np.zeros((len(days), 4), dtype="object")
			days_array[:, 0] = days
			days_array[:, 1] = manoeuvre_type
			days_array[:, 2] = DeltaV_H
			days_array[:, 3] = I_sp
			return days_array
	
	def Mp_tsiolkovsky(M_1: float, DeltaV: float, I_sp: float):
		return M_1 * (1 - np.exp(- DeltaV / (const.g_0 * I_sp)))
	
	def Mp_angular_momentum(DeltaH: float, I_sp):
		return DeltaH / (const.g_0 * I_sp)
	
	mission_days = int(365 * mission_duration)
	
	days_insertion = create_day_array("insertion", sat_object.DeltaV_insertion, sat_object.I_sp_insertion)
	days_desaturation = create_day_array("desaturation", sat_object.DeltaH_desaturation, sat_object.I_sp_desaturation, sat_object.freq_desaturation)
	days_maintenance = create_day_array("maintenance", sat_object.DeltaV_maintenance, sat_object.I_sp_maintenance, sat_object.freq_maintenance)
	days_deorbit = create_day_array("deorbit", sat_object.DeltaV_deorbit, sat_object.I_sp_deorbit)
	
	# Day   Type    Delta V/H     I_sp    M_0     M_1     M_p
	output = np.append(days_insertion, days_desaturation, axis=0)
	output = np.append(output, days_maintenance, axis=0)
	output = np.append(output, days_deorbit, axis=0)
	# Revert the sorting of the days, so we start with the end of the mission and calculate backwards
	output = output[output[:, 0].argsort()[::-1]]
	
	output_other_cols = np.zeros((len(output), 5), dtype="object")
	output = np.append(output, output_other_cols, axis=1)
	
	# Set last day of the mission the final value
	output[0, 5] = dry_mass
	
	for idx in range(len(output)):
		M_1 = output[idx, 5]
		I_sp = output[idx, 3]
		DeltaV_DeltaH = output[idx, 2]
		if output[idx, 1] != "desaturation":
			M_p = Mp_tsiolkovsky(M_1, DeltaV_DeltaH, I_sp)
		else:
			M_p = Mp_angular_momentum(DeltaV_DeltaH, I_sp)
		
		output[idx, 6] = M_p
		output[idx, 4] = M_1 + M_p
		if idx != len(output) - 1:
			output[idx + 1, 5] = M_1 + M_p
	
	M_p_tot = np.sum(output[:, 6])
	M_p_insertion = np.sum(output[np.where(output[:, 1] == "insertion"), 6])
	M_p_desaturation = np.sum(output[np.where(output[:, 1] == "desaturation")[0], 6])
	M_p_maintenance = np.sum(output[np.where(output[:, 1] == "maintenance"), 6])
	M_p_deorbit = np.sum(output[np.where(output[:, 1] == "deorbit"), 6])
	
	if __name__ == "__main__":
		
		print(
			f"###### Propulsion ######:\n"
			f"    Total propellant mass:        {M_p_tot} [kg]\n"
			f"    Insertion propellant mass:    {M_p_insertion} [kg]\n"
			f"    Desaturation propellant mass: {M_p_desaturation} [kg]\n"
			f"    Maintenance propellant mass:  {M_p_maintenance} [kg]\n"
			f"    De-orbit propellant mass:     {M_p_deorbit} [kg]\n")
	
		plt.scatter(output[:, 0], output[:, 5], s=1)
		plt.show()
	
	return M_p_tot


def min_thrust_main_engine(mass: float, DeltaV: float, t_pulse: float):
	acc = DeltaV / t_pulse
	Thrust = mass * acc
	
	print(f"Minimum thrust required for the main engine:\n"
	      f"T: {Thrust} [N]\n")
	


### Attitude control
def max_MIB_attitude_thruster(angular_momentum_per_desat: float, max_arm: float):
	N_engines = 2
	
	MIB = angular_momentum_per_desat / N_engines / max_arm
	
	print(f"Maximum minimum impulse bit for attitude thrusters: {MIB} [Ns]\n"
	      f"Assuming:   {N_engines} Burning at the time\n"
	      f"            {angular_momentum_per_desat} [Nms] desaturation\n"
	      f"            {max_arm} [m] as maximum arm\n")


if __name__ == "__main__":
	SAT = satellite.load_sat(1)
	
	dry_mass = 1254  # kg
	DeltaV_main = 179.108  # m/s
	t_pulse = 3600*24
	
	min_thrust_main_engine(dry_mass, DeltaV_main, t_pulse)
	
	angular_momentum_per_desat = 40  # [Nms]
	max_arm = 2  # [m]
	max_MIB_attitude_thruster(angular_momentum_per_desat, max_arm)
	
	mission_duration = 7.5  # [Year]
	
	# DeltaV_insertion = 0  # [m/s]
	SAT.DeltaH_desaturation = 3 * 31  # [Nms]
	SAT.DeltaV_maintenance = 3.75  # [m/s]
	# DeltaV_deorbit = 179.108  # [m/s]
	
	SAT.freq_desaturation = 1 / 2  # [/day]
	SAT.freq_maintenance = 1 / 84  # [/day]
	
	M_p_tot = determine_propellant_mass(dry_mass, mission_duration, SAT)
	
	