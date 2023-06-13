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
		return M_1 * (-1 + np.exp(DeltaV / (const.g_0 * I_sp)))
	
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
	
	N_desaturation = len(output[np.where(output[:, 1] == "desaturation")[0]])
	N_main_engine = len(output) - N_desaturation
	
	if __name__ == "__main__":
		
		print(
			f"###### Propulsion ######:\n"
			f"    Total propellant mass:         {M_p_tot} [kg]\n"
			f"    Insertion propellant mass:     {M_p_insertion} [kg]\n"
			f"    Desaturation propellant mass:  {M_p_desaturation} [kg]\n"
			f"    Maintenance propellant mass:   {M_p_maintenance} [kg]\n"
			f"    De-orbit propellant mass:      {M_p_deorbit} [kg]\n"
			f"    Amount of ADCS thruster burns: {N_desaturation} [-]\n"
			f"    Amount of main thruster burns: {N_main_engine} [-]\n")
	
		# plt.scatter(output[:, 0], output[:, 5], s=1)
		# plt.show()
	
	return M_p_tot


def min_thrust_main_engine(mass: float, DeltaV: float, t_pulse: float):
	acc = DeltaV / t_pulse
	Thrust = mass * acc
	
	print(f"Minimum thrust required for the main engine:\n"
	      f"T: {Thrust} [N]\n")
	

### Attitude control
def max_MIB_attitude_thruster(angular_momentum_per_desat: float, max_arm: float):
	N_engines = 6
	
	MIB = angular_momentum_per_desat / N_engines / max_arm
	
	print(f"Maximum minimum impulse bit for attitude thrusters: {MIB} [Ns]\n"
	      f"Assuming:   {N_engines} Burning at the time\n"
	      f"            {angular_momentum_per_desat} [Nms] desaturation\n"
	      f"            {max_arm} [m] as maximum arm\n")
	
	
def design_prop_tank(M_p: float):
	def calc_fuel_to_ox_ratio():
		M_N = 14.067    # [g/mol]
		M_O = 15.9994   # [g/mol]
		M_C = 12.011    # [g/mol]
		M_H = 1.00797   # [g/mol]
		
		M_N2O = 2 * M_N + M_O
		M_C3H6 = 3 * M_C + 6 * M_H
		
		fuel_to_ox = M_C3H6 / (9 * M_N2O)
		return fuel_to_ox
	
	def my_ceil(a, precision=0):
		return np.true_divide(np.ceil(a * 10 ** precision), 10 ** precision)

	def calc_radius_tank(volume: float):
		return (volume / (4 / 3 * np.pi)) ** (1 / 3)
	
	def calc_t(P: float, r: float, sigma_max: float):
		t_min = 2e-4  # [m] https://www.spacematdb.com/spacemat/manudatasheets/TITANIUM%20ALLOY%20GUIDE.pdf
		safety_factor = 1.5
		return max(my_ceil(1 / 2 * (P * r) / (sigma_max * safety_factor), 4), t_min)
	
	def calc_mass_tank(r: float, t: float, rho: float):
		return 4 / 3 * ((r + t) ** 3 - r ** 3) * rho
	
	fuel_to_ox = calc_fuel_to_ox_ratio()
	
	m_N2O = M_p / (1 + fuel_to_ox)
	m_C3H6 = M_p - m_N2O
	
	rho_N2O = 1.22e6     # [g/m^3] https://pubchem.ncbi.nlm.nih.gov/compound/Nitrous-Oxide
	rho_C3H6 = 613.9e3   # [g/m^3] https://en.wikipedia.org/wiki/Propylene
	
	V_N2O = m_N2O / rho_N2O
	V_C3H6 = m_C3H6 / rho_C3H6
	
	r_N2O = calc_radius_tank(V_N2O)
	r_C3H6 = calc_radius_tank(V_C3H6)
	
	sigma_max = 880e6   # [Pa] https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mtp641
	rho_Ti6Al4V = 4430  # [kg/m^3] https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mtp641
	
	P_N2O = 7.5e6       # [Pa] https://arc.aiaa.org/doi/pdf/10.2514/1.47131
	P_C3H6 = 14.7e5     # [Pa] https://satcatalog.s3.amazonaws.com/components/1275/SatCatalog_-_Dawn_Aerospace_-_B20_Thruster_-_Datasheet.pdf?lastmod=20220429194642
	
	t_N2O = calc_t(P_N2O, r_N2O, sigma_max)
	t_C3H6 = calc_t(P_C3H6, r_C3H6, sigma_max)
	
	m_tank_N2O = calc_mass_tank(r_N2O, t_N2O, rho_Ti6Al4V)
	m_tank_C3H6 = calc_mass_tank(r_C3H6, t_C3H6, rho_Ti6Al4V)
	## N2 tank design
	# Assume that the maximum N2 tank pressure is the maximum operating pressure of the selected valve
	# Assume that, when all C3H6 is gone, N2 tank and C3H6 tanks have the same pressure
	P_N2 = 310e5  # [Pa] https://satcatalog.s3.amazonaws.com/components/1367/SatCatalog_-_Rafael_-_Latch_Valve_-_Datasheet.pdf?lastmod=20230109090051
	T = 273.15 + 20  # [K]
	
	V_N2 = (P_C3H6 * V_C3H6) / (P_N2 - P_C3H6)
	r_N2 = calc_radius_tank(V_N2)

	t_N2 = calc_t(P_N2, r_N2, sigma_max)
	m_tank_N2 = calc_mass_tank(r_N2, t_N2, rho_Ti6Al4V)
	
	n_N2 = P_N2 * V_N2 / (const.R * T)
	M_N2 = 2 * 14.067
	m_N2 = n_N2 * M_N2
	
	if __name__ == "__main__":
		print(
			f"Fuel to oxygen mass ration: {fuel_to_ox} [-]\n"
			f"Mass of N2O               : {m_N2O} [g]\n"
			f"Mass of C3H6              : {m_C3H6} [g]\n"
			f"Mass of N2                : {m_N2} [g]\n"
			f"Volume of N2O             : {V_N2O} [m^3]\n"
			f"                            {V_N2O * 1e3} [L]\n"
			f"Volume of C3H6            : {V_C3H6} [m^3]\n"
			f"                            {V_C3H6 * 1e3} [L]\n"
			f"Volume of N2 tank         : {V_N2} [m^3]\n"
			f"                          : {V_N2 * 1e3} [L]\n"
			f"Radius of N2O tank        : {r_N2O} [m]\n"
			f"Radius of C3H6 tank       : {r_C3H6} [m]\n"
			f"Radius of N2 tank         : {r_N2} [m]\n"
			f"Thickness of N2O tank     : {t_N2O} [m]\n"
			f"                            {t_N2O*1e3} [mm]\n"
			f"Thickness of C3H6 tank    : {t_C3H6} [m]\n"
			f"                          : {t_C3H6*1e3} [mm]\n"
			f"Thickness of N2 tank      : {t_N2} [m]\n"
			f"                          : {t_N2*1e3} [mm]\n"
			f"Mass of the N2O tank      : {m_tank_N2O} [kg]\n"
			f"Mass of the C3H6 tank     : {m_tank_C3H6} [kg]\n"
			f"Mass of the N2 tank       : {m_tank_N2} [kg]\n"
		)


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
	
	SAT.DeltaV_insertion = 0  # [m/s]
	# SAT.DeltaH_desaturation = 3 * 31.45  # [Nms]
	SAT.DeltaV_maintenance = 17.68  # [m/s]
	# SAT.DeltaV_deorbit = 179.108  # [m/s]
	
	SAT.freq_desaturation = 1 / 2  # [/day]
	SAT.freq_maintenance = 1 / 423  # [/day]
	
	M_p_tot = determine_propellant_mass(dry_mass, mission_duration, SAT)
	design_prop_tank(M_p_tot*10**3)
	
	