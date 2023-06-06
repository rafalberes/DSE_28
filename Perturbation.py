import satellite
import tudat
import numpy as np
import matplotlib.pyplot as plt

#Define Classed from Tudat and Satellite
tudat_env = tudat.tudat_environment(epochstart=2464328.500000)
sat = satellite.Satellite(i=82.9, e=0.0359,
	             FoV_hor=np.deg2rad(20), FoV_ver=np.deg2rad(10), mass=1254, frontal_area=2,
	             Res_spac=10, Res_temp=1, reference_area = 6.59, drag_coefficient=1.17, radiation_reference_area=19.76, solar_pressure_coefficient=1.2)

class Perturbations:
    def __init__(self, satlist):
        self.satlist = satlist

    def Setup(self, satlist):
        # Create satellite;
        tudat_env.create_sat(sat.mass, "Taking Control") #Dry mass in kg
        tudat_env.set_up_aerodynamics("Taking Control", 6.59,1.17) #Volume to cube: 16.9m^3-->6.59m^2 with Cd= 1.17 (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
        tudat_env.set_up_solar_pressure("Taking Control",19.76, 1.2) #Radiation reference area: 3*6.95m^2 with Qpr= [0,2] https://www.quora.com/What-are-the-physics-behind-radiation-pressure
        tudat_env.propagation_setup(["Taking Control"])
        tudat_env.set_up_acceleration(["Taking Control"])
        tudat_env.set_up_initial_states(["Taking Control"])
        tudat_env.set_initial_state()

if __name__ == "__main__":
	Sat1 = Perturbations([sat])