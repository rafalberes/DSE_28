import satellite
import tudat
import numpy as np
import matplotlib.pyplot as plt

#Define Classed from Tudat
tudat_env = tudat.tudat_environment(epochstart=2464328.500000)

class Perturbations:
    def __init__(self, satlist):
        self.satlist = satlist
        self.sat = satlist[0]
        #print(self.sat.__dict__)

    def Setup(self):
        # Create satellite;
        tudat_env.create_sat(self.sat, self.sat.name) #Dry mass in kg
        tudat_env.set_up_aerodynamics(self.sat.name, self.sat.reference_area, self.sat.drag_coefficient) #Volume to cube: 16.9m^3-->6.59m^2 with Cd= 1.17 (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
        tudat_env.set_up_solar_pressure(self.sat.name, self.sat.radiation_reference_area, self.sat.solar_pressure_coefficient) #Radiation reference area: 3*6.95m^2 with Qpr= [0,2] https://www.quora.com/What-are-the-physics-behind-radiation-pressure
        tudat_env.propagation_setup(self.satlist)
        tudat_env.set_up_acceleration(self.satlist)
        tudat_env.set_up_initial_states(self.satlist)
        tudat_env.set_initial_state()

if __name__ == "__main__":
    Sat1 = satellite.Satellite(name="Taking Control", i=82.9, e=0.0359,
                                  FoV_hor=np.deg2rad(20), FoV_ver=np.deg2rad(10), mass=1254, frontal_area=2,
                                  Res_spac=10, Res_temp=1, reference_area=6.59, drag_coefficient=1.17,
                                  radiation_reference_area=19.76, solar_pressure_coefficient=1.2)
    print(Sat1.__dict__)
    satlist = [Sat1]
    Perb = Perturbations(satlist)
    Perb.Setup()