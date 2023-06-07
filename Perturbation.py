import satellite
import tudat
import numpy as np
import matplotlib.pyplot as plt
import constants as const

#Define Class from Tudat

class Perturbations:
    def __init__(self, satlist):
        self.satlist = satlist
        self.sat = satlist[0]
        #print(self.sat.__dict__)

    def Groundtrack(self):
        # Create satellite;
        tudat_env = tudat.tudat_environment(epochstart=0)   # 2464328.500000
        tudat_env.create_sat(self.sat, self.sat.name) #Dry mass in kg
        tudat_env.set_up_aerodynamics(self.sat.name, self.sat.reference_area, self.sat.drag_coefficient) #Volume to cube: 16.9m^3-->6.59m^2 with Cd= 1.17 (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
        tudat_env.set_up_solar_pressure(self.sat.name, self.sat.radiation_reference_area, self.sat.solar_pressure_coefficient) #Radiation reference area: 3*6.95m^2 with Qpr= [0,2] https://www.quora.com/What-are-the-physics-behind-radiation-pressure
        tudat_env.propagation_setup([self.sat.name])
        tudat_env.set_up_acceleration([self.sat.name]) #Need to update this to multiple satelite list where necessary!
        tudat_env.set_up_initial_states([self.sat.name]) #Same here!!!
        tudat_env.set_initial_state(self.sat, self.sat.name, [self.sat.name])
        tudat_env.finalise_initial_states()
        tudat_env.define_vars_to_store(self.sat.name)
        tudat_env.define_propagator_settings()
        tudat_env.simulate()
        tudat_env.plot_ground_track(show_plot=True)

    def ThirdBody(self,i,j): #Third body perturbations from Sun and Moon with error of O(e^2)
        OMEGAdotmoon = -0.00338 * np.cos(i)/j
        print("Moon effect on Ascending Note=",OMEGAdotmoon, "deg/day")
        OMEGAdotsun = -0.00154 * np.cos(i)/j
        print("Sun effect on Ascending Note=", OMEGAdotsun, "deg/day")

        omegadotmoon = 0.00169*(4-5*np.sin(i**2))/j
        print("Moon effect on Argument of Perigee=", omegadotmoon, "deg/day")
        omegadotsun = 0.00077 * (4 - 5 * np.sin(i ** 2)) / j
        print("Moon effect on Argument of Perigee=", omegadotsun, "deg/day")

    def NonSphericalEarth(self,a,i,e):
        OMEGAdotj2 = -2.06474e14*a**(-7/2)*np.cos(i)*(1-e**2)**-2
        print("J2 effect on Ascending Note=", OMEGAdotj2, "deg/day")
        omegadotj2 = 1.03237e14 * a ** (-7 / 2) * (4-5*np.sin(i)**2) * (1 - e ** 2) ** -2
        print("J2 effect on Ascending Note=", omegadotj2, "deg/day")

    def Aerodynamics(selfs):
        DeltaaRev = -2*np.pi*(self.sat.Cd)



if __name__ == "__main__":
    Sat1 = satellite.Satellite(name="Taking Control", i=86.27, e=0.04136,
                                  FoV_hor=np.deg2rad(40), FoV_ver=np.deg2rad(10), mass=1254, frontal_area=2,
                                  Res_spac=10, Res_temp=1, reference_area=6.59, drag_coefficient=1.17,
                                  radiation_reference_area=19.76, solar_pressure_coefficient=1.2, a = 7253137)
    print(Sat1.__dict__)
    satlist = [Sat1]
    Perb = Perturbations(satlist)
    #Perb.Groundtrack()
    Perb.ThirdBody(Sat1.i,Sat1.j)
    Perb.NonSphericalEarth()