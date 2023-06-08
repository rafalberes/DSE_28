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

    def Groundtrack(self, sat, name, reference_area, drag_coefficient, radiation_reference_area, solar_pressure_coefficient): #TBU
        # Create satellite;
        tudat_env = tudat.tudat_environment(epochstart=0)   # 2464328.500000
        tudat_env.create_sat(sat, name) #Dry mass in kg
        tudat_env.set_up_aerodynamics(name, reference_area, drag_coefficient) #Volume to cube: 16.9m^3-->6.59m^2 with Cd= 1.17 (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
        tudat_env.set_up_solar_pressure(name, radiation_reference_area, solar_pressure_coefficient) #Radiation reference area: 3*6.95m^2 with Qpr= [0,2] https://www.quora.com/What-are-the-physics-behind-radiation-pressure
        tudat_env.propagation_setup([name])
        tudat_env.set_up_acceleration([name]) #Need to update this to multiple satelite list where necessary!
        tudat_env.set_up_initial_states([name]) #Same here!!!
        tudat_env.set_initial_state(sat, name, [name])
        tudat_env.finalise_initial_states()
        tudat_env.define_vars_to_store(name)
        tudat_env.define_propagator_settings()
        tudat_env.simulate()
        tudat_env.plot_ground_track(show_plot=True)

    def ThirdBody(self,i,j): #Third body perturbations from Sun and Moon with error of O(e^2)
        self.OMEGAdotmoon = -0.00338 * np.cos(i)/j
        print("Moon effect on Ascending Note=",self.OMEGAdotmoon, "deg/day")
        self.OMEGAdotsun = -0.00154 * np.cos(i)/j
        print("Sun effect on Ascending Note=", self.OMEGAdotsun, "deg/day")

        self.omegadotmoon = 0.00169*(4-5*np.sin(i**2))/j
        print("Moon effect on Argument of Perigee=", self.omegadotmoon, "deg/day")
        self.omegadotsun = 0.00077 * (4 - 5 * np.sin(i ** 2)) / j
        print("Sun effect on Argument of Perigee=", self.omegadotsun, "deg/day")

    def NonSphericalEarth(self,a,i,e):
        self.OMEGAdotj2 = -2.06474e14*a**(-7/2)*np.cos(i)*(1-e**2)**-2
        print("J2 effect on Ascending Note=", self.OMEGAdotj2, "deg/day")

        self.omegadotj2 = 1.03237e14 * a ** (-7 / 2) * (4-5*np.sin(i)**2) * (1 - e ** 2) ** -2
        print("J2 effect on Argument of Perigee=", self.omegadotj2, "deg/day")

    def Aerodynamics(self, Cd, reference_area, mass, a, rho,j):
        # Assumption circular orbit for simplification
        # Acceleration due to drag:
        V = np.sqrt(const.mu_E / a)  # Circular velocity of satellite
        self.dragacceleration = -(1/2)*rho*(Cd*reference_area/mass)*V**2
        print("Acceleration due to drag=", self.dragacceleration, "m/s^2")

        DeltaaRev = -2*np.pi*(Cd*reference_area/mass)*a**2*rho
        self.DeltaaDay = DeltaaRev*j
        print("Drag effect on Semi-Major Axis=", self.DeltaaDay, "m/day")

        DeltaVRev = np.pi*(Cd*reference_area/mass)*rho*a*V #This should be considered for the orbit maintenance
        self.DeltaVDay = DeltaVRev*j
        print("Drag effect on Velocity=", self.DeltaVDay, "DeltaV/day")

    def SolarRadiation(self, reference_area, mass):
        self.solaracceleration = -4.5e-8*reference_area/mass
        print("Magnitude of acceleration due to solar radiation=", self.solaracceleration, "m/s^2")

    def OrbitMaintenance(self):
        deltavdrag = self.DeltaVDay
        deltaOMEGA = self.OMEGAdotsun + self.OMEGAdotsun + self.OMEGAdotj2



if __name__ == "__main__":
    # Sat1 = satellite.Satellite(name="Taking Control", i=86.27, e=0.04136,
    #                                        FoV_hor=np.deg2rad(40), FoV_ver=np.deg2rad(10), mass=1254, frontal_area=2,
    #                               Res_spac=10, Res_temp=1, reference_area=6.59, drag_coefficient=1.17,
    #                               radiation_reference_area=19.76, solar_pressure_coefficient=1.2, a = 7253137, rho=3.7045e-13)
    Sats = satellite.load_sat(1)
    Sat1 = Sats
    satref = 25
    Sat1.reference_area = 6.59 + satref
    Sat1.drag_coefficient = 1.17
    Sat1.radiation_reference_area = 19.76 + satref
    Sat1.solar_pressure_coefficient = 1.2
    Sat1.rho = 1.5645e-13
    Sat1.mass = 1254
    print(Sat1.__dict__)
    satlist = [Sat1]
    Perb = Perturbations(satlist)
    Perb.Groundtrack(Sat1,Sat1.name, Sat1.reference_area, Sat1.drag_coefficient, Sat1.radiation_reference_area, Sat1.solar_pressure_coefficient)
    Perb.ThirdBody(Sat1.i,Sat1.j)
    Perb.NonSphericalEarth(Sat1.a, Sat1.i, Sat1.e)
    Perb.Aerodynamics(Sat1.drag_coefficient, Sat1.reference_area, Sat1.mass, Sat1.a, Sat1.rho, Sat1.j)
    Perb.SolarRadiation(Sat1.reference_area, Sat1.mass)
