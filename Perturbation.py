import satellite
import tudat
import numpy as np
import matplotlib.pyplot as plt
import constants as const


# Define Class from Tudat

class Perturbations:
    def __init__(self, satlist):
        self.satlist = satlist
        self.sat = satlist[0]
        # # print(self.sat.__dict__)

    def Groundtrack(self, sat, name, reference_area, drag_coefficient, radiation_reference_area,
                    solar_pressure_coefficient):  # TBU
        # Create satellite;
        tudat_env = tudat.tudat_environment(epochstart=0)  # 2464328.500000
        tudat_env.create_sat(sat, name)  # Dry dry_mass in kg
        tudat_env.set_up_aerodynamics(name, reference_area,
                                      drag_coefficient)  # Volume to cube: 16.9m^3-->6.59m^2 with Cd= 1.17 (https://www.engineeringtoolbox.com/drag-coefficient-d_627.html)
        tudat_env.set_up_solar_pressure(name, radiation_reference_area,
                                        solar_pressure_coefficient)  # Radiation reference area: 3*6.95m^2 with Qpr= [0,2] https://www.quora.com/What-are-the-physics-behind-radiation-pressure
        tudat_env.propagation_setup([name])
        tudat_env.set_up_acceleration([name])  # Need to update this to multiple satelite list where necessary!
        tudat_env.set_up_initial_states([name])  # Same here!!!
        tudat_env.set_initial_state(sat, name, [name])
        tudat_env.finalise_initial_states()
        tudat_env.define_vars_to_store(name)
        tudat_env.define_propagator_settings()
        tudat_env.simulate()
        subset = tudat_env.create_plot(hours=None)
        tudat_env.plot_ground_track(sat_objects=[sat], show_plot=True, subset=subset)

    def ThirdBody(self, i: float, j: int):  # Third body perturbations from Sun and Moon with error of O(e^2)
        self.OMEGAdotmoon = -0.00338 * np.cos(i) / j
        # print("Moon effect on Ascending Note=", self.OMEGAdotmoon, "deg/day")
        self.OMEGAdotsun = -0.00154 * np.cos(i) / j
        # print("Sun effect on Ascending Note=", self.OMEGAdotsun, "deg/day")

        self.omegadotmoon = 0.00169 * (4 - 5 * np.sin(i ** 2)) / j
        # print("Moon effect on Argument of Perigee=", self.omegadotmoon, "deg/day")
        self.omegadotsun = 0.00077 * (4 - 5 * np.sin(i ** 2)) / j
        # print("Sun effect on Argument of Perigee=", self.omegadotsun, "deg/day")

    def NonSphericalEarth(self, a, i, e):
        self.OMEGAdotj2 = -2.06474e14 * a ** (-7 / 2) * np.cos(i) * (1 - e ** 2) ** -2
        # print("J2 effect on Ascending Note=", self.OMEGAdotj2, "deg/day")

        self.omegadotj2 = 1.03237e14 * a ** (-7 / 2) * (4 - 5 * np.sin(i) ** 2) * (1 - e ** 2) ** -2
        # print("J2 effect on Argument of Perigee=", self.omegadotj2, "deg/day")

    def Aerodynamics(self, Cd, reference_area, mass, a, rho, j):
        # Assumption circular orbit for simplification
        # Acceleration due to drag:
        V = np.sqrt(const.mu_E / a)  # Circular velocity of satellite
        self.dragacceleration = -(1 / 2) * rho * (Cd * reference_area / mass) * V ** 2
        # print("Acceleration due to drag=", self.dragacceleration, "m/s^2")

        DeltaaRev = -2 * np.pi * (Cd * reference_area / mass) * a ** 2 * rho
        self.DeltaaDay = DeltaaRev * j
        # print("Drag effect on Semi-Major Axis=", self.DeltaaDay, "m/day")

        DeltaVRev = np.pi * (Cd * reference_area / mass) * rho * a * V  # This should be considered for the orbit maintenance
        self.DeltaVDay = DeltaVRev * j
        # print("Drag effect on Velocity=", self.DeltaVDay, "DeltaV/day")

    def SolarRadiation(self, reference_area, mass, r):
        self.solaracceleration = -4.5e-6 * (1 + r) * reference_area / mass
        # print("Magnitude of acceleration due to solar radiation=", self.solaracceleration, "m/s^2")

    def OrbitMaintenance(self, a, e, lifetime, r_p, Ld):
        deltavdrag = self.DeltaVDay
        deltaOMEGA = self.OMEGAdotsun + self.OMEGAdotsun + self.OMEGAdotj2
        deltaomega = self.omegadotsun + self.omegadotmoon + self.omegadotj2
        omega = deltaomega * 365 * lifetime
        ## print(omega)
        DeltaVomega = 2 * np.sqrt(const.mu_E / (a * (1 - e ** 2))) * np.sin(np.deg2rad(deltaomega / 2))
        V = np.sqrt(const.mu_E / a)
        DeltaVomega2 = np.sqrt(2 * V ** 2 * (1 - np.cos(np.deg2rad(deltaomega))))  ##Cross-check the formula
        # print("DeltaV to maintain omega:", DeltaVomega)
        self.deltavdrag = deltavdrag * 365 * lifetime
        # print("DeltaV to maintain altitude due to drag:", self.deltavdrag)
        frac = self.solaracceleration / self.dragacceleration
        # print("Fraction solar/drag:", frac * 100, "[%]")
        self.totaldeltaV = (1 + frac) * self.deltavdrag
        # print("Total DeltaV to maintain altitude:", self.totaldeltaV)
        self.deltaadot = self.DeltaaDay * (1 + frac)
        # print("Semi-major axis reduction per day:", self.deltaadot)
        # self.daystomaintain = (r_p - const.R_E - 540e3)/-self.deltaadot #540 from the minimum height for the saturation of the torque storage
        # self.deltavpermaintian = self.totaldeltaV/(365*lifetime/self.daystomaintain)
        # print("Amount of days until maintenance:", self.daystomaintain, "with a deltaV of:", self.deltavpermaintian, "for ADCS")

        ##Max deviation
        Ld= 2*Ld #max two sided deviation
        self.deltaa  =np.sqrt(-4*a*self.deltaadot*Ld/(3*(const.OMEGA_e-deltaOMEGA)*const.R_E))
        # print("Every decrease of a=", self.deltaa, "orbit maintenance needs to be performed for payload req")
        self.Tsep = np.sqrt((-16*a*Ld)/(3*(const.OMEGA_e-deltaOMEGA)*const.R_E*self.deltaadot))
        self.Tsep = -2*self.deltaa/self.deltaadot
        self.deltaVsep = self.totaldeltaV/(365*lifetime/self.Tsep)
        # print("Every", self.Tsep, "days orbit maintenance is required with a deltaV of:", self.deltaVsep)
        
        return self.deltaVsep, self.Tsep

    def Maintenanceplots(self,a):
        hours = np.arange(0,10*24,1)
        deltaadot = self.deltaadot / 24
        semimajor = a + deltaadot*hours - const.R_E
        limit = a - self.deltaa - const.R_E
        for i in range(len(semimajor)):
            if semimajor[i] < limit:
                semimajor[i:] = semimajor[i:] + self.deltaa*2
        semimajor = semimajor*1e-3
        aline = [(a-const.R_E)*1e-3,(a-const.R_E)*1e-3]
        tline = [0, hours[-1]]
        plt.plot(hours,semimajor,c="blue", label= "Satellite")
        plt.plot(tline, aline, c="red", label="Target Semi-Major Axis")
        #plt.title("Semi-major axis over time [h]")
        plt.xlabel("Time [h]", fontsize=15)
        plt.ylabel("Semi-major axis altitude [km]", fontsize=15)
        plt.legend(fontsize=15)
        plt.show()




if __name__ == "__main__":
    # Sat1 = satellite.Satellite(name="Taking Control", i=86.27, e=0.04136,
    #                                        FoV_hor=np.deg2rad(40), FoV_ver=np.deg2rad(10), dry_mass=1254, frontal_area=2,
    #                               Res_spac=10, Res_temp=1, reference_area=6.59, drag_coefficient=1.17,
    #                               radiation_reference_area=19.76, solar_pressure_coefficient=1.2, a = 7253137, rho=3.7045e-13)
    Sats = satellite.load_sat(1)
    Sat1 = Sats
    satref = 3.5*2.35  # Initial maximum estimate of solar area
    satfrontal = 3.5*0.0155 #Initial width times height of the solar panels
    Sat1.reference_area = const.Struc_H * const.Struc_W + satfrontal
    # print(Sat1.reference_area)
    Sat1.drag_coefficient = 2  # Box with arrays average drag coefficient from SMAD, https://issfd.org/ISSFD_2017/paper/ISTS-2017-d-002__ISSFD-2017-002.pdf
    Sat1.radiation_reference_area = const.Struc_W*const.Struc_L + satref #+ const.Struc_H*const.Struc_L + const.Struc_H * const.Struc_W
    # print(Sat1.radiation_reference_area)
    Sat1.reflectivity = 0.58  # r = 0  for absorption; r = I  for specular and r := 0.4 for diffuse reflection. from SMAD, absorb=
    Sat1.solar_pressure_coefficient = Sat1.reflectivity + 1
    Sat1.lifetime = 7.5
    meanrho550 = 2.14e-13
    meanrho600 = 9.89e-14
    Sat1.rho = np.average([meanrho600, meanrho550])
    Sat1.rho = 1.188e-12  # Worst case density due to solar maximum
    Sat1.mass = 1363
    Sat1.Ld = 500
    # print(Sat1.__dict__)
    satlist = [Sat1]
    Perb = Perturbations(satlist)
    # Perb.Groundtrack(Sat1,Sat1.name, Sat1.reference_area, Sat1.drag_coefficient, Sat1.radiation_reference_area, Sat1.solar_pressure_coefficient)
    Perb.ThirdBody(Sat1.i, Sat1.j)
    Perb.NonSphericalEarth(Sat1.a, Sat1.i, Sat1.e)
    Perb.Aerodynamics(Sat1.drag_coefficient, Sat1.reference_area, Sat1.mass, Sat1.a, Sat1.rho, Sat1.j)
    Perb.SolarRadiation(Sat1.reference_area, Sat1.mass, Sat1.reflectivity)
    Perb.OrbitMaintenance(Sat1.a, Sat1.e, Sat1.lifetime, Sat1.r_p, Sat1.Ld)
    Perb.Maintenanceplots(Sat1.a)
    # print(Sat1.e**2*100)
