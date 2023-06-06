import numpy as np

J_2 = 1.08262668e-3     # [-] https://ai-solutions.com/_freeflyeruniversityguide/j2_perturbation.htm
R_E = 6378.137e3        # [m]
mu_E = 3.986004418e14   # [m^3/s^2]
T_E = 23 * 60 * 60 + 56 * 60 + 4.0916    # [s] universe today <- source
AU = 149_597_870.7e3    # [m] https://www.britannica.com/science/astronomical-unit

## Conversions
m_to_deg = 1 / 111_139  # Conversion to find longitude/latitude degrees from surface distance https://sciencing.com/convert-latitude-longtitude-feet-2724.html

## Orbital parameters
j = 14                          # Number of orbits before ground repeat [-]
k = 1                           # Number of days before ground repeat [-]
DeltaL = np.deg2rad(25.714)     # Ground track shift per orbit [Rad]
r_p = 610e3                     # Peri-centre radius [m]
a = 870e3 + R_E                 # Semi-major axis [m]
i = np.deg2rad(82.9)            # Inclination [Rad]
N = 4                           # Number of satellites in the constellation [-]
e = 0.916                       # Eccentricity [-]
r_a = 3783e3                    # Apo-centre radius [m]
