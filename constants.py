import numpy as np

J_2 = 1.08262668e-3     # [-] https://ai-solutions.com/_freeflyeruniversityguide/j2_perturbation.htm
R_E = 6378.137e3        # [m]
mu_E = 3.986004418e14   # [m^3/s^2]
T_E = 23 * 60 * 60 + 56 * 60 + 4.0916    # [s] universe today <- source
AU = 149_597_870.7e3    # [m] https://www.britannica.com/science/astronomical-unit
g_0 = 9.80665           # [m/s]
R = 8.31446261815324    # [J/mol K] https://www.britannica.com/science/decay-constant

## Conversions
m_to_deg = 1 / 111_139  # Conversion to find longitude/latitude degrees from surface distance https://sciencing.com/convert-latitude-longtitude-feet-2724.html

## Payload
FoV_PL = np.deg2rad(80)            # Field of view of the chosen payload [Rad]
Temporal_Res = 1                  # How often measurements are performed [day]

## Structure (Initial Estimates) - assuming a box shape and evenly distributed dry_mass
Struc_L = 4.2                       # Structure Length (m)
Struc_W = 2.0                       # Structure Width (m)
Struc_H = 2.0                       # Structure Height (m)
Struc_Ixx = 3750                    # Structure MMOI through lateral axis (kgm^2)
Struc_Iyy = 3750                    # Structure MMOI through other lateral axis (kgm^2)
Struc_Izz = 2535                    # Structure MMOI through longitudinal axis (kgm^2)
