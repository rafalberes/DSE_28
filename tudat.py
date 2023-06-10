"""
Based on: https://docs.tudat.space/en/latest/_src_getting_started/_src_examples/notebooks/propagation/perturbed_satellite_orbit.html?highlight=ground#Ground-track
"""

# Load standard modules
import numpy as np

import matplotlib
from matplotlib import pyplot as plt

# Load tudatpy modules
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel import constants
from tudatpy.util import result2array

import satellite


class tudat_environment:
    def __init__(self, epochstart: float = 0):
        self.initial_state = None
        self.time_hours = None
        self.dep_vars_array = None
        self.states_array = None
        self.propagator_settings = None
        self.acceleration_models = None
        self.dependent_variables_to_save = []
        self.initial_states_settings = None
        self.initial_states = None
        self.acceleration_settings = None
        self.central_bodies = None
        self.bodies_to_propagate = None
        
        # Load spice kernels
        spice.load_standard_kernels()
        
        # Set simulation start and end epochs
        self.simulation_start_epoch = epochstart   ## January 1 2035
        self.simulation_end_epoch = constants.JULIAN_DAY
        
        # Define string names for bodies to be created from default.
        self.bodies_to_create = ["Sun", "Earth", "Moon"] #"Mars", "Venus"
        
        # Use "Earth"/"J2000" as global frame origin and orientation.
        self.global_frame_origin = "Earth"
        self.global_frame_orientation = "J2000"
        
        # Create default body settings, usually from `spice`.
        self.body_settings = environment_setup.get_default_body_settings(
            self.bodies_to_create,
            self.global_frame_origin,
            self.global_frame_orientation)
        
        # Create system of selected celestial bodies
        self.bodies = environment_setup.create_system_of_bodies(self.body_settings)
    
    def create_sat(self, sat_object, name: str):
        # Create vehicle objects.
        self.bodies.create_empty_body(name)
        self.bodies.get(name).dry_mass = sat_object.dry_mass

    def set_up_aerodynamics(self,
                            sat: str,
                            reference_area: float, ## Initial vol
                            drag_coefficient: float  ## Initital drag coefficient
                            ):
        # Create aerodynamic coefficient interface settings, and add to vehicle
        aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
            reference_area, [drag_coefficient, 0, 0]
        )
        environment_setup.add_aerodynamic_coefficient_interface(
            self.bodies, sat, aero_coefficient_settings)
        
    def set_up_solar_pressure(self,
                              sat: str,
                              reference_area_radiation: float,
                              radiation_pressure_coefficient: float
                              ):
        # Create radiation pressure settings, and add to vehicle
        occulting_bodies = ["Earth"]
        radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
            "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies
        )
        environment_setup.add_radiation_pressure_interface(
            self.bodies, sat, radiation_pressure_settings)
        
    def propagation_setup(self, sats: list):
        self.bodies_to_propagate = sats
        self.central_bodies = ["Earth"]
        
    def set_up_acceleration(self, sats: list):
        # Define accelerations acting on satellite by Sun and Earth.
        accelerations_settings = dict(
            #Sun=[
                #propagation_setup.acceleration.cannonball_radiation_pressure(),
            #    propagation_setup.acceleration.point_mass_gravity()
            #],
            Earth=[
                propagation_setup.acceleration.spherical_harmonic_gravity(1, 5),
                #propagation_setup.acceleration.aerodynamic()
            ]#,
            #Moon=[
            #    propagation_setup.acceleration.point_mass_gravity()
            #]#, ## Should be removed due to insignificance compared to the others:
            # Mars=[
            #     propagation_setup.acceleration.point_mass_gravity()
            # ],
            # Venus=[
            #     propagation_setup.acceleration.point_mass_gravity()
            # ]
        )
        
        # Create global accelerations settings dictionary.
        self.acceleration_settings = {}
        for sat in sats:
            self.acceleration_settings[sat] = accelerations_settings

        # Create acceleration models.
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies,
            self.acceleration_settings,
            self.bodies_to_propagate,
            self.central_bodies)

    def set_up_initial_states(self, sats: list):
        self.initial_states = np.empty((len(sats), 6))
    
    def set_initial_state(self, sat_object, sat: str, sats: list):
        # Set initial conditions for the satellite that will be
        # propagated in this simulation. The initial conditions are given in
        # Keplerian elements and later on converted to Cartesian elements
        earth_gravitational_parameter = self.bodies.get("Earth").gravitational_parameter
        
        self.initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=earth_gravitational_parameter,
            semi_major_axis=sat_object.a,
            eccentricity=sat_object.e,
            inclination=sat_object.i,
            argument_of_periapsis=sat_object.omega,
            longitude_of_ascending_node=sat_object.OMEGA,
            true_anomaly=sat_object.nu,
        )

        self.initial_states[np.where(sats == sat)] = self.initial_state
        
    def finalise_initial_states(self):
        self.initial_states_settings = np.concatenate(self.initial_states)
        
    def define_vars_to_store(self, sat: str):
        # Define list of dependent variables to save
        self.dependent_variables_to_save = [
            propagation_setup.dependent_variable.total_acceleration(sat),
            propagation_setup.dependent_variable.keplerian_state(sat, "Earth"),
            propagation_setup.dependent_variable.latitude(sat, "Earth"),
            propagation_setup.dependent_variable.longitude(sat, "Earth"),
            #propagation_setup.dependent_variable.single_acceleration_norm(
            #    propagation_setup.acceleration.point_mass_gravity_type, sat, "Sun"
            #),
            #propagation_setup.dependent_variable.single_acceleration_norm(
            #    propagation_setup.acceleration.point_mass_gravity_type, sat, "Moon"
            #),
            # propagation_setup.dependent_variable.single_acceleration_norm(
            #     propagation_setup.acceleration.point_mass_gravity_type, sat, "Mars"
            # ),
            # propagation_setup.dependent_variable.single_acceleration_norm(
            #     propagation_setup.acceleration.point_mass_gravity_type, sat, "Venus"
            # ),
            propagation_setup.dependent_variable.single_acceleration_norm(
                propagation_setup.acceleration.spherical_harmonic_gravity_type, sat, "Earth"
            )#,
            #propagation_setup.dependent_variable.single_acceleration_norm(
            #    propagation_setup.acceleration.aerodynamic_type, sat, "Earth"
            #)#,
            #propagation_setup.dependent_variable.single_acceleration_norm(
            #    propagation_setup.acceleration.cannonball_radiation_pressure_type, sat, "Sun"
            #)
        ]
        
    def define_propagator_settings(self):
        # Create termination settings
        termination_condition = propagation_setup.propagator.time_termination(self.simulation_end_epoch)
        
        # Create numerical integrator settings
        fixed_step_size = 1.0
        integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)
        
        # Create propagation settings
        self.propagator_settings = propagation_setup.propagator.translational(
            self.central_bodies,
            self.acceleration_models,
            self.bodies_to_propagate,
            self.initial_state,
            self.simulation_start_epoch,
            integrator_settings,
            termination_condition,
            output_variables=self.dependent_variables_to_save
        )
    
    def simulate(self):
        # Create simulation object and propagate the dynamics
        dynamics_simulator = numerical_simulation.create_dynamics_simulator(
            self.bodies, self.propagator_settings
        )
        
        # Extract the resulting state and dependent variable history and convert it to an ndarray
        states = dynamics_simulator.state_history
        self.states_array = result2array(states)
        dep_vars = dynamics_simulator.dependent_variable_history
        self.dep_vars_array = result2array(dep_vars)
        
        self.time_hours = self.dep_vars_array[:, 0]/3600
        
    def plot_ground_track(self, sat_objects: list | np.ndarray, show_plot: bool = True, subset: int | None = None):
        for sat in sat_objects:
            if subset is None:
                latitude = np.rad2deg(self.dep_vars_array[:, 10]) + np.rad2deg(sat.OMEGA)
                longitude = np.rad2deg(self.dep_vars_array[:, 11])
            else:
                latitude = np.rad2deg(self.dep_vars_array[:, 10])[:subset] + np.rad2deg(sat.OMEGA)
                longitude = np.rad2deg(self.dep_vars_array[:, 11])[:subset]

            plt.scatter(longitude, latitude, s=1)
        
        if show_plot:
            plt.show()
        
    def create_plot(self, hours: float | int | None = None):
        if hours is None:
            subset = None
        else:
            subset = int(len(self.time_hours) / 24 * hours)
        plt.figure(figsize=(9, 5))
        plt.title(f"Ground track of satellite over {hours} hours")
        plt.xlabel('Longitude [deg]')
        plt.ylabel('Latitude [deg]')
        
        plt.yticks(np.arange(-91, 91, step=45))
        plt.xlim([-180, 180])
        plt.grid()
        plt.tight_layout()
        
        return subset


def set_up_and_run_tudat(sat_objects):
    sat1 = sat_objects[0]
    tudat_env = tudat_environment()
    tudat_env.create_sat(sat1, sat1.name)
    tudat_env.set_up_aerodynamics(sat1.name, sat1.reference_area, sat1.drag_coefficient)
    tudat_env.set_up_solar_pressure(sat1.name, sat1.radiation_reference_area, sat1.solar_pressure_coefficient)
    tudat_env.propagation_setup([sat1.name])
    tudat_env.set_up_acceleration([sat1.name])
    tudat_env.set_up_initial_states([sat1.name])
    tudat_env.set_initial_state(sat1, sat1.name, [sat1.name])
    tudat_env.finalise_initial_states()
    tudat_env.define_vars_to_store(sat1.name)
    tudat_env.define_propagator_settings()
    tudat_env.simulate()
    return tudat_env
    

if __name__ == "__main__":
    sat = satellite.load_sat(1)
    sat.reference_area = 2
    sat.drag_coefficient = 1
    sat.radiation_reference_area = 3
    sat.solar_pressure_coefficient = 0.5
    sat.mass = 1000
    print(sat.__dict__)
    tudat_env = set_up_and_run_tudat([sat])







