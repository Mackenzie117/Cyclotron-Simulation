from Simulation import escp_vel
import numpy as np



def compare_escape_velocity(escp_vel, final_speed, final_position,final_mass):
    """
    Compares the escape velocity of the particle from the simulation
    to the theoretical escape velocity of a cyclotron.
    
    Parameters:
        escp_vel (float): The escape velocity of the particle from the simulation in meters per second.
        final_speed (float): The final speed of the particle in meters per second.
        final_position (float): The final position of the particle in the cyclotron in meters.
        final_mass (float): The final mass of the particle in kilograms.

    Returns:
        None
    """
    # Constants for theoretical calculation (you may need to adjust these based on your cyclotron setup)
    B_field = 1.5  # Magnetic field strength in Tesla
    r_Cycl = 0.5  # Cyclotron radius in meters
    q = 1.6e-19  # Charge of the particle in Coulombs
    mass = 1.67e-27  # Mass of the particle in kilograms


    if escp_vel != 0:
        # Calculate theoretical escape velocity using the correct formula
        theoretical_escape_velocity = (q * B_field * r_Cycl) / final_mass
        mod_diff = abs(escp_vel - theoretical_escape_velocity)

        # Compare escape velocities and print the results
        if escp_vel > theoretical_escape_velocity:
            print("The escape velocity from the simulation ({:.2e} m/s) is greater than the theoretical escape velocity ({:.2e} m/s).".format(escp_vel, theoretical_escape_velocity))
            print("The difference is ({:.2e}m/s). ".format(mod_diff))
        elif escp_vel < theoretical_escape_velocity:
            print("The escape velocity from the simulation ({:.2e} m/s) is less than the theoretical escape velocity ({:.2e} m/s).".format(escp_vel, theoretical_escape_velocity))
            print("The difference is ({:.2e}m/s). ".format(mod_diff))
        else:
            print("The escape velocity from the simulation ({:.2e} m/s) matches the theoretical escape velocity ({:.2e} m/s).".format(escp_vel, theoretical_escape_velocity))
    
    
    else:

        theoretical_speed = (q*B_field*final_position)/final_mass
        mod_diff = abs(final_speed - theoretical_speed)
        if final_speed > theoretical_speed:
            print("The speed from the simulation ({:.2e} m/s) is greater than the theoretical speed ({:.2e} m/s).".format(final_speed, theoretical_speed))
            print("The difference is ({:.2e}m/s). ".format(mod_diff))
        elif final_speed < theoretical_speed:
            print("The speed from the simulation ({:.2e} m/s) is less than the theoretical escape velocity ({:.2e} m/s).".format(final_speed, theoretical_speed))
            print("The difference is ({:.2e}m/s). ".format(mod_diff))
        else:
            print("The escape velocity from the simulation ({:.2e} m/s) matches the theoretical escape velocity ({:.2e} m/s).".format(final_speed, theoretical_speed)) 






