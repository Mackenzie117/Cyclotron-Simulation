import numpy as np
from matplotlib import pyplot as plt
from Particle import particle
from Simulation import filename, deltaT, theoretical_period
from Simulation import method
from Simulation import use_relativistic
from Checks import compare_escape_velocity
from Simulation import escp_vel
from Simulation import final_position
from Simulation import final_speed
from Simulation import final_mass

def load_simulation_data():
    """
    Loads simulation data from a NumPy file.

    Returns:
        tuple: A tuple containing time, x-position, y-position, D-count, and D-time data.
    """
    x1, y1 = [], []
    data = np.load(filename, allow_pickle=True)
    for time_step in data:
        Particle1 = time_step[1]
        x1.append(Particle1.position[0])
        y1.append(Particle1.position[1])
    return x1, y1

def calculate_period(x1, y1):
    """
    Calculates the oscillation period of the particle trajectory.

    Args:
        x1 (list): List of x-positions.
        y1 (list): List of y-positions.

    Returns:
        float: The average period.
        list: List of periods.
        list: List of indices.
    """
    average_period, periods, index = particle.find_oscillation_period(x1, y1, deltaT)
    print(periods)
    return average_period, periods, index

def plot_particle_trajectory(x1, y1, method_name):
    """
    Plots the trajectory of the particle in 2D space.

    Args:
        x1 (list): List of x-positions.
        y1 (list): List of y-positions.
        method_name (str): The name of the numerical integration method.
    """
    fig = plt.figure()
    plt.rcParams['font.size'] = '18'
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r'$X position$ (m)')
    ax.set_ylabel(r'$Y position$ (m)')
    plt.plot(x1, y1, label="Proton")
    ax.legend()
    max_abs_x = max(abs(min(x1)), abs(max(x1)))
    max_abs_y = max(abs(min(y1)), abs(max(y1)))
    ax.set_xlim(-max_abs_x, max_abs_x)
    ax.set_ylim(-max_abs_y, max_abs_y)
    plt.title(f"Particle Trajectory ({method_name})")
    plt.show()









def plot_period(periods, index, method_name):
    """
    Plots the period as a function of index.

    Args:
        periods (list): List of periods.
        index (list): List of indices.
        method_name (str): The name of the numerical integration method.
    """
    plt.plot(index[1:], periods, marker='o', linestyle='-', color='b')
    plt.xlabel('Timestep Index', fontsize=30)
    plt.ylabel('Time Difference (s)', fontsize=30)
    plt.title(f'Period as a Function of Index ({method_name})',fontsize=30)
    plt.grid(True)
    plt.show()
    plt.yticks(np.arange(min(periods), max(periods), step=0.1))
    plt.xticks(fontsize=20)  # Adjust fontsize as needed
    plt.yticks(fontsize=20)  # Adjust fontsize as needed












def visualise_simulation(method, use_relativistic):
    """
    Visualizes the results of the particle simulation.

    Args:
        method (str): The chosen numerical integration method.
        use_relativistic (bool): True if the simulation is relativistic, False otherwise.
    """
    x1, y1 = load_simulation_data()
    calc_period, periods, index = calculate_period(x1, y1)
    mod_diff = abs(theoretical_period - calc_period)
    if use_relativistic:
        method_name = f"{method} (Relativistic)"
    else:
        method_name = f"{method} (Non-Relativistic)"  # Corrected here
    print(f"Period analysis for {method_name}:")
    if calc_period > theoretical_period:
        print("The period from the simulation ({:.2e} s) is greater than the theoretical period ({:.2e} s).".format(calc_period, theoretical_period))
        print("The difference is ({:.2e} m/s).".format(mod_diff))
    elif calc_period < theoretical_period:
        print("The period from the simulation ({:.2e} s) is less than the theoretical period ({:.2e} s).".format(calc_period, theoretical_period))
        print("The difference is ({:.2e} s).".format(mod_diff))
    else:
        print("The period from the simulation ({:.2e} s) matches the theoretical period ({:.2e} s).".format(calc_period, theoretical_period))
    plot_particle_trajectory(x1, y1, method_name)
    plot_period(periods, index, method_name)

# Call the function to visualize the simulation results
visualise_simulation(method, use_relativistic)
compare_escape_velocity(escp_vel, final_speed, final_position, final_mass)


