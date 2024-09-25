import numpy as np
import matplotlib.pyplot as plt
from Simulation import method

def plot_velocity_comparison(method):
    """
    Plots the relativistic and non-relativistic speeds against the time index i.

    Args:
        method (str): Chosen numerical integration method.
    """
    t_rel = []
    speed_rel = []
    t_nonrel = []
    speed_nonrel = []

    # Construct file names based on the chosen numerical method
    relativistic_file = f"my_test_{method}_rel.npy"
    non_relativistic_file = f"my_test_{method}_nonrel.npy"

    # Load relativistic data
    relativistic_data = np.load(relativistic_file, allow_pickle=True)
    # Load non-relativistic data
    non_relativistic_data = np.load(non_relativistic_file, allow_pickle=True)

    niter_rel = len(relativistic_data)
    niter_nonrel = len(non_relativistic_data)

    # Extract data for relativistic case
    for i in range(niter_rel):
        time = relativistic_data[i][0]  # Extract time from simulation data
        t_rel.append(time)

        Rel_Particle1 = relativistic_data[i][1]  # Extract particle object from simulation data
        speed_rel.append(np.linalg.norm(Rel_Particle1.velocity))

    # Extract data for non-relativistic case
    for i in range(min(niter_rel, niter_nonrel)):
        time = non_relativistic_data[i][0]  # Extract time from simulation data
        t_nonrel.append(time)

        non_rel_Particle1 = non_relativistic_data[i][1]  # Extract particle object from simulation data
        speed_nonrel.append(np.linalg.norm(non_rel_Particle1.velocity))

   # Plot the data
    plt.plot(t_rel, speed_rel, label="Relativistic")
    plt.plot(t_nonrel, speed_nonrel, label="Non-relativistic")
    
    # Set font size for axis labels, title, and legend
    plt.xlabel("Time (s)", fontsize=30)  # Adjust fontsize as needed
    plt.ylabel("Speed (m/s)", fontsize=30)  # Adjust fontsize as needed
    plt.title(f"Comparison of Relativistic and Non-Relativistic Speeds (Method: {method})", fontsize=30)  # Adjust fontsize as needed
    plt.legend(fontsize=18)  # Adjust fontsize as needed

    # Set font size for tick labels on both axes
    plt.xticks(fontsize=20)  # Adjust fontsize as needed
    plt.yticks(fontsize=20)  # Adjust fontsize as needed

    plt.grid(True)
    plt.show()  

# Example usage
if __name__ == "__main__":
    chosen_method = method  
    plot_velocity_comparison(chosen_method)
