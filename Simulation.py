import numpy as np
import copy
from Particle import particle
from  Validation import Bz_check

use_relativistic = True    #Specify if you want to use relativistic functions or not 

"""
Cyclotron Simulation Script

This script performs a simulation of particle motion in a cyclotron. It initializes a particle, specifies simulation parameters such as time step and total simulation time, and executes the simulation loop. The particle's motion is updated based on the net force acting on it, and simulation data is saved to a NumPy file.

The script uses the particle class defined in the Particle module to represent the particle, and it relies on the numpy and copy modules for numerical computations and data manipulation.

Usage:
    Run this script to perform the cyclotron simulation and save the results to a NumPy file.


"""

# Create a proton particle object with specified parameters
Proton = particle(name="Proton", B= np.array([0,0,1.5]))


Proton.B = Bz_check(Proton.B)
theoretical_period = particle.calc_theoretical_period(Proton)


# Store the particle object in a list called Particles
Particles = [Proton]

# Initialize an empty list to store simulation data
Data = []
F_D_time = []
F_D_count=[]
D_time_data = []
D_count_data = []

# Set simulation parameters
deltaT = 6.958e-10 # Time step duration in seconds, calculated based on the time it takes a particle with escape velocity to travel the gap
tend =  1e-4                 # Total simulation time in seconds 
Ttotal = int(tend / deltaT)  # Total number of time steps

# Specify the numerical integration method
method = "Euler-Cromer"

# Choose filename based on the integration method
if use_relativistic == True:
    if method == "Euler":
        filename = "my_test_Euler_rel.npy"
    elif method == "Euler-Cromer":
        filename = "my_test_EulerCromer_rel.npy"
    elif method == "Verlet":
        filename = "my_test_Verlet_rel.npy"
    else:
        raise ValueError("No valid method has been selected")
else:
    if method == "Euler":
        filename = "my_test_Euler_nonrel.npy"
    elif method == "Euler-Cromer":
        filename = "my_test_Euler-Cromer_nonrel.npy"
    elif method == "Verlet":
        filename = "my_test_Verlet_nonrel.npy"
    else:
        raise ValueError("No valid method has been selected")

# Initialize time step index
i = 0

# Perform simulation loop until particle escapes or maximum time steps reached
while Particles[0].Escaped == False and i < Ttotal:
    # Print current time step, particle position, and velocity
    print(i, Particles[0].position, Particles[0].velocity)

    for p in Particles:

    #calculates the relativistic correction for mass and velocity
        p.calculate_relativistic_properties()
    
        # Calculate net force acting on each particle
        if use_relativistic == True:
            p.F_Net_Relativistic(Particles, deltaT, i)
        else:
            p.F_Net(Particles, deltaT, i)

        # Update position and velocity of each particle
        if use_relativistic == True:
            p.Update_Relativistic(deltaT, method, i, Particles)
        else:
            p.Update(deltaT, method, i, Particles)

        

    # # Append simulation data to Data list at regular intervals
   

   
   
    if i % 1 == 0:
        temp_list = [deltaT * (i + 1)]
        for p in Particles:
            temp_list.append(copy.deepcopy(p))
        Data.append(temp_list)

      


    
    
    # Increment time step index
    i = i + 1






################################################################
#################escaped the while loop#########################

# Check if simulation terminated due to particle escape or maximum time steps reached
final_mass = Particles[0].rel_mass
if i < Ttotal:
    print("The particle has escaped the cyclotron with a speed of:", np.linalg.norm(Particles[0].velocity))
    
    escp_vel = np.linalg.norm(Particles[0].velocity)
    final_speed = np.linalg.norm(Particles[0].velocity)
    final_position = np.linalg.norm(Particles[0].position)
    
else:
    print("The simulation has run over all time steps, the particle has NOT escaped")
    escp_vel = 0
    final_speed = np.linalg.norm(Particles[0].velocity)
    final_position = np.linalg.norm(Particles[0].position)

# Save simulation data to a NumPy file
np.save(filename, Data, allow_pickle=True)



print("File has been saved")
