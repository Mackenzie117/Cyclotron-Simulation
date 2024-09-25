# Cyclotron Simulation Particle Class, Particle.py

The Cyclotron Simulation Particle Class is a Python class representing a charged particle in a cyclotron simulation.

## Usage

To use the Cyclotron Simulation Particle Class, simply include the `particle.py` file in your Python project directory.

```python
from Particle import particle

# Create a particle instance
particle_instance = particle()

# Perform simulations and calculations
particle_instance.F_Net(Particles, deltaT, i)
particle_instance.Update(deltaT, method)

For detailed documentation on the class attributes and methods, please refer to the docstrings within the particle.py file.

#############################################################################

# Cyclotron Simulation Script, Simulation.py

The Cyclotron Simulation Script is a Python script that performs a simulation of particle motion in a cyclotron. It initializes a particle, specifies simulation parameters such as time step and total simulation time, and executes the simulation loop. The particle's motion is updated based on the net force acting on it, and simulation data is saved to a NumPy file.

## Usage

To use the Cyclotron Simulation Script, simply run the script. Make sure you have the required dependencies installed (NumPy).

The script allows you to choose whether to use relativistic functions or not by setting the use_relativistic variable to True or False.

You can specify the numerical integration method (method) and the filename for saving the simulation data (filename). By default, the script saves the data to a NumPy file with a filename based on the integration method chosen.

########################################################################


# Particle Simulation Visualization Analysis.py

The Particle Simulation Visualization script is a Python script that visualizes the results of a particle simulation performed in a cyclotron. It loads simulation data from a NumPy file, calculates the oscillation period of the particle trajectory, and plots the trajectory and period as a function of index.

## Usage

To use the Particle Simulation Visualization script, simply run the script. Make sure you have the required dependencies installed (NumPy, Matplotlib).

###########################################################################

# Velocity Comparison, Rel_comp.,py

The Velocity Comparison Plot script is a Python script that compares the relativistic and non-relativistic speeds of a particle over time in a cyclotron simulation. It plots the speeds against the time index `i`.

## Usage

To use the Velocity Comparison Plot script, simply run the script. Make sure you have the required dependencies installed (NumPy, Matplotlib).

Ensure that the simualation has been run for both the rel and non_rel versions otherwise the plot wont be useful.


###########################################################################

# Checks.py, Escape Velocity Comparison

The Escape Velocity Comparison script is a Python script that compares the escape velocity of a particle from a cyclotron simulation to the theoretical escape velocity of a cyclotron. It calculates the theoretical escape velocity based on given parameters and compares it to the escape velocity obtained from the simulation.

## Usage

To use the Escape Velocity Comparison script, simply import it and call the `compare_escape_velocity` function with the required parameters. Make sure you have the required dependencies installed (NumPy).

##############################################################################

#Validation.py

The Bz_check function ensures that the magnetic field vector has only a z-component. It takes a numpy array representing the magnetic field vector as input and returns a modified array with only the z-component if necessary.

################################################################################

Unittest.py

The unit tests in this module verify the behavior of the particle class in the Particle module under various conditions. These tests ensure that the particle class handles invalid inputs and boundary cases appropriately, raising exceptions when necessary.

Test Cases
Negative Mass: Ensures that a ValueError is raised when attempting to initialize a particle with a negative mass.
Invalid Position: Verifies that a ValueError is raised when providing an invalid position argument.
Invalid Velocity: Checks if a ValueError is raised for an invalid velocity argument.
Invalid Force: Tests if a ValueError is raised for an invalid force argument.
Invalid Charge: Verifies that a ValueError is raised for an invalid charge argument.
Invalid Name: Ensures that a ValueError is raised when providing an invalid name argument.
Invalid Magnetic Field: Checks if a ValueError is raised for an invalid magnetic field argument.
Negative Cyclotron Radius: Verifies that a ValueError is raised for a negative cyclotron radius.
Negative Gap: Ensures that a ValueError is raised for a negative gap.
Gap Greater Than Radius: Checks if a ValueError is raised when the gap is greater than the cyclotron radius.

Contributors: Aiden Mackenzie Owen 