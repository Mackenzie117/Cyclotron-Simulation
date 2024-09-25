import numpy as np
import unittest
import scipy.constants as constants
import statistics
class particle:
    """
    A class representing a charged particle in a cyclotron simulation.

    Attributes:
        position (numpy.ndarray): The position vector of the particle in meters.
        velocity (numpy.ndarray): The velocity vector of the particle in meters per second.
        Force (numpy.ndarray): The force vector acting on the particle in Newtons.
        mass (float): The mass of the particle in kilograms. Default is 1.67e-27 kg (mass of a proton).
        q (float): The charge of the particle in coulombs. Default is 1.6e-19 C (charge of a proton).
        name (str): The name of the particle. Default is 'Particle'.
        B (numpy.ndarray): The magnetic field vector in Tesla.
        V (float): The voltage applied to the dees in the cyclotron in volts. Default is 50000 V.
        d (float): The separation distance between the dees in meters. Default is 0.05 m.
        r_Cycl (float): The radius of the cyclotron in meters. Default is 10 m.

   Methods:
        F_Net(Particles, deltaT, i):
            Calculates the net force acting on the particle at the current time step.
            Parameters:
                Particles (list): A list containing all particles in the simulation.
                deltaT (float): The time step duration in seconds.
                i (int): The current time step index.
        
        Update(deltaT, method):
            Updates the position and velocity of the particle based on the net force acting on it.
            Parameters:
                deltaT (float): The time step duration in seconds.
                method (str): The numerical integration method ('Euler' or other).

        F_Net_Relativistic(Particles, deltaT, i):
            Calculates the net force acting on the particle at the current time step using relativistic corrections.
            Parameters:
                Particles (list): A list containing all particles in the simulation.
                deltaT (float): The time step duration in seconds.
                i (int): The current time step index.

        Update_Relativistic(deltaT, method, i, Particles):
            Updates the position and velocity of the particle based on the net force acting on it using relativistic corrections.
            Parameters:
                deltaT (float): The time step duration in seconds.
                method (str): The numerical integration method ('Euler' or other).
                i (int): The current time step index.
                Particles (list): A list containing all particles in the simulation.

        find_oscillation_period(x, y, deltaT, tolerance=1e-2, d=0.05):
            Finds the period of oscillation from particle trajectory data.
            Parameters:
                x (list): List of x positions of the particle over time.
                y (list): List of y positions of the particle over time.
                deltaT (float): Time step duration in seconds.
                tolerance (float): Tolerance level for considering positions equal.
                d (float): Separation distance between the dees in meters. Default is 0.05 m.
    """


    def __init__(
    self,
    position = np.array([0,0,0], dtype= float),
    velocity=np.array([0, 0, 0], dtype=float),
    Force = np.array([0,0,0], dtype= float),
    mass = 1.67e-27,
    q =  1.6e-19,
    name = 'Particle',
    B = np.array([0.0,0.0,1.5]),
    V = 50000,
    d = 0.05,
    r_Cycl = 0.5,
    

    
    
    
    
    ):
        # Initialize attributes and run validation on the inputs:
        if position is None:
            self.position = position = np.array([0,0,0], dtype= float)

        else:
           
            if not isinstance(position,(list,tuple,np.ndarray)):
                raise ValueError("Position should be a list tuple or numpy array")
            # if len(position) !=3:
            #     raise ValueError("Position should have exactly 3 elements")

            self.position = np.array(position, dtype=float)

        if velocity is None:
            self.velocity = velocity=np.array([0, 0, 0], dtype=float)

        else:
            if not isinstance(velocity,(list,tuple,np.ndarray)):
                raise ValueError("Velocity should be a list tuple or numpy array")
            if len(velocity) !=3:
                raise ValueError("Velocity should have exactly 3 elements")
            self.velocity = np.array(velocity, dtype=float)

        if Force is None:
            self.Force = np.array([0,0,0], dtype= float)
        else:
            if not isinstance(Force,(list,tuple,np.ndarray)):
                raise ValueError("Force should be a list tuple or numpy array")
           
            self.Force = np.array(Force, dtype=float)
        
        if mass is None:
            raise ValueError("A mass must be chosen")
        
        
        if mass <0:
            raise ValueError("The mass must be positive")

        self.mass = mass

        if q is None:
            raise ValueError("A charge must be selected")
        
        self.q = q

        if not isinstance(name,str) or not name:
            raise ValueError("Name must be a non-empty string")
        
        self.name = name


        if not isinstance(B,(list,tuple,np.ndarray)):
            raise ValueError("B should be a list tuple or numpy array")
        if len(position) !=3:
            raise ValueError("B should have exactly 3 elements")
            
        self.B = np.array(B, dtype = float)
       

        self.V = V
        if r_Cycl < 0:
            raise ValueError("The radius of the cyclotron needs to be greater than zero")
        self.r_Cycl = r_Cycl

        if d<0: 
            raise ValueError("The gap needs to be greater than zero")
        if d > self.r_Cycl:
            raise ValueError("The gap cannot be wider the radius of the cyclotron")
        self.d = d

        self.Escaped = False

        self.rel_mass = 0

        self.D_count = 0
        self.D_time = 0
        

    def calculate_relativistic_properties(self):
        """
        Calculate relativistic properties: Lorentz factor, relativistic mass, and relativistic velocity.
        """
        # Calculate Lorentz factor
        v_squared = np.linalg.norm(self.velocity) ** 2
        gamma = 1 / np.sqrt(1 - v_squared / constants.c ** 2)

       
        # Combine relativistic velocity components into a numpy array
        

        # Update other relativistic properties
        self.gamma = gamma
        self.rel_mass = gamma * self.mass

        
    def F_Net(self, Particles, deltaT, i):

        """
        Calculates the net force acting on the particle at the current time step.

        Parameters:
            Particles (list): A list containing all particles in the simulation.
            deltaT (float): The time step duration in seconds.
            i (int): The current time step index.
        """

       
        
        self.Force = np.array([0,0,0], dtype=float)
        


        E_0 = np.array([self.V/self.d, 0 ,0], dtype= float)  #initialises the electric field as a vector and then calculates it this is done for the vector addition on line 87
       

        Angular_freq = self.q*np.linalg.norm(self.B)/self.mass # calculates the Cyclotron freq based on w =qB/m = 2pi/T

        if np.linalg.norm(self.position ) < self.r_Cycl: # this makes sure that the particle is still inside the cyclotron 
            self.Escaped = False
            if np.absolute(self.position[0]) < self.d/2: # This checks the position of the particle to see if it is between the two dees
                 
                self.Force = self.q*E_0*np.cos(Angular_freq*(deltaT*i)) + self.q*np.cross(self.velocity, self.B) # if particle is between the dees the force experienced is  from the electric field and magnetic

            else:
                self.Force = self.q*np.cross(self.velocity, self.B) # if the particle is in the dees the force is from the magnetic field 

        else:
            self.Escaped = True # This boolean indicates whether the particle has escaped 

        

        

            
        
    def Update(self, deltaT, method,i, Particles):
        """
        Updates the position and velocity of the particle based on the net force acting on it.

        Parameters:
            deltaT (float): The time step duration in seconds.
            method (str): The numerical integration method ('Euler' or other).
        """


        
        
        if method == "Euler":
            # Update position and velocity using Euler method
            self.position = self.position + self.velocity * deltaT
            self.velocity = self.velocity + self.Force/self.mass * deltaT  


        elif method == "Euler-Cromer":
            # Update velocity and position using Euler-Cromer method 
            self.velocity = self.velocity + self.Force/self.mass *deltaT 
            self.position = self.position + self.velocity * deltaT

        elif method == "Verlet":
            acceleration_current = self.Force / self.mass

            # Update position using Verlet method
            self.position += self.velocity * deltaT + 0.5 * acceleration_current * deltaT**2

            # Calculate acceleration at updated position
            self.Force = np.array([0,0,0], dtype=float)  # Reset force to zero before recalculating
            self.F_Net(Particles, deltaT,i)  # Recalculate force at new position
            acceleration_new = self.Force / self.mass

            # Update velocity using Verlet method
            self.velocity += 0.5 * (acceleration_current + acceleration_new) * deltaT
            

    def F_Net_Relativistic(self, Particles, deltaT, i):
        """
        Calculates the net force acting on the particle at the current time step using relativistic corrections.

        Parameters:
            Particles (list): A list containing all particles in the simulation.
            deltaT (float): The time step duration in seconds.
            i (int): The current time step index.
        """

        self.rel_Force = np.array([0, 0, 0], dtype=float)
        E_0 = np.array([self.V / self.d, 0, 0], dtype=float)

        # Calculate Lorentz factor
        v_squared = np.linalg.norm(self.velocity) ** 2
        gamma = 1 / np.sqrt(1 - v_squared / constants.c ** 2)

        Angular_freq = self.q * np.linalg.norm(self.B) / (gamma * self.rel_mass)

        if np.linalg.norm(self.position) < self.r_Cycl:
            self.Escaped = False
            if np.absolute(self.position[0]) < self.d / 2:
                self.rel_Force = self.q * E_0 * np.cos(Angular_freq * (deltaT * i)) + self.q * np.cross(self.velocity,self.B)

            else:
                self.rel_Force = self.q * np.cross(self.velocity, self.B)
        else:
            self.Escaped = True

    def Update_Relativistic(self, deltaT, method, i, Particles):
        """
        Updates the position and velocity of the particle based on the net force acting on it using relativistic corrections.

        Parameters:
            deltaT (float): The time step duration in seconds.
            method (str): The numerical integration method ('Euler' or other).
            i (int): The current time step index.
            Particles (list): A list containing all particles in the simulation.
        """

        self.calculate_relativistic_properties()

        if method == "Euler":
            self.position = self.position + self.velocity * deltaT
            self.velocity = self.velocity + self.rel_Force / self.rel_mass * deltaT

        elif method == "Euler-Cromer":
            self.velocity = self.velocity + self.rel_Force / self.rel_mass * deltaT
            self.position = self.position + self.velocity * deltaT

        elif method == "Verlet":
            acceleration_current = self.rel_Force / self.rel_mass

            self.position += self.velocity * deltaT + 0.5 * acceleration_current * deltaT ** 2

            self.rel_Force = np.array([0, 0, 0], dtype=float)
            self.F_Net_Relativistic(Particles, deltaT, i)
            acceleration_new = self.rel_Force / self.rel_mass

            self.velocity += 0.5 * (acceleration_current + acceleration_new) * deltaT















        
            


    
    def find_oscillation_period(x, y, deltaT, tolerance=1e-2, d=0.05):
        """
        Finds the period of oscillation from particle trajectory data.

        Parameters:
        x (list): List of x positions of the particle over time.
        y (list): List of y positions of the particle over time.
        deltaT (float): Time step duration in seconds.
        tolerance (float): Tolerance level for considering positions equal.

        Returns:
        float: The period of oscillation in seconds.
        list: List of time differences between consecutive occurrences of the particle returning to the same position.
        list: List of indices where the particle returns to the same position.
        """
        # Initialize lists to store the time indices where the particle reaches the same position twice
        time = []
        index = []

        # Iterate over the x and y positions to find the indices where the particle returns to the same position
        for i in range(1, len(y)):
            if x[i] > d/2 and abs(y[i]) < tolerance:  # if the particle is on the line y=0 for some tolerance and in the right hand dee
                index.append(i)
                time.append(i * deltaT)  # append the time that the particle is there

        # Calculate the time differences between consecutive occurrences of the particle returning to the same position
        periods = [time[i] - time[i - 1] for i in range(1, len(time))]

        if not periods:
            print("There are no instances of overlap, try decreasing your tolerance.")
            return None, None, None

        # Calculate the average time difference, which corresponds to the period of oscillation
        average_period = statistics.mean(periods)

        return average_period, periods, index

    
    def calc_theoretical_period(self):
        """" Calculate the theoretical period of oscillation based on cyclotron frequency.

        This method calculates the theoretical period of oscillation of a charged particle
        in a cyclotron based on the cyclotron frequency formula.

        Returns:
            float: The theoretical period of oscillation in seconds.
            
        Raises:
            ValueError: If the calculated theoretical period is negative."""

        theoretical_period = 2*np.pi*self.rel_mass/(self.q*np.linalg.norm(self.B))
        if theoretical_period < 0:
            raise ValueError("The time period cannot be negative")
        
        

        return theoretical_period


   



