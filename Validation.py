
import numpy as np


def Bz_check(B):
    """
    Ensures that the magnetic field vector has only a z-component.

    Parameters:
        B (numpy.ndarray): The magnetic field vector.

    Returns:
        numpy.ndarray: The modified magnetic field vector with only the z-component.
    """
    # Check if the magnetic field vector already has only a z-component
    if B[0] == 0.0 and B[1] == 0.0:
        print("Mag field is valid")
        return B  # Return the input vector as it already has only a z-component
       

    # If not, create a new array with only the z-component of the magnetic field
    B_z = np.array([0.0, 0.0, B[2]] , dtype= float)
    print("Mag field has been amended to only have a z component")

    return B_z
