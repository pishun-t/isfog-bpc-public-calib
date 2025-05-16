"""
Functions to calculate the in-situ states with depth
"""

import numpy as np
from elastic import fun_void_ratio, fun_mean_stress, fun_Gmax

def fun_sigeff_v(depths, gamma, option="uniform-dry"):
    """
    Function to calculate the effective vertical stress with depth
    """
    if option == "uniform-dry":
        # Uniform dry density
        sig_eff_v = gamma * depths
    else:
        raise NotImplementedError(f"Option {option} not implemented")
    return sig_eff_v


def fun_Gmax_with_depth(depths, void_ratio, Gref, K0, unit_weight):
    """
    Calculate the maximum elastic shear modulus with depth
     assuming constant void ratio and a constant K0
    Parameters
    ----------
    depths : np.1darray
        Depths in m
    void_ratio : float
        Void ratio of the material
    Gref : float
        Reference shear modulus
    K0 : float

    """
    p0_with_depth = (2*K0 + 1)/3 * fun_sigeff_v(depths, unit_weight)
    Gmax_with_depth = fun_Gmax(Gref, void_ratio, p0_with_depth)
    return Gmax_with_depth