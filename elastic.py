"""
Functions describing the constitutive response of hypo- or non-linear elastic materials.
"""
import numpy as np

def fun_void_ratio(void_ratio, constants=None, option="HardinRichart1963_v1"):
    """
    Function to calculate void ratio dependency of a material
    Parameters
    ----------
    void_ratio : float
        Void ratio of the material
    constants : list
        Material constants
    option : str
        Option for the calculation
    Returns
    -------
    fe : float
        Void ratio dependency of the material
    """
    if option == "HardinRichart1963_v1":
        # Hardin and Richart (1963)
        k = 2.17
        fe = ((k - void_ratio) ** 2) / (1 + void_ratio)
    elif option == "HardinRichart1963_v2":
        k = 2.97
        fe = ((k - void_ratio) ** 2) / (1 + void_ratio)
    else:
        raise NotImplementedError(f"Option {option} not implemented")
    return fe

def fun_mean_stress(mean_stress, constants=None, option="sand_default"):
    """
    Function to calculate mean stress dependency of a material
    Parameters
    ----------
    mean_stress : float
        Mean effective stress of the material
    constants : list
        Material constants
    option : str
        Option for the calculation
    Returns
    -------
    fp : float
        Mean stress dependency of the material
    """
    if option == "sand_default":
        fp = (mean_stress / 100) ** 0.5
    elif option == "clay_default":
        fp = (mean_stress / 100) ** 1.0
    else:
        raise NotImplementedError(f"Option {option} not implemented")
    return fp

def fun_ICG3S_R(eps, a, b, rmin, n=1.0):
    """
    Function to calculate the ICG3S-R function.
    Reduction ratio of any elastic modulus
    Parameters
    ----------
    eps : float
        Strain
    a : float
        Strain level at which the reduction ratio is 0.5
    b : float
        Slope of the reduction ratio curve in logarithmic strain scale
    rmin : float
        Minimum reduction ratio
    Returns
    -------
    r : float
        Reduction ratio
    """
    # calculate reduction ratio
    r = rmin + (1 - rmin) / (1 + (eps / (n * a)) ** b)
    return r

def fun_Gmax(Gref, void_ratio, mean_stress, vr_constants=None, p_constants=None, vr_option="HardinRichart1963_v1", p_option="sand_default"):
    """
    Function to calculate the maximum elastic shear modulus of a material.
    Parameters
    ----------
    Gref : float
        Reference shear modulus
    void_ratio : float
        Void ratio of the material
    mean_stress : float
        Mean effective stress of the material
    vr_constants : list
        Material constants for void ratio dependency
    p_constants : list
        Material constants for mean stress dependency
    vr_option : str
        Option for void ratio dependency calculation
    p_option : str
        Option for mean stress dependency calculation
    Returns
    -------
    Gmax : float
        Maximum elastic shear modulus of the material
    """
    fe = fun_void_ratio(void_ratio, constants=vr_constants, option=vr_option)
    fp = fun_mean_stress(mean_stress, constants=p_constants, option=p_option)
    Gmax = Gref * fe * fp
    return Gmax