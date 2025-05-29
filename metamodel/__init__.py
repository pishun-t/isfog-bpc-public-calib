import pickle
import numpy as np
import os

### --- Metamodel --- ###
model_folder = r"metamodel"
with open(os.path.join(model_folder, r"pce_model_fitted_2_param.pkl"), 'rb') as f:
    pce_model_2_param = pickle.load(f)

def mm_output(ux_mm, Gref=65, a=5e-5, b=1.1, K0=0.4598, intf_delta=24):
    input_params = np.array([Gref, np.log10(a), b, K0, intf_delta])
    p1, p2 = np.array(pce_model_2_param(input_params))
    Fx_kN = p1 * (ux_mm / ( 1 + ux_mm / p2) )
    return Fx_kN