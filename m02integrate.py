# txc-integrator-M02-M09
#
# This script integrates IC MAGE M02 & M09 for drained constant-p triaxial compression tests.
#
# Reference:
# Taborda & Pirrone (2021, September 22) TXC-Integrator for IC MAGE M02 & M09 (Version 1.0). Zenodo. https://doi.org/10.5281/zenodo.5521316.
#
# --------------------------------------------------------------
# Date       | Version | Log changes
# --------------------------------------------------------------
# 01/10/2021 | 1.0     | DMG Taborda | Initial version for TXC
# --------------------------------------------------------------
#
# Import packages
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
#
def fun_m02_txc(material_parameters, initial_conditions):
    """
    Function to integrate IC MAGE M02 for drained constant-p triaxial compression tests.
    
    Parameters
    ----------
    material_parameters : dict
        Dictionary containing material parameters for the model.
    initial_conditions : dict
        Dictionary containing initial conditions for the model.
    
    Returns
    -------
    output : dict
        Dictionary containing the results of the integration.

    """
    # Define parameters
    ecsref = material_parameters['e_csref']
    l = material_parameters['l_cs']
    csi = material_parameters['csi_cs']
    pref = material_parameters['pref']
    mcs = material_parameters['mcs']
    k1 = material_parameters['k1']
    k2 = material_parameters['k2']
    l1 = material_parameters['l1']
    l2 = material_parameters['l2']

    # Small strain stiffness model parameters (see Taborda et al. 2016 and model manuals)
    g0    = material_parameters['g0']
    mg    = material_parameters['mg']
    a0    = material_parameters['a0']
    b     = material_parameters['b']
    rgmin = material_parameters['rgmin']
    pr    = material_parameters['pr']

    # Initial conditions
    e0 = initial_conditions['e0']
    p0 = initial_conditions['p0']
    q0 = initial_conditions['q0']

    psi0  = e0 - (ecsref-l*(p0/pref)**csi)
    if psi0 > 0:
        mt = mcs
    else:
        mt    = mcs + k1*((-psi0)**(k2+1.0))

    ppeak = p0
    for i in range(0,100):
        ppeak = ppeak/3.0*mcs+k1/3.0*ppeak*(ecsref - l*(ppeak/pref)**csi-e0)**(k2+1.0)+p0-q0/3.0

    qpeak = ppeak*mt

    pcs = (3.0*p0-q0)/(3.0 - mcs)
    qcs = pcs*mcs
    #
    eax  = []
    erad = []
    evol = []
    ed   = []
    p    = []
    q    = []
    vr   = []
    pvr  = []   # plastic void ratio
    psi  = []
    mtxc = []
    mdil = []
    dp   = []
    dq   = []
    fy   = []
    gtan = []
    ktan = []
    #
    # Initialise quantities
    eax.append(0)
    erad.append(0)
    evol.append(0)
    ed.append(0)
    p.append(p0)
    q.append(q0)
    vr.append(e0)
    pvr.append(e0)
    #
    psi0 = e0 - (ecsref-l*(p0/pref)**csi)
    psi.append(psi0)
    #
    mt = mcs + k1*((-psi0)**(k2+1.0))
    mtxc.append(mt)
    #
    if psi0 > 0: # loose
        mt = 1e-6
    else:
        mt = mcs - l1*((-psi0)**(l2+1.0))
    mdil.append(mt)
    #
    gt = g0*((p0/pref)**mg)
    gtan.append(gt)
    kt = 2.0*gt*(1.0+pr)/(3.0*(1.0-2.0*pr))
    ktan.append(kt)
    #
    # Yield surface
    fy.append(q0 - p0*mtxc[0])
    #
    # Start elastic integration
    for i in range(1, 1001):
        pstep = (ppeak - p0)/1000
        qstep = 3.0*pstep
        #
        devole = pstep/ktan[i-1]
        dede   = qstep/(math.sqrt(3)*gtan[i-1])
        #
        # update quantities
        x = evol[i-1]+devole
        evol.append(x)
        #
        x = ed[i-1]+dede
        ed.append(x)
        #
        x = eax[i-1]+devole/3.0+dede/math.sqrt(3)
        eax.append(x)
        #
        x = erad[i-1]+devole/3+dede/(2.0*math.sqrt(3))
        erad.append(x)
        #
        x = p[i-1]+pstep
        p.append(x)
        #
        x = q[i-1]+qstep
        q.append(x)
        #
        x = vr[i-1]-(1+e0)*devole
        vr.append(x)
        #
        pvr.append(pvr[i-1])
        #
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#    
        # Preparation for next step
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        x = pvr[i] - (ecsref-l*(p[i]/pref)**csi)
        #
        psi.append(x)
        #
        x = g0*((p[i]/pref)**mg)*(rgmin + (1.0-rgmin)/(1.0+(ed[i]/a0)**b))
        gtan.append(x)
        #
        kt = 2.0*x*(1.0+pr)/(3.0*(1.0-2.0*pr))
        ktan.append(kt)
        #
        mt = mcs + k1*((-psi[i])**(k2+1.0))
        mtxc.append(mt)
        #
        if psi[i] > 0: # loose
            mt = 1e-9
        else:
            mt = l1*((-psi[i])**(l2+1.0))
        mdil.append(mt)    
        #
        fy.append(q[i] - p[i]*mtxc[i])
        #
    #
    # Start plastic integration
    for i in range(1001, 2001):
        pstep = (pcs - ppeak)/1000
        qstep = 3.0*pstep
        #
        devole = pstep/ktan[i-1]
        dede   = qstep/(math.sqrt(3)*gtan[i-1])
        #
        dfdp   = -mtxc[i-1]+p[i-1]*k1*(k2+1.0)*((-psi[i-1])**k2)*(l*csi/pref)*((p[i-1]/pref)**(csi-1.0))
        dfdq   = 1.0
        dpdp   = -mdil[i-1]
        dpdq   = 1.0
        dfdep  = p[i-1]*k1*(k2+1.0)*((-psi[i-1])**k2)
        Qval   = dfdp*dpdp*ktan[i-1]+3.0*dfdq*dpdq*gtan[i-1]
        Rval   = Qval + dfdep*(1+e0)*dpdp
        D11    = ktan[i-1]-ktan[i-1]*ktan[i-1]*dfdp*dpdp/Rval
        D12    = -3.0*ktan[i-1]*gtan[i-1]*dfdq*dpdp/Rval
        D21    = -3.0*ktan[i-1]*gtan[i-1]*dfdp*dpdq/Rval
        D22    =  3.0*gtan[i-1]-(9.0*gtan[i-1]*gtan[i-1])*dfdq*dpdq/Rval
        detD   = D11*D22 - D12*D21
        Di11   = D22/detD
        Di12   = -D12/detD
        Di21   = -D21/detD
        Di22   = D11/detD
        devol  = Di11*pstep+Di12*qstep
        ded    = math.sqrt(3)*(Di21*pstep+Di22*qstep)   # Multiply by sqrt(3) to convert to Ed
        devolp = devol - devole
        dedp   = ded   - dede
        dpvr   = -(1+e0)*devolp
        #
        # Update quantities
        x = evol[i-1]+devol
        evol.append(x)
        #
        x = ed[i-1]+ded
        ed.append(x)
        #
        x = eax[i-1]+devol/3.0+ded/math.sqrt(3.0)
        eax.append(x)
        #
        x = erad[i-1]+devol/3.0+ded/(2.0*math.sqrt(3.0))
        erad.append(x)
        #
        x = p[i-1]+pstep
        p.append(x)
        #
        x = q[i-1]+qstep
        q.append(x)
        #
        x = vr[i-1]-(1+e0)*devol
        vr.append(x)
        #
        x = pvr[i-1]+dpvr
        pvr.append(x)
        #
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#    
        # Preparation for next step
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        #
        x = pvr[i] - (ecsref-l*(p[i]/pref)**csi)
        # print(pvr[i-1], devol, devolp, dpvr, pvr[i], p[i])
        #
        if x>-1E-6:
            x = -1E-6
        #
        psi.append(x)
        #
        x = g0*((p[i]/pref)**mg)*(rgmin + (1.0-rgmin)/(1.0+(ed[i]/a0)**b))
        gtan.append(x)
        #
        kt = 2.0*x*(1.0+pr)/(3.0*(1.0-2.0*pr))
        ktan.append(kt)
        #
        mt = mcs + k1*((-psi[i])**(k2+1.0))
        mtxc.append(mt)
        #
        mt = l1*((-psi[i])**(l2+1.0))
        mdil.append(mt)    
        #
        fy.append(q[i] - p[i]*mtxc[i])    

    # for i in range(0, 2001):
    #     eax[i] = eax[i]*100.0
    #     evol[i] = evol[i]*100.0
    # #
    # df = pd.DataFrame(columns=['eax0%','evol%','eax1%','q0','p0','q1','p1','e','Ed0','q2','Ed1','Gtan'])
    # df['eax0%'] = eax
    # df['evol%'] = evol
    # df['eax1%'] = eax
    # df['q0'] = q
    # df['p0'] = p
    # df['q1'] = q
    # df['p1'] = p
    # df['e'] = vr
    # df['Ed0'] = ed
    # df['q2'] = q
    # df['Ed1'] = ed
    # df['Gtan'] = gtan
    output = {
        'eax_prc': np.array(eax) * 100,
        'evol_prc': np.array(evol) * 100,
        'erad_prc': np.array(erad) * 100,
        'ed_prc': np.array(ed) * 100,
        'p': np.array(p),
        'q': np.array(q),
        'vr': np.array(vr),
        'pvr': np.array(pvr),
        'psi': np.array(psi),
        'mtxc': np.array(mtxc),
        'mdil': np.array(mdil),
        'dpvr': np.array(dpvr),
        'gtan': np.array(gtan),
        'ktan': np.array(ktan),
    }
    return output