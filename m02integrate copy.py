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
import math
import matplotlib.pyplot as plt
import pandas as pd
#
# Define the number of tests to integrate when running
ntest = 1
#
# The variable path determines the location of where to write the output files. 
# The name of the file (defined in labels) will be added to this path, followed by the extension csv.
path = "C:\\Users\\dtaborda\\OneDrive - Imperial College London\\Documents\\Github\\IC-MAGE-M02-tools\\M02_"
path = "D:\\Users\\dmgta\\OneDrive - Imperial College London\\Documents\\Github\\IC-MAGE-M02-tools\\M02_"
#
# Parameter sets for Nevada sand. New sets can be added or introduced by changing those below.
# Parameter names can be found in IC MAGE M02 and M09 manuals.
set = 1
if set == 1:
    # Nevada sand (original formulation)
    ecsref = 0.887
    l      = 0.079
    csi    = 0.250
    pref   = 100.0
    #
    mcs    = 1.2872
    #
    k1 = 2.4
    k2 = 0.0
    l1 = 3.7
    l2 = 0.0
elif set == 2:
    # Nevada sand (dr-based formulation)
    ecsref = 0.774
    l      = 0.0
    csi    = 1.0
    pref   = 100.0
    #
    mcs    = 1.2872
    #
    k1 = 1.8
    k2 = -0.3
    l1 = 2.4
    l2 = -0.35
#
# Small strain stiffness model parameters (see Taborda et al. 2016 and model manuals)
g0    = 51860.0
mg    = 0.5
a0    = 1.09E-4
b     = 1.13
rgmin = 0.1
pr    = 0.2
#
# Initial conditions
labels = ['CIDC40-107', 'CIDC40-100', 'CIDC40-106', 'CIDC60-82', 'CIDC60-75', 'CIDC60-81', 'CADC40-108', 'CADC60-70']
e0 = [0.728, 0.726, 0.718, 0.657, 0.652, 0.651, 0.723, 0.653]
p0 = [ 40.0,  80.0, 160.0,  40.0,  80.0, 160.0,  80.0,  80.0]
q0 = [  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  44.0,  44.0]
#
for itest in range(0,ntest):
    # Counter
    print('*** Processing test:'+str(itest+1))
    # Peak conditions
    psi0  = e0[itest] - (ecsref-l*(p0[itest]/pref)**csi)
    mt    = mcs + k1*((-psi0)**(k2+1.0))
    #
    ppeak = p0[itest]
    for i in range(0,100):
        ppeak = ppeak/3.0*mcs+k1/3.0*ppeak*(ecsref - l*(ppeak/pref)**csi-e0[itest])**(k2+1.0)+p0[itest]-q0[itest]/3.0
    #
    qpeak = ppeak*mt
    #
    # Critical state
    pcs = (3.0*p0[itest]-q0[itest])/(3.0 - mcs)
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
    p.append(p0[itest])
    q.append(q0[itest])
    vr.append(e0[itest])
    pvr.append(e0[itest])
    #
    psi0 = e0[itest] - (ecsref-l*(p0[itest]/pref)**csi)
    psi.append(psi0)
    #
    mt = mcs + k1*((-psi0)**(k2+1.0))
    mtxc.append(mt)
    #
    mt = mcs - l1*((-psi0)**(l2+1.0))
    mdil.append(mt)
    #
    gt = g0*((p0[itest]/pref)**mg)
    gtan.append(gt)
    kt = 2.0*gt*(1.0+pr)/(3.0*(1.0-2.0*pr))
    ktan.append(kt)
    #
    # Yield surface
    fy.append(q0[itest] - p0[itest]*mtxc[0])
    #
    # Start elastic integration
    for i in range(1, 1001):
        pstep = (ppeak - p0[itest])/1000
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
        x = vr[i-1]-(1+e0[itest])*devole
        vr.append(x)
        #
        pvr.append(pvr[i-1])
        #
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#    
        # Preparation for next step
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        #
        x = pvr[i] - (ecsref-l*(p[i]/pref)**csi)
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
        Rval   = Qval + dfdep*(1+e0[itest])*dpdp
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
        dpvr   = -(1+e0[itest])*devolp
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
        x = vr[i-1]-(1+e0[itest])*devol
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

    #plt.plot(eax,evol)
    #
    #plt.show()
    #
    #plt.plot(eax,q)
    #
    #plt.show()
    #
    for i in range(0, 2001):
        eax[i] = eax[i]*100.0
        evol[i] = evol[i]*100.0
    #
    df = pd.DataFrame(columns=['eax0%','evol%','eax1%','q0','p0','q1','p1','e','Ed0','q2','Ed1','Gtan'])
    df['eax0%'] = eax
    df['evol%'] = evol
    df['eax1%'] = eax
    df['q0'] = q
    df['p0'] = p
    df['q1'] = q
    df['p1'] = p
    df['e'] = vr
    df['Ed0'] = ed
    df['q2'] = q
    df['Ed1'] = ed
    df['Gtan'] = gtan
    #
    path_out = path+labels[itest]+'.csv'
    df.to_csv(path_out)
    #