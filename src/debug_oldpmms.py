# for debug only
import numpy as np
import pandas as pd
import os
import json
rho = 1025
Lpp = 340.5
beam =53.5
meanDraft = 17.0
underWaterLateralArea = 0.9124 * Lpp * meanDraft
SeparationPoint = 45.0
serviceSpeed = 10.8548
from scipy.constants import g


if os.path.isfile(r'data\PMMdata\acoef.json'):
    with open(r'data\PMMdata\acoef.json', 'r') as file:
        acoef = json.load(file)
    with open(r'data\PMMdata\bcoef.json', 'r') as file:
        bcoef = json.load(file)
        
def getForceCoefficient(multiPlier,force,utot2):
    '''
    returns forceCoefficient for a given force and speed
    multiPlier: dict of multipliers in a given baseTable .eg. mulitipler={'WHRD2':512.5,'UR_D':26.460736,'S':...}
    force: actual force calculated by pmm coefs
    utot2: squared speed
    '''
    factor = 1.0
    for key in multiPlier.keys():
        factor *= multiPlier[key]
    return force*utot2/factor          
    
def speedFactor(coef,U,V,R):
    '''
    based on standard PMM notation e.g UIVI read as U abs V
    the function calculate the speed factor a give coeeficients needs to be multiplied with
    so the pmmecoefficeint XUIVI with say value 123 then needs to be multiplied with U and abs(U) 
    '''
    st = list(coef)
    factor = 1
    absval = False
    istart = False
    iend = False
    for val in st:
      if val in ['X','Y','N']:
        continue # SKIP the Force 'direction'
      if val == 'I' and not istart:
        istart = True
        absval = True
        continue
      if val == 'I' and istart:
        istart = False
        continue
      if val =='0':
        factor *= 1
      if val == 'U':
        if absval:
          factor *=np.abs(U)
        else:
          factor *=U
      if val == 'V':
        if absval:
          factor *=np.abs(V)
        else:
          factor *=V
      if val == 'R':
        if absval:
          factor *=np.abs(R)
        else:
          factor *=R   
    return factor         

def pmmfor(pmmCoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh):
    '''
    this function calculates the forces on the hull using the PMM coefficient/hydrodynamic derivatives
    the function is translate to python by BTJ originally part of shipYard1 see tfs
    $/SimFlex Classic/src/lib/core/msdat/lib/pmm/pmmfor.f
    '''

    fdim=np.array([0.,0.,0.])
    fdimu2=np.array([0.,0.,0.])
    ## SKIPPING D,Q,P derivatives as we dont have them
    xderivatives=['X0','XU','XUU','XUUU','XV','XVV','XR','XRR','XVR','XUDOT','XRDOT','XUIUI']
    yderivatives=['Y0','Y0U','YUU','YV','YVV','YVIVI','YVVV','YVU','YR','YRR','YRIRI','YRRR','YRIVI','YRVV','YVIRI','YVRR','YVDOT','YRDOT']
    nderivatives=['N0','N0U','NUU','NV','NVV','NVIVI','NVVV','NVU','NR','NRR','NRIRI','NRRR','NRIVI','NRVV','NVIRI','NVRR','NVDOT','NRDOT ']
    lpp = Lpp

    if coftyp == 1:
        utot2 = udim*udim + vdim*vdim
        utot  = np.sqrt(utot2)
        U   = (udim - uo)/utot
        V   = vdim/utot
        R   = rdim/(utot/lpp)
        Q  = qdim
        P   = pdim/(utot/lpp)
        D   = ddim
        usign = 1.0
    else:
        utot2 = 1.
        utot  = 1.0
        U   = udim
        V   = vdim
        R   = rdim*lpp
        Q   = qdim
        P   = pdim*lpp
        D  = ddim
    #     if vdim and udim is zero we jump over the next because when gamma is -90 or 90
    #     this would lead to not a number (NAN) for BETAMSMFQ, which are used in the XVVSPEC coef.
        if (np.sqrt(vdim**2+udim**2) < 0.000001) :
            usign = np.sign(U)
            '''
            if U == 0:
                usign = 0
            elif U > 0:
                usign = 1
            elif U < 0:
                usign = -1
            '''
        else:
            BETAMSMFQ=np.sign(udim)*(-np.arcsin(vdim/np.sqrt(vdim**2+udim**2)))
            usign = np.sign(U)
            pass
    # X --- force
    xcoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    #inserting coefs and speed in arrays so we can do simple array multiplication to calculate force
    # check for key in dict needs to be added
    for val in xderivatives:
        if val in pmmCoefs.keys():
            if 'DOT' in val:
                continue  # in the fortran PMMFOR no DOT derivatives appears
            xcoeff=np.append(xcoeff,pmmCoefs[val])
            factor = speedFactor(val,U,V,R)
            speedvector=np.append(speedvector,factor)
            if val in acoef.keys():
                a = acoef[val]
                b = bcoef[val]
            else:
                a = 0
                b = 1
            TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))
        
    fdimu2[0]= np.sum(xcoeff * TOHCorrection * speedvector) * 0.5 * rho * lpp**2
    fdim[0] = fdimu2[0]*utot2
    # Y-force
    ycoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    for val in yderivatives:
        if val in pmmCoefs.keys():
            if 'DOT' in val:
                continue  # in the fortran PMMFOR no DOT derivatives appears
            ycoeff=np.append(ycoeff,pmmCoefs[val])
            speedvector=np.append(speedvector,speedFactor(val,U,V,R)) 
            if val in acoef.keys():
                a = acoef[val]
                b = bcoef[val]
            else:
                a = 0
                b = 1
            TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))       

    fdimu2[1]= np.sum(ycoeff * TOHCorrection * speedvector) * 0.5 * rho * lpp**2
    fdim[1] = fdimu2[1]*utot2    
    # N-moment
    ncoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    for val in nderivatives:
        if val in pmmCoefs.keys():
            if 'DOT' in val:
                continue  # in the fortran PMMFOR not DOT derivatives appears
            ncoeff=np.append(ncoeff,pmmCoefs[val])
            speedfactor = speedFactor(val,U,V,R)
            speedvector=np.append(speedvector,speedfactor) 
            # comment from fortran code
            #It is believed that N must change its sign for ship going astern
            #Therefore, USIGN is multiplied to the NvIvI
            if 'VIVI' in val:
                speedfactor *= usign
            if val in acoef.keys():
                a = acoef[val]
                b = bcoef[val]
            else:
                a = 0
                b = 1            
            TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b)*speedfactor)        

    fdimu2[2]= np.sum(ncoeff * TOHCorrection * speedvector) * 0.5 * rho * lpp**3
    fdim[2] = fdimu2[2]*utot2      

    return fdim,fdimu2,utot2  

def updateMultiplierWithActualSpeed(multiPlierDict,speed_dict):
    '''
    if a key from multiplierDict is found in speed_dic then 
    update the actual multiplier dict with the actual speed
    '''
    for key in multiPlierDict.keys():  
        if key in speed_dict.keys():
            multiPlierDict[key] = speed_dict[key]
    return multiPlierDict

def hluref(udim,vdim,rdim,pdim):
    '''
    calculated speeds used for dimensionalising
    '''

    lpp = Lpp
    b = beam
    dm = meanDraft

    speed_dict={}

    speed_dict['FN']    = udim / np.sqrt( lpp* g )
    speed_dict['GC']    = dm**2 + (b/2.)**2
    speed_dict['UR']    = udim**2
    speed_dict['UR_R']  = udim**2 + pdim**2*speed_dict['GC']
    speed_dict['UR_H']  = udim**2
    speed_dict['UR_D']  = udim**2 + vdim**2
    speed_dict['UR_Y']  = udim**2 +           (rdim*lpp/2.)**2
    speed_dict['UR_YD'] = udim**2 + vdim**2 + (rdim*lpp/2.)**2
    speed_dict['UR_DH'] = speed_dict['UR_D']
    speed_dict['UR_YH'] = speed_dict['UR_Y']
    return speed_dict

def pmmcar(uo,betad,coftyp):
    '''
    BTJ assume pmmcar is short for pmm carriage
    uo      : speed of carriage
    betad   : drift angle
    coftyp  : coefficient type 1 -> pmm  else msm
    
    '''
    if coftyp == 1 :
      ucar = uo/np.cos(betad)
    else:
      ucar = uo
    return ucar
def pmmmot(ucar,betad,gamma,delta,heel,epsil):
    '''
    BTJ: think the pmmmot is short for pmm motion
        ucar could be speed of carriage in the towingtank
    the routine returns 6 dof dimensional speeds depending on speed and driftangle,turnrate..
    '''
    lpp = Lpp
    udim  = ucar * np.cos(-betad)
    vdim  = ucar * np.sin(-betad)
    if (np.abs(np.abs(gamma)-1.570796) < .01) :
        rdim = 1.
        udim = 0.
        vdim = 0.
        if (ucar != 0):
          rdim = 2*ucar/lpp
        if (gamma < 0): 
          rdim = -np.abs(rdim)
    else:
        rdim  = ucar/lpp * 2. * np.tan(gamma)
  
    qdim  = heel
  #Comment : To be modified later ! ! ! ! ! ! ! ! ! ! ! !
    pdim  = epsil
    ddim  = delta
    return udim,vdim,rdim,qdim,pdim,ddim
def pmmmsm(uo,betad,gamma,pmmtyp,SepPoint):
    '''
    translated (modified and cleaned) to python by BTJ originaly inside tfs
      $/SimFlex Classic/src/lib/core/msdat/lib/pmm/pmms.f
    
    the functions determines the coefficient type depending on separationpoint SepPoint  
    
    cbla     If IMO maneouvres are most important then SepPoint should be 60.5
    cbla     If a 4 quadrant ship is most important then SepPoint should be 45.5
    cbla     The seppoint has to be sligthly over the value for which you want 
    cbla     pmm derivatives to be used.  
    '''
    coftyp = pmmtyp
    if np.abs(betad)-np.radians(SepPoint) > 0.00001 or  \
        np.abs(gamma)-np.radians(SepPoint) > 0.00001:
          coftyp=3

    if np.abs(betad) >  np.radians(185.0)-np.radians(SepPoint) or  \
        np.abs(gamma) >  np.radians(185.0)-np.radians(SepPoint):
        coftyp=2

    ucar = uo
    if np.abs(np.abs(gamma)-np.radians(90.)) < .01:
        ucar = 0.
    return coftyp,ucar   

def pmms(pmmcoefs,motion,icoty,uo,ipmms):
    '''
    the function is highly modified compared to original fortran version
    the purpose is to convert coefficient model to table model
    c              |                icoty
    c       -------|----------------------------------------
    c       motion |       0           |       1
    c       -------|-------------------|--------------------
    c         1    |    beta           |    beta,toh
    c         2    |    gamma          |    gamma,toh
    c         3    |    beta,gamma     |    toh
    c         4    |    phi            |    phi,toh
    c         5    |    beta,phi       |    toh
    c         6    |    gamma,phi      |    toh
    c         8    |    epsi           |    ----
    c       -------|-------------------|--------------------
    '''

   
    if motion == 1:
        absc1=[-180.0,-177.0,-175.0,-170.0,-135.0,-90.0,
                -70.0,-45.0,-20.0,-18.0,-15.0,-10.0,-5.0,-2.0,-1.0,
                    0.0,1.0,2.0,5.0,10.0,15.0,18.0,20.0,45.0,70.0,90.0,
                    135.0,170.0,175.0,177.0,180.0]
    if motion == 2:
        absc1=[-90.,-26.,-24.,-22.,-20.,-15.,-10., -5., -2., -1., 
              0.0,  1.,  2.,  5., 10., 15., 20., 22., 24., 26., 90.]
    absc2 = None

    MultiPliers={'DRIFT':{'X_HL':{'WHRD2':rho/2.,'UR_D':1,'ALW':underWaterLateralArea},
                           'Y_HL':{'WHRD2':rho/2.,'UR_D':1,'ALW':underWaterLateralArea},
                           'K_HL':{'WHRD2':rho/2.,'UR_D':1,'ALW':underWaterLateralArea,'DM':meanDraft},
                           'N_HL':{'WHRD2':rho/2.,'UR_D':1,'ALW':underWaterLateralArea,'LPP':Lpp}
                          },
                          'YAW':{'X_HL':{'WHRD2':rho/2.,'UR_Y':1,'ALW':underWaterLateralArea},
                           'Y_HL':{'WHRD2':rho/2.,'UR_Y':1,'ALW':underWaterLateralArea},
                           'N_HL':{'WHRD2':rho/2.,'UR_Y':1,'ALW':underWaterLateralArea,'LPP':Lpp}
                          }}
    if motion == 1: # Drift
      X_HL_multiPliers = MultiPliers['DRIFT']['X_HL']
      Y_HL_multiPliers = MultiPliers['DRIFT']['Y_HL']
      N_HL_multiPliers = MultiPliers['DRIFT']['N_HL']
      X_HL_Drift ={}
      Y_HL_Drift ={}
      N_HL_Drift ={}
      if absc2 != None:
        correctionTableX_HL = np.zeros((len(absc1),len(absc2)))
        correctionTableY_HL = np.zeros((len(absc1),len(absc2)))
        correctionTableN_HL = np.zeros((len(absc1),len(absc2)))
      for ixbeta,betad in enumerate(absc1):
        dummy = betad
        if dummy == -135.:
            idum =0
        betad = np.radians(betad)
        gamma = 0
        pmmtyp = 1
        SepPoint = SeparationPoint
        coftyp,ucar = pmmmsm(uo,betad,gamma,pmmtyp,SepPoint)
        uoo = ucar
        ucar = pmmcar(uoo,betad,coftyp)
        gamma = delta = heel = epsil = 0.0  

        udim,vdim,rdim,qdim,pdim,ddim = pmmmot(ucar,betad,gamma,delta,heel,epsil)
        speed_dict = hluref(udim,vdim,rdim,pdim)
        # update the multipliers with speed values if they are available
        X_HL_multiPliers = updateMultiplierWithActualSpeed(X_HL_multiPliers,speed_dict)
        Y_HL_multiPliers = updateMultiplierWithActualSpeed(Y_HL_multiPliers,speed_dict)        
        N_HL_multiPliers = updateMultiplierWithActualSpeed(N_HL_multiPliers,speed_dict)

        toh=0
        uo = serviceSpeed
        fdim,fdimu2,utot2 = pmmfor(pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh)
       
        if icoty == 0: #BaseTable
          X_HL = getForceCoefficient(X_HL_multiPliers,fdimu2[0],utot2)
          Y_HL = getForceCoefficient(Y_HL_multiPliers,fdimu2[1],utot2)
          N_HL = getForceCoefficient(N_HL_multiPliers,fdimu2[2],utot2)
    
          X_HL_Drift[dummy] = X_HL
          if np.isclose(np.abs(dummy),90.0) :
            X_HL_Drift[dummy] = 0.0
          if np.abs(dummy) >= 170:
              print(" Dirty hardcoding as I c'ant figure out where the sign change in SY1")
              X_HL_Drift[dummy] *= -1
          
          Y_HL_Drift[dummy] = Y_HL
          N_HL_Drift[dummy] = N_HL
            
        else: ## now we create a correction table
          FBASE = fdimu2
          for ix,toh in enumerate(absc2):
            if ix == 0:
              correctionTableX_HL[ixbeta,ix] = 1.0
              correctionTableY_HL[ixbeta,ix] = 1.0
              correctionTableN_HL[ixbeta,ix] = 1.0
              continue
            fdim,fdimu2,utot2 = pmmfor(pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh)
            correctionTableX_HL[ixbeta,ix] = fdimu2[0]/FBASE[0]
            correctionTableY_HL[ixbeta,ix] = fdimu2[1]/FBASE[1]
            correctionTableN_HL[ixbeta,ix] = fdimu2[2]/FBASE[2]
          pass # just used for debug after looping over toh's
      if icoty == 0:  # construct baseTables
        driftTables = {}
        dfx = pd.DataFrame(index=X_HL_Drift.keys(),data=X_HL_Drift.values())
        dfx.index.name = 'BETAD'
        dfx = dfx.rename(columns={0:'X_HL'})
      
        # These comments and the following code transfered and translated from fortran code    
        #cbla    the next if block are made to make the X_HL for drift look right,
        #cbla    meaning that we multiply the value at 20 degree driftangle with 1.5 
        #cbla    to get the value at 45 degree and make the table symmetrical. 
        #cbla    The 1.5 factor is taken from ship3005. 
          
        #index135 = np.where(np.abs(df['BETAD'].values) == 135.0)   
        #index70 =  np.where(np.abs(df['BETAD'].values) ==  70.0) 
        #index45 =  np.where(np.abs(df['BETAD'].values) ==  45.0)     
        #index20 =  np.where(np.abs(df['BETAD'].values) ==  20.0) 
        dfx.loc[-135.0]  = -1.5* dfx.loc[-20.0]
        dfx.loc[-70.0]   = dfx.loc[-20.0]
        dfx.loc[-45.0]   = 1.5* dfx.loc[-20.0]
        dfx.loc[ 45.0]   = 1.5* dfx.loc[-20.0]
        dfx.loc[ 70.0]   = dfx.loc[-20.0]
        dfx.loc[ 135.0]  = -1.5* dfx.loc[-20.0]
        
        dfy = pd.DataFrame(index=Y_HL_Drift.keys(),data=Y_HL_Drift.values())
        dfy.index.name = 'BETAD'
        dfy = dfy.rename(columns={0:'Y_HL'})
    
        dfn = pd.DataFrame(index=Y_HL_Drift.keys(),data=N_HL_Drift.values())
        dfn.index.name = 'BETAD'
        dfn = dfn.rename(columns={0:'N_HL'})    

        driftTables['X_HL = DRIFT'] = dfx
        driftTables['Y_HL = DRIFT'] = dfy
        driftTables['N_HL = DRIFT'] = dfn    
      else: # construct correction tables
        driftTables = {}
        dfx = pd.DataFrame(correctionTableX_HL,index=absc1,columns=absc2)
        dfy = pd.DataFrame(correctionTableY_HL,index=absc1,columns=absc2)
        dfn = pd.DataFrame(correctionTableN_HL,index=absc1,columns=absc2)
        driftTables['X_HL(BETAD,TOH)']=dfx
        driftTables['Y_HL(BETAD,TOH)']=dfy
        driftTables['N_HL(BETAD,TOH)']=dfy
      pass # end if motion == 1 aka DRIFT
      return driftTables
    elif motion == 2: # YAW
      X_HL_multiPliers = MultiPliers['YAW']['X_HL']
      Y_HL_multiPliers = MultiPliers['YAW']['Y_HL']
      N_HL_multiPliers = MultiPliers['YAW']['N_HL']
      X_HL_Yaw ={}
      Y_HL_Yaw ={}
      N_HL_Yaw ={}
      if absc2 != None:
        correctionTableX_HL = np.zeros((len(absc1),len(absc2)))
        correctionTableY_HL = np.zeros((len(absc1),len(absc2)))
        correctionTableN_HL = np.zeros((len(absc1),len(absc2)))
      for ixgamma,gamma in enumerate(absc1):
        dummy = gamma
        gamma = np.radians(gamma)
        betad = 0
        pmmtyp = 1
        SepPoint = SeparationPoint
        coftyp,ucar = pmmmsm(uo,betad,gamma,pmmtyp,SepPoint)
        
        uoo = ucar
        ucar = pmmcar(uoo,betad,coftyp)
        delta = heel = epsil = 0.0  

        udim,vdim,rdim,qdim,pdim,ddim = pmmmot(ucar,betad,gamma,delta,heel,epsil)
        speed_dict = hluref(udim,vdim,rdim,pdim)
        # update the multipliers with speed values if they are available
        X_HL_multiPliers = updateMultiplierWithActualSpeed(X_HL_multiPliers,speed_dict)
        Y_HL_multiPliers = updateMultiplierWithActualSpeed(Y_HL_multiPliers,speed_dict)        
        N_HL_multiPliers = updateMultiplierWithActualSpeed(N_HL_multiPliers,speed_dict)

        toh=0
        uo = serviceSpeed
        fdim,fdimu2,utot2 = pmmfor(pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh)
        if np.abs(dummy) == 90.0 :
          print('Can not see how this works in the fortan code but Y and X (gamma) seems to be zero so hardcoded here !!')
          fdimu2[0] = 0.0
          fdimu2[1] = 0.0
       
        if icoty == 0: #BaseTable
          X_HL = getForceCoefficient(X_HL_multiPliers,fdimu2[0],utot2)
          Y_HL = getForceCoefficient(Y_HL_multiPliers,fdimu2[1],utot2)
          N_HL = getForceCoefficient(N_HL_multiPliers,fdimu2[2],utot2)
    
          X_HL_Yaw[dummy] = X_HL       
          Y_HL_Yaw[dummy] = Y_HL
          N_HL_Yaw[dummy] = N_HL
            
        else: ## now we create a correction table
          FBASE = fdimu2
          if np.abs(dummy) == 90.0 :
            print (' need to set FBASE to avoid div 0')
            FBASE[0] = FBASE[1] = 1.0
          for ix,toh in enumerate(absc2):
            if ix == 0:
              correctionTableX_HL[ixgamma,ix] = 1.0
              correctionTableY_HL[ixgamma,ix] = 1.0
              correctionTableN_HL[ixgamma,ix] = 1.0
              continue
            fdim,fdimu2,utot2 = pmmfor(pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh)
            correctionTableX_HL[ixgamma,ix] = fdimu2[0]/FBASE[0]
            correctionTableY_HL[ixgamma,ix] = fdimu2[1]/FBASE[1]
            correctionTableN_HL[ixgamma,ix] = fdimu2[2]/FBASE[2]
          pass # just used for debug after looping over toh's
      if icoty == 0:  # construct baseTables
        yawTables = {}
        dfx = pd.DataFrame(index=X_HL_Yaw.keys(),data=X_HL_Yaw.values())
        dfx.index.name = 'GAMMA'
        dfx = dfx.rename(columns={0:'X_HL'})
        
        dfy = pd.DataFrame(index=Y_HL_Yaw.keys(),data=Y_HL_Yaw.values())
        dfy.index.name = 'GAMMA'
        dfy = dfy.rename(columns={0:'Y_HL'})
    
        dfn = pd.DataFrame(index=Y_HL_Yaw.keys(),data=N_HL_Yaw.values())
        dfn.index.name = 'GAMMA'
        dfn = dfn.rename(columns={0:'N_HL'})    

        yawTables['X_HL = YAW'] = dfx
        yawTables['Y_HL = YAW'] = dfy
        yawTables['N_HL = YAW'] = dfn    
      else: # construct correction tables
        yawTables = {}
        dfx = pd.DataFrame(correctionTableX_HL,index=absc1,columns=absc2)
        dfy = pd.DataFrame(correctionTableY_HL,index=absc1,columns=absc2)
        dfn = pd.DataFrame(correctionTableN_HL,index=absc1,columns=absc2)
        yawTables['X_HL(GAMMA,TOH)']=dfx
        yawTables['Y_HL(GAMMA,TOH)']=dfy
        yawTables['N_HL(GAMMA,TOH)']=dfy
      pass # end if motion == 1 aka DRIFT
      return yawTables

if __name__ == '__main__':
    with open(r'data\PMMdata\3949NeuPMMDe_005_001.txt','r') as f:
        pmmCoefs = {}
        for line in f:
            splt = line.split()
            pmmCoefs[splt[0]] = float(splt[1])/1.0E5  
    motion = 2
    icoty = 0
    ipmms = 1
    uo = serviceSpeed
    table = pmms(pmmCoefs,motion,icoty,uo,ipmms)
    pass