

# %%
import numpy as np
from scipy.constants import g
from scipy.optimize import fmin_bfgs,minimize
import json
import os

import pandas as pd

import Thumbs as Thumbs


def FUNCAB(XAB,dfTOH,index):
  '''
  The following is as understod by code inspection of shipyard1 code and traslated to python by BTJ june 2024
  
  Shallowater correction of manoeuvrecoefficients are typically calculated as
  (1 + a*TOH**b)
  the function is evaluated over a range of a and b's and the (a,b) which gives the minimum FUNCV is used as 
  the correction factors for the given hydrodynamic derivative
  Parameters:
              XAB:          a length 2 array with a,b
            dfTOH:          dataframe with manoeuvre coefficient names as index (rows) and TOH as columns
                            each item in dfTOH is calculated by the dmimix function, which combines different methods
                            to find shallow water corrections for the different derivatives
            index:          the actual manoeuvring coefficient
            
              
  '''
  TOH = [1./10., 1./5., 1./2., 1./1.5, 1./1.2]
  # Initialize FUNCV
  FUNCV = 0
  NTOH = 5
  # Iterate over NTOH
  # Note: In Python, you would typically have NTOH as an input to the function, or it would be a global variable.
  # Here, I'm assuming it's a global variable for the sake of this translation.
  for I in range(NTOH):
      # Check the condition
      if XAB[1] > -87.0/np.log(TOH[I]):
          RDUMMY = 0.0
      else:
          RDUMMY = XAB[0]*TOH[I]**XAB[1]

      # Update FUNCV
      FUNCV += (RDUMMY - (dfTOH.loc[index,dfTOH.columns[I]] - 1.0))**2
      
  # Return FUNCV
  #print(f"FUNCAB {FUNCV}")
  return FUNCV 

class HullThumbs(Thumbs.Thumbs):
  '''
  class used to hold the different functions used in shipyard1 to create hull tables
  '''
    
  def __init__(self,shipDataDict):
    super().__init__(shipDataDict)
#
# FNH values for squat calculation
#
    self.FNH_values = [-1.90,-1.20,-1.10,-1.05,-1.00, 
                  -0.95,-0.90,-0.80,-0.70,-0.60,
                  -0.50,-0.40,-0.30,-0.20,-0.10, 0.00,
                    0.10, 0.20, 0.30, 0.40, 0.50,
                    0.60, 0.70, 0.80, 0.90, 0.95,
                    1.00, 1.05, 1.10, 1.20, 1.90]
      
    self.TOH=[1./10.,1./5.,1./2.,1./1.5,1./1.2]
    self.linearDamping()
    self.calculateSquat()
    # define the names of the hull force tables X,Y,N K not considered yet
    self.X_HLBaseTables=['DRIFT_HEEL','YAW_HEEL','ROLL','FN','DRIFT','YAW','YAW_DRIFT']
    self.Y_HLBaseTables=self.X_HLBaseTables.copy()
    self.Y_HLBaseTables.remove('FN')
    self.N_HLBaseTables= self.Y_HLBaseTables.copy()
    
    # for each drift table set the variable name and the scale factors the speed parameter will be updated later
    self.DriftMultiPliers={'X_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea},
                           'Y_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea},
                           'K_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea,'DM':self.meanDraft},
                           'N_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea,'LPP':self.Lpp}
                          }
    print('!!! To save time during development we read acoef and bcoef from files !!!')
    #the dmimix calculate a and b coefficients for shallow water correction
    #it is time consuming and for debug purpose they are for now read from files
    if os.path.isfile(r'data\PMMdata\acoef.json'):
      with open(r'data\PMMdata\acoef.json', 'r') as file:
        self.acoef = json.load(file)
      with open(r'data\PMMdata\bcoef.json', 'r') as file:
        self.bcoef = json.load(file)
    else:
      self.acoef,self.bcoef = self.dmimix()
      with open(r'data\PMMdata\acoef.json','w') as file:
        json.dump(self.acoef, file)
      with open(r'data\PMMdata\bcoef.json','w') as file:
        json.dump(self.bcoef,file)
    pass
    
    

  
  def pmmmot(self,ucar,betad,gamma,delta,heel,epsil):
    '''
    BTJ: think the pmmmot is short for pmm motion
        ucar could be speed of carriage in the towingtank
    the routine returns 6 dof dimensional speeds depending on speed and driftangle
    '''
    lpp = self.Lpp
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
      
  
  
  def pmms(self,pmmcoefs,motion,icoty,uoo,ipmms):
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
    absc1 , absc2 = self.defval(motion,icoty)
    
    if motion == 1: # Drift
      X_HL_multiPliers = self.DriftMultiPliers['X_HL']
      X_HL_Drift =[]
      for betad in absc1:
        betad = np.radians(betad)
        gamma = 0
        pmmtyp = 1
        SepPoint = np.radians(45.5)
        coftyp,ucar = self.pmmmsm(uoo,betad,gamma,pmmtyp,SepPoint)
        ucar = self.pmmcar(uoo,betad,coftyp)
        gamma = delta = heel = epsil = 0.0  
       

        udim,vdim,rdim,qdim,pdim,ddim = self.pmmmot(ucar,betad,gamma,delta,heel,epsil)
        speed_dict = self.hluref(udim,vdim,rdim,pdim)
        # update the multipliers with speed values if they are available
        for key in X_HL_multiPliers.keys():  
          if key in speed_dict.keys():
            X_HL_multiPliers[key] = speed_dict[key]
        toh=0
        uo = self.serviceSpeed
        fdim,fdimu2,utot2 = self.pmmfor(pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh)
        factor = 1
        if icoty == 0: #Base table
          for key in X_HL_multiPliers.keys():
            factor *= X_HL_multiPliers[key]
          val = fdimu2[0]*utot2/factor  
          X_HL_Drift.append(val)
        else:
          print ('Correction tables not implemented yet !!!')
    
          pass  
    d = {'BETAD': absc1, 'X_HL': np.array(X_HL_Drift)}
    df = pd.DataFrame(data=d)
    return df,X_HL_multiPliers.keys()
  

  def getBaseTable(self,pmmcoefs,tableName):
    '''
    returns a dict of baseTables given a set of pmmCoefficients
    
    '''
    baseTable = None
    if tableName == 'DRIFT':
      motion = 1
      icoty = 0
      ipmms = 1
      uo = self.serviceSpeed
      baseTable,multiPliers = self.pmms(pmmcoefs,motion,icoty,uo,ipmms)
    return baseTable,multiPliers
    

  def linearDamping(self):
    '''
    internal function to create linear damping tables
    '''
    RHO=self.rho
    grav = g
    L5 = self.Lpp
    B9 = self.beam
    TXAFT = self.draftAft
    TXFORE = self.draftFore
    FMFORE = TXFORE
    FMAFT = TXAFT
    T9 = self.meanDraft
    V0 = self.blockCoefficient   

    C2 = V0  
    V0 = self.displacement
    COG = self.COG
    ALW=self.waterPlaneArea
    GYRATION = self.gyration
    VCOB =self.verticalCenterOfBoyancy
    IXX = GYRATION[0]**2*V0*RHO
    IYY = GYRATION[1]**2*V0*RHO
    GMT = self.GMT
    GML = self.GML
    omerol = np.sqrt(grav*V0*RHO*GMT/IXX)
    omepit = np.sqrt(grav*V0*RHO*GML/IYY)
    omehea = np.sqrt(ALW*grav/V0)
    print (f'resonans periods, roll  {1/omerol*2*np.pi}, pitch {1/omepit*2*np.pi} heave {1/omehea*2*np.pi} ')

    #	Very rough assumptions are applied here. one section shape, 
    #	pitch equal to roll
    #	calculate the Ap coefficients (dp assumed = 0.5)

    DPROLL = 0.5
    DPPITCH = 0.1*B9/L5
    APROLL = DPROLL*(omerol**2*B9/(2.0*grav))**2
    APPITC = DPPITCH*(omepit**2*L5/(2.0*grav))**2
    APHEAV = 2*np.sin(omehea**3*B9/(2*grav))*np.exp(-omehea**2*T9/grav)
    BSECROLL = RHO*grav**2/omerol**3*(B9/2)**2*APROLL**2
    ROLLDAMP = BSECROLL*L5
    BSECPITC = RHO*grav**2/omepit**3*(L5/2)**2*APPITC**2
    PITCDAMP = BSECPITC*B9
    BSECHEAV = RHO*grav**2/omehea**3*APHEAV**2
    HEAVDAMP = BSECHEAV*L5
    print (f'roll , pitch, heave, ',ROLLDAMP, PITCDAMP, HEAVDAMP)


    #	NOW MAKE THE TABLES
    #------------------------------------------------------------------------
    #	HEAVE TABLES
    #------------------------------------------------------------------------


    self.HeaveDamping=np.zeros((3,2))
    self.HeaveDamping[0,0] = -10
    self.HeaveDamping[0,1] = 10*HEAVDAMP
    self.HeaveDamping[1,0] = 0
    self.HeaveDamping[1,1] = 0
    self.HeaveDamping[2,0] = 10
    self.HeaveDamping[2,1] = -10*HEAVDAMP


    #------------------------------------------------------------------------
    #	PITCH TABLES
    #------------------------------------------------------------------------
    self.PitchDamping=np.zeros((3,2))
    self.PitchDamping[0,0] = -1
    self.PitchDamping[0,1] = PITCDAMP
    self.PitchDamping[1,0] = 0
    self.PitchDamping[1,1] = 0
    self.PitchDamping[2,0] = 1
    self.PitchDamping[2,1] = -PITCDAMP


    #------------------------------------------------------------------------
    #	ROLL TABLES
    #------------------------------------------------------------------------
    self.RollDamping=np.zeros((3,2))
    self.RollDamping[0,0] = -1
    self.RollDamping[0,1] = ROLLDAMP
    self.RollDamping[1,0] = 0
    self.RollDamping[1,1] = 0
    self.RollDamping[2,0] = 1
    self.RollDamping[2,1] = -ROLLDAMP


    #     ----------------------------------------------------------------
    #	Now calculate the high frequency dampings in sway, yaw and surge
    #     ----------------------------------------------------------------


    self.SurgeDamping=np.zeros((3,2))
    self.SurgeDamping[0,0] = -10
    self.SurgeDamping[0,1] = 20
    self.SurgeDamping[1,0] = 0
    self.SurgeDamping[1,1] = 0
    self.SurgeDamping[2,0] = 10
    self.SurgeDamping[2,1] = -20

    self.SwayDamping=np.zeros((3,2))
    self.SwayDamping[0,0] = -10
    self.SwayDamping[0,1] = 0
    self.SwayDamping[1,0] = 0
    self.SwayDamping[1,1] = 0
    self.SwayDamping[2,0] = 10
    self.SwayDamping[2,1] = -10      
  
  def getSurgeDamping(self):
    return self.SurgeDamping
  def getSwayDamping(sel):
    return self.SwayDamping
  def getRollDamping(sel):
    return self.RollDamping
  def getPitchDamping(sel):
    return self.PitchDamping
  def getHeaveDamping(sel):
    return self.HeaveDamping     


  def S_mid(self,CB, T, L, FNH):
      '''
      helper function for squat calculation
      '''
      rVal = 0.01 * ((38 * CB * T / L) * FNH**2) / (np.sqrt(1 - FNH**2))
      return rVal

  def trim(self,CB, T, L, FNH):
      '''
      helper function for squat calculation
      '''    
      rVal = 0.01 * (47.4 * CB * T / L - 1.2) * FNH**2 / (np.sqrt(1 - FNH**2))
      return rVal


  def calculateSquat(self):
    '''
    calculate squat tables
    '''
    FNH_values = self.FNH_values
    CB = self.blockCoefficient
    T = self.meanDraft
    L = self.Lpp
    squaFPlist = []
    squaAPlist = []
    for fnh in FNH_values:
        if np.abs(fnh) == 1.9 or np.abs(fnh) == 1.2:
            squaFP = -np.sign(fnh)*np.tan(np.radians(1.4)) / 2.0
            squaAP = np.sign(fnh)*np.tan(np.radians(1.4)) / 2.0
        if np.abs(fnh) == 1.1:
            squaFP = -np.sign(fnh)*np.tan(np.radians(1.35)) / 2
            squaAP = np.sign(fnh)*np.tan(np.radians(1.35)) / 2
        if np.abs(fnh) == 1.05:
            squaFP = -self.S_mid(CB, T, L, np.sign(fnh)*0.9) / 3 -np.sign(fnh)*np.tan(np.radians(1.3)) / 2
            squaAP = -self.S_mid(CB, T, L, np.sign(fnh)*0.9) / 3 + np.sign(fnh)*np.tan(np.radians(1.3)) / 2
        if np.abs(fnh) == 1.0:
            squaFP = -np.sign(fnh)*np.tan(np.radians(2.0)) / 2
            squaAP = np.sign(fnh)*np.tan(np.radians(2.0)) / 2
        if np.abs(fnh) == 0.95:
            squaFP = 2 * self.S_mid(CB, T, L, np.sign(fnh)*0.9) / 3 - np.sign(fnh)*np.tan(np.radians(1.8)) / 2
            squaAP = 2 * self.S_mid(CB, T, L, np.sign(fnh)*0.9) / 3 + np.sign(fnh)*np.tan(np.radians(1.8)) / 2
        if fnh >= -0.9 and fnh <= 0.9:
            sign = np.sign(fnh)
            squaFP = self.S_mid(CB, T, L, fnh) + sign * 0.5 * self.trim(CB, T, L, fnh)
            squaAP = self.S_mid(CB, T, L, fnh) - sign * 0.5 * self.trim(CB, T, L, fnh)

        squaFPlist.append(squaFP)
        squaAPlist.append(squaAP)

    self.SquatFP=np.transpose(np.array((self.FNH_values,squaFPlist)))
    self.SquatAP=np.transpose(np.array((self.FNH_values,squaAPlist)))


  def getSquatFP(self):
    return self.SquatFP

  def getSquatAP(self):
    return self.SquatAP

  def addedMassTOH(self,tableName):
    '''
    This method returns the added mass TOH correction curves 
    according to Leif Wagner Smidth (LWS)
    original fortran code from ShipYard1 
    admtoh.f (added mass toh)
    '''
    if tableName in ['XDCOR','XHCOR','NDCOR','NHCOR']:
      acoef = 2.0
      bcoef = 3.0
    else:
      acoef = 4.0
      bcoef = 3.0
    xvals=np.arange(0,1.1,0.2)
    yvals =np.zeros(len(xvals))
    for ix in range(len(xvals)):
      yvals[ix]= 1 + acoef*xvals[ix]**bcoef
    res = np.array(list(zip(xvals,yvals)))
    return res
    
  def speedfactor(self,coef,U,V,R):
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
        
       
    
  def FINDAB(self,dfTOH):
    '''
    dfTOH is a dataframe with manoeuvre coefficents and the corresponding shallow water correction
    the routine find the a and b's to be used for shallow water correction, 
    the a's and b's that gives the minimum FUNCAB value is used
    it is assumed the double loop arround the minimize function is to try to ensure that a global minimum is found
    '''
    acoef={}
    bcoef={}
    for index,row in dfTOH.iterrows():
      # Initial guess
      X0 = [3, 1]
      ab=[3,1]
      XMIN = 100
      for I in range(3):
        ab[1]=2
        for J in range(4):
          minima = minimize(FUNCAB, X0,args=(dfTOH,index),
                    method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
          if minima.fun < XMIN:
            XMIN = minima.fun
            amin = minima.x[0]
            bmin = minima.x[1]
          X0[1] += 3.5
        X0[0]+=10
        acoef[index] = amin
        bcoef[index] = bmin
    return acoef,bcoef
    
  def HoltropResistance(self,velocitiesocities):
    '''
    calculate resistance according to Holtrop
    '''
    assert(True,'HoltropResistance not implementet fully yet')
    W5 = 0.0
    T5 = 0.0
    L5 = self.Lpp
    Frouden=np.zeros(len(velocitiesocities))
    CTotal = np.zeros(len(velocitiesocities))
    for ix in range(len(velocitiesocities)):
        V5 = velocities[I]
        if V5 < 0.0 :
            V5 = nb.abs(V5)     
        Frouden[I] = V5/np.sqrt(g *L5)
    if V5 != 0.0 :
        self.resis(SHAL,J,LER)
        pass
    else:
        R0 = 0.0
        
    if np.abs(V5) < 2.0E-3 :
        if  V5 <= 0.0 :
          C0 = 0.0
        else:
            C0 = np.abs(R0/(.5*1.025*S9*V5**2))
    else: 
        C0 = R0/(.5*1.025*S9*V5**2)
	     

#C-SHT	 *******	When Froude No. is less than .1 
##C				the resistance coefficient is a combination of 
#C				OCIMF and Holtrop:
#C				
#C				When abs(Fn)<=0.1:
#C				Cr = Cocimf*Focimf + Choltrop*(1-Focimf)
#C				Where 
#C				Focimf=abs((0.1-Fn)*(1/.1))*((0.1-Fn)*(1/.1))
#C				
#C				When abs(Fn)>=0.1:#
#C				Cr = Choltrop
#C 
    Cocimf=0.036*L5*T9/S9 #0.036 See OCIMF 1994
    Focimf=abs((0.1 - Frouden[I])*1/0.1)*((0.1 - Frouden[I])*1/0.1)

    if (Focimf > 0) :
      if (C0 != 0) :
        C0 = Focimf*Cocimf + (1-Focimf)*C0

    #C-SHT ******************************************			 
    if (velocities[I] < 0.0 ) :
      C0 = -1.1*C0
      V5 = -V5
      FROUDEN[I] = -FROUDEN[I]
    Ctotal[I] = -C0
    return FROUDEN,Ctotal



  
      
  def dmimix(self):
    '''
    translated from ShipYard1 msdat/lib/pmm/dmimix.f 
    by btj
    calcualte a and b shallow water corrections with the soccalled mix method
    that is a mix between different methods by Ankudinov,Kobayashi, Clarke
    '''
    
    TOH = [1./10., 1./5., 1./2., 1./1.5, 1./1.2]
    NTOH = len(TOH)
    #	Ankudinov
    CB = self.blockCoefficient
    BEAM = self.beam
    LENGTH = self.Lpp 
    DRAFT = self.meanDraft
    PI = np.pi
    B1OT = CB*BEAM*(1+BEAM/LENGTH)**2/DRAFT
    B2OT = B1OT*0.83/CB**2
    BEODR = BEAM/DRAFT
    SIGMA0 = 0.3
    #	Kobayashi
    BEODR = BEAM/DRAFT
    PKOBA = 1.0
    KVAR = 2.0*DRAFT/LENGTH
    NAEVNER = 0.5*PI*KVAR+PKOBA*CB*BEAM/LENGTH
    q1yv = 3.0		# $D$9
    q2yr = 1.2      # $E$11
    q2nv = 1.4      # $F$11
    q2nr = 0.5      # $G$11
    q3yvdot = 0.21  # $$H$12
    q3nrdot = 0.15  # $I$12
    q4yvdot = 1.2   # $H$13
    q4nrdot = 1.2   # $I$13
    #	Clarke
    a0 = 0.0775		# 0.0774
    b0 = -0.0110    # -0.0151
    a1 = -0.0643    # -0.0125
    b1 = 0.0742     # 0.1674
    c1 = -0.0113    # -0.0199
    a2 = 0.0242     # 0.0431
    
    index=['XUDOT','XVR','XVV','XRR',
            'YVDOT','YRDOT','YV','YVL','YR','YVIVI','YRIRI','YRRR','YVRR','YVIRI','YRIVI','YRVV',
            'NVIVI','NVVV','NRIRI','NRRR','NVIRI','NRVV','NRIVI'] 
    dfTOH=pd.DataFrame(np.zeros((len(index),len(TOH))),index=index,columns=TOH)
    
    for I in range(len(TOH)):
      #	Ankudinov
      FRATIO = 1.0/TOH[I]-1.0
      DEPTH = (1.0+FRATIO)*DRAFT
      K0I = 1+0.0775/FRATIO**2-0.011/FRATIO**3+0.000068/FRATIO**5
      K1I = -0.0643/FRATIO+0.0724/FRATIO**2- 0.0113/FRATIO**3+0.0000767/FRATIO**5
      if BEODR > 0.4 :
        K2I = 0.137*DRAFT/BEAM/FRATIO
      else:
        K2I = 0.0342/FRATIO
      
      SIGMAD = DRAFT/DEPTH-SIGMA0
      SIGMA = 1/(1+FRATIO)
      FNV = K0I+K1I*B1OT+K2I*B1OT**2
      FNR = K0I+1.0/2.0*K1I*B1OT+1/3*K2I*B1OT**2
      FYV = 1.5*FNV-0.5
      FYR = K0I+2.0/5.0*K1I*B1OT+24.0/105.0*K2I*B1OT**2
      GV = K0I+2.0/3.0*K1I*B1OT+8.0/15.0*K2I*B1OT**2
      GNR = K0I+8.0/15.0*K1I*B1OT+40.0/105.0*K2I*B1OT**2
      FRND = 1.0/12.0*FNR+11.0/12.0
      #	Kobayashi
      DUMMY = PI/2.0*DRAFT/DEPTH
      BIGPAREN = 1/(0.5*DRAFT/DEPTH*(KVAR+PI*1.0/np.tan(DUMMY)))
      LASTYV = PKOBA*CB*BEAM/LENGTH
      DOTPAR = np.tan(PI/2.0*DRAFT/DEPTH)
      #	Clarke
      K0 = 1+a0/FRATIO**2+b0/FRATIO**3
      K1 = a1/FRATIO+b1/FRATIO**2+c1/FRATIO**3
      K2 = a2/FRATIO
      #	In the following Ankudinov is used excepth where other 
      #	method is specified
      # Now we take all the Y derivatives:
      dfTOH.loc['YVDOT',dfTOH.columns[I]] = GV
      dfTOH.loc['YRDOT',dfTOH.columns[I]] = GV
      dfTOH.loc['YV',dfTOH.columns[I]] = FYV
      dfTOH.loc['YVL',dfTOH.columns[I]] = FYV
      #	Yr is Kobayshi
      dfTOH.loc['YR',dfTOH.columns[I]] = BIGPAREN**q2yr	
      dfTOH.loc['YVIVI',dfTOH.columns[I]] = 2.25*FNV-1.25
      dfTOH.loc['YRIRI',dfTOH.columns[I]] = FNR
      # NOTE that the same correction has been taken for YRIRI and YRRR
      dfTOH.loc['YRRR',dfTOH.columns[I]] = FNR
      # NOTE that the same correction has been taken for YVRR AND YVIRI
      dfTOH.loc['YVRR',dfTOH.columns[I]] = FYV
      dfTOH.loc['YVIRI',dfTOH.columns[I]] = FYV
      # NOTE that the same correction has been taken for YRVV AND YRIVI
      dfTOH.loc['YRIVI',dfTOH.columns[I]] = FYV
      dfTOH.loc['YRVV',dfTOH.columns[I]]  = FYV
      #
      # Now we take all the N derivatives:
      #
      dfTOH.loc['NVDOT',dfTOH.columns[I]] = GV
      #	Nrdot is Kobayashi
      dfTOH.loc['NRDOT',dfTOH.columns[I]] = 1+q3nrdot*DOTPAR**q4nrdot
      dfTOH.loc['NV',dfTOH.columns[I]] = FNV
      dfTOH.loc['NVL',dfTOH.columns[I]] = FNV
      #	Nr is CLarke   
      dfTOH.loc['NR',dfTOH.columns[I]] = K0+0.5*K1*BEODR+1.0/3.0*K2*BEODR*2

      # NOTE that the same correction has been taken for NVIV and NVVV
      dfTOH.loc['NVIVI',dfTOH.columns[I]] = 2.25*FNV-1.25
      dfTOH.loc['NVVV',dfTOH.columns[I]] = 2.25*FNV-1.25
      # NOTE that the DMI default has been taken for NRIRI and NRRR
      #        FNCTOH(I,iof + INRIRI) = GV
      #	  FNCTOH(I,iof + INRRR) = GV
      #NOTE that the same correction has been taken for NVRR and NVIRI
      dfTOH.loc['NVRR',dfTOH.columns[I]] = GNR
      dfTOH.loc['NVIRI',dfTOH.columns[I]] = GNR
      # NOTE that the same correction has been taken for NRVV and NRIVI
      dfTOH.loc['NRVV',dfTOH.columns[I]] = GNR
      dfTOH.loc['NRIVI',dfTOH.columns[I]] = GNR
      #
      # Now we take all the X derivatives:
      #
      dfTOH.loc['XUDOT',dfTOH.columns[I]] = GV
      dfTOH.loc['XVR',dfTOH.columns[I]] = FNV
      dfTOH.loc['XVV',dfTOH.columns[I]] = FNV
      dfTOH.loc['XRR',dfTOH.columns[I]] = GNR

    acoef,bcoef = self.FINDAB(dfTOH)
    return acoef,bcoef
      
  def hluref(self,udim,vdim,rdim,pdim):
    '''
    calculated speeds used for dimensionalising
    '''
    
    lpp = self.Lpp
    b = self.beam
    dm = self.meanDraft
    
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
      
  def pmmcar(self,uo,betad,coftyp):
    '''
    BTJ assume pmmcar is short for pmm carriage
    uo      : speed of carriage
    betad   : drift angle
    coftyp  : coefficient type 1 -> pmm  else msm
    
    '''
    assert(betad >3.15, "pmmcar called with degrees ?")
    if coftyp == 1 :
      ucar = uo/np.cos(betad)
    else:
      ucar = uo
    return ucar

  def defval(self,motion,icoty):
    '''
    translated (modified and cleaned) to python by BTJ originaly inside tfs
      $/SimFlex Classic/src/lib/core/msdat/lib/pmm/pmms.f
    the function returns the angle series to be used in the tables generated
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

    cbla*******************************************************************************
    cbla comment number001
    cbla*******************************************************************************
    cbla   The changes in points below are made to make use of the right 4 quadrant 
    cbla   special coefficients in the interval outside of pmm range. The idea is
    cbla   to use pmm coefs in a range around 0 degrees both in betad and gamma, 
    cbla   outside of the range where we use predicted pmm coefs, we use some special
    cbla   four quadrant coefs, these are however only valid around betad=90 degrees,
    cbla   and gamma=90 degrees, to prevent using these coefs in the whole range where 
    cbla   these special coefs are used, some points have been moved into the pmm range.
    cbla   For betad in between 180.0 and 185.0-SepPoint, some old MSM coefs are used.
    cbla
    cbla   If changes are made in the points below care should be taken, to avoid using 
    cbla   the special fixed four quadrants coefficients, for betads and gammas 
    cbla   different than a narrow band around 90 degrees.
    cbla*******************************************************************************        

    '''
    absc1 = None
    absc2 = None
    ndef1 = None
    ndef2 = None
    if motion == 1 and icoty == 0: 
      absc1=[-180.0,-177.0,-175.0,-170.0,-95.0,-90.0,
              -85.0,-22.0,-20.0,-18.0,-15.0,-10.0,-5.0,-2.0,-1.0,
                0.0,1.0,2.0,5.0,10.0,15.0,18.0,20.0,22.0,85.0,90.0,
                95.0,170.0,175.0,177.0,180.0]
      ndef1 = len(absc1)
      idof  =  0    
#       betad and toh
    if motion ==  1 and icoty == 1 : 
      absc1= [-180.,-170., -90., -20., -15.,  -5.,   5.,  15.,  20.,  90., 170., 180.]
      ndef1 = len(absc1)
      idof  =  2
#cbla   see also 'cbla comment number001'.
    if motion == 2 and icoty == 0: 
      absc1=[-90.,-26.,-24.,-22.,-20.,-15.,-10., -5., -2., -1., 
              0.0,  1.,  2.,  5., 10., 15., 20., 22., 24., 26., 90.]
      ndef1 =  len(absc1)
      idof  =  0
    
    if motion == 2 and icoty == 1:
      absc1=[-90.,  -20.,  -15.,   -5.,    5.,   15.,   20.,   90.]
      ndef1 =  len(absc1)
      idof  =  2
#cbla   see also 'cbla comment number001'.
    if motion == 3 and icoty == 0:
      absc1 =[-180.0,-165.0, -90.0, -25.0, -22.0,-20.0, -10.0 -5.0, -2.0,
                0.0, 2.0,5.0,10.0,20.0,22.0,25.0,90.0,165.0,180.0]
      ndef1 = len(absc1)
      absc2= [-90.0, -25.0, -20.0, -15.0, -10.0,  -5.0,  -2.0,   
                0.0,   2.0,   5.0,  10.0,  15.0,  20.0,  25.0,  90.0]
      ndef2 = len(absc2)
      idof  =  0
      
    if motion == 3 and icoty == 1 : 
      idof  =  1

    if motion == 4 :
      absc1 = [-45.,-20. ,-10. , -5. ,  0.,  5., 10., 20., 45.]
      ndef1 =  len(absc1)
      if  icoty == 0 :
        idof  =  0
      else:
        idof  =  2
      
    if motion == 5: 
      absc1=[-180.0, -160.0, -135.0,  -90.0,  -45.0,  -20.0,  -10.0,
                0.0,   10.0,   20.0,   45.0,   90.0,  135.0,  160.0,  180.0]
      ndef1 =  len(absc1)
      absc2 =[-45., -20. , -10. ,  -5. ,   0.,   5.,  10.,  20.,  45.]
      ndef2 =  len(absc2)
      if icoty == 0:
        idof = 0
      else :
        idof = 1

    if motion == 6 :
      absc1 = [-90.0,-60.0,-45.0,-25.0,-10.0,  0.0,  10.0,  25.0,  45.0,  60.0 , 90.0]
      ndef1 =  len(absc1)
      absc2=[-45.,-20. ,-10. , -5. ,  0.,  5., 10., 20., 45.]
      ndef2 =len(absc2)
      if icoty == 0: 
        idof  =  0
      else:
        idof  =  1

    if motion == 8 and icoty == 0: 
      absc1 =[ -180.0,-177.0, -90.0, -45.0, -38.0, -30.0, -25.0, -20.0, -18.0, -15.0,
              -11.0,  -8.0,  -5.0,  -2.0,  -1.0,   0.0,   1.0,   2.0,   5.0,   8.0,
                11.0,  15.0,  18.0,  20.0,  25.0,  30.0,  38.0,  45.0,  90.0, 177.0, 180.0]
      ndef1 =  len(absc1)
      idof  =  0
      
    if idof == 1 : 
      absc1 = [0.0, 0.3, 0.5, 0.7, 1.0]
      ndef1 =  len(absc1)
    if idof == 2 : 
      absc2 = [0.0, 0.3, 0.5, 0.7, 1.0]
      ndef2 = len(absc2)
    return absc1,absc2
      
    

      
  def pmmmsm(self,uo,betad,gamma,pmmtyp,SepPoint):
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
    if np.abs(np.abs(gamma)-np.radians(90.) < .01):
        ucar = 0.
    return coftyp,ucar   
      
      
      
      
  def pmmfor(self,pmmcoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh):
    '''
    this function calculates the forces on the hull using the PMM coefficient/hydrodynamic derivatives
    the function is translate to python by BTJ originally part of shipYard1 see tfs
    $/SimFlex Classic/src/lib/core/msdat/lib/pmm/pmmfor.f
    '''
    acoef = self.acoef
    bcoef = self.bcoef
    fdim=np.array([0.,0.,0.])
    fdimu2=np.array([0.,0.,0.])
    ## SKIPPING D,Q,P derivatives as we dont have them
    xderivatives=['X0','XU','XUU','XUUU','XV','XVV','XR','XRR','XVR','XUDOT','XRDOT','XUIUI']
    yderivatives=['Y0','Y0U','YUU','YV','YVV','YVIVI','YVVV','YVU','YR','YRR','YRIRI','YRRR','YRIVI','YRVV','YVIRI','YVRR','YVDOT','YRDOT']
    nderivatives=['N0','N0U','NUU','NV','NVV','NVIVI','NVVV','NVU','NR','NRR','NRIRI','NRRR','NRIVI','NRVV','NVIRI','NVRR','NVDOT','NRDOT ']
    lpp = self.Lpp
    



        
    if coftyp == 1:
      utot2 = udim*udim + vdim*vdim
      utot  = np.sqrt(utot2)
      U   = (udim - uo)/utot
      V   = vdim/utot
      R   = rdim/(utot/lpp)
      Q  = qdim
      P   = pdim/(utot/lpp)
      D   = ddim
      USIGN = 1.0
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
      usign = np.sign(1.0,U)
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
          speedvector=np.append(speedvector,self.speedfactor(val,U,V,R))
          if val in acoef.keys():
            a = acoef[val]
            b = bcoef[val]
          else:
            a = 0
            b = 1
          TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))
          
      fdimu2[0]= np.sum(xcoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**2
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
          speedvector=np.append(speedvector,self.speedfactor(val,U,V,R)) 
          TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))       
      
      fdimu2[1]= np.sum(ycoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**2
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
          speedvector=np.append(speedvector,self.speedfactor(val,U,V,R)) 
          TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))        
      
      fdimu2[2]= np.sum(ncoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**3
      fdim[2] = fdimu2[2]*utot2      
      
      return fdim,fdimu2,utot2      
        
        
# %%
import matplotlib.pyplot as plt
if __name__ == '__main__':
  
  shipDatadict={}
  shipDatadict['shipnr'] = 3949
  shipDatadict['lpp'] = 312
  shipDatadict['Beam'] = 53.5
  shipDatadict['wettedSurface'] = 23243.8
  shipDatadict['waterPlaneArea'] = 0.893628117 * shipDatadict['lpp'] * shipDatadict['Beam']
  shipDatadict['propellerType'] ='FP'
  shipDatadict['PropellerPitch'] = 0.71505
  shipDatadict['displacement'] =218220.0
  shipDatadict['propellerDiameter'] = 7.0
  shipDatadict['draftAft']  = 17.
  shipDatadict['draftFore'] = 17.0
  shipDatadict['meanDraft'] = (shipDatadict['draftAft'] + shipDatadict['draftFore'] )/2.0
  shipDatadict['underWaterLateralArea'] = 0.9124 * shipDatadict['lpp'] * shipDatadict['meanDraft']
  shipDatadict['blockCoefficient'] = 0.704661757
  shipDatadict['CenterofGravity'] = np.array([-0.002290749, 0, -0.415058824]) * np.array([shipDatadict['lpp'] ,shipDatadict['Beam'],shipDatadict['meanDraft']])
  shipDatadict['verticalCenterOfBoyancy'] = 0.453647058824 * shipDatadict['meanDraft'] 
  shipDatadict['GyrationArms'] = np.array([0.4, 0.25, 0.25]) * np.array([shipDatadict['Beam'],shipDatadict['lpp'],shipDatadict['lpp']])

  hull_Thumbs = HullThumbs(shipDatadict)
  # !!! dmimix takes a bit of time !!!
  #a,b =hull_Thumbs.dmimix()  
  with open('data\PMMdata\PMMCoefs3686.dat','r') as f:
    pmmCoefs = {}
    for line in f:
      splt = line.split()
      pmmCoefs[splt[0]] = float(splt[1])/1.0E5
  baseTable,multiPliers = hull_Thumbs.getBaseTable(pmmCoefs,'DRIFT','X_HL')
  plt.plot(baseTable['BETAD'],baseTable['X_HL'],label='X_HL(BETAD)')
  plt.legend()
  plt.grid()
  plt.show()
  pass

           
          

        
        
        
        
        
      
    



# %%

# %%
