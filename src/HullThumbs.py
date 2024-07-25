

# %%
import numpy as np
from scipy.constants import g
from scipy.optimize import fmin_bfgs,minimize
from scipy.interpolate import interp1d
import json
import os
import glob

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
    self.MultiPliers={'DRIFT':{'X_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea},
                           'Y_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea},
                           'K_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea,'DM':self.meanDraft},
                           'N_HL':{'WHRD2':self.rho/2.,'UR_D':1,'ALW':self.underWaterLateralArea,'LPP':self.Lpp}
                          },
                      'YAW':{'X_HL':{'WHRD2':self.rho/2.,'UR_Y':1,'ALW':self.underWaterLateralArea},
                           'Y_HL':{'WHRD2':self.rho/2.,'UR_Y':1,'ALW':self.underWaterLateralArea},
                           'N_HL':{'WHRD2':self.rho/2.,'UR_Y':1,'ALW':self.underWaterLateralArea,'LPP':self.Lpp}
                          },
                      'YAW_DRIFT':{'X_HL':{'WHRD2':self.rho/2.,'UR_YD':1,'ALW':self.underWaterLateralArea},
                           'Y_HL':{'WHRD2':self.rho/2.,'UR_YD':1,'ALW':self.underWaterLateralArea},
                           'N_HL':{'WHRD2':self.rho/2.,'UR_YD':1,'ALW':self.underWaterLateralArea,'LPP':self.Lpp}
                          }
                      }
    self.forceIndex={'X_HL':0,'Y_HL':1,'N_HL':2}
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
  
  def RESIS(self,shallow,velocity):
    '''
    This function calculates the resistance according to Holtrop Mennen , converted to python from 
    $/SimFlex Classic/src/lib/core/msdat/lib/thmb/resis.f
    shallow: True -> shallow water resistance calculated else deep
    velocity: speed through water in m/s
    
    note: indirectly tested through test_resistance
    '''
    V5 = velocity
    M9 = 1.026 ## as in fortran code
    S4 = M9 - 1
    T6 = 15.0   #! Water temperature
    S9 = self.wettedSurface
    K8 = self.WettedSurfaceAppendage
    K8 = 1 + K8/S9	# Wetted surface correction(total ship) (K8 was
	                  # before this statement the wetted surface for appendages!)
    L9 = self.lengthWaterLine
    B9 = self.beam
    T9 = self.meanDraft
    V0 = self.displacement
    C3 = self.midshipSection
    FORMAL = self.Formal
    O1 = self.LCB_ratio
    C8 = self.waterLineBlock
    C9 = self.prismaticCoefficeint    
    CWAT = self.waterPlaneArea/( self.Lpp * self.beam)
    ATRANS = self.ship_data['transomArea']
    M = self.someCoefficient
    GRAV = g
    K6R = self.verticalCenterBulb
    FMFORE = self.draftFore
    ABULB = self.Abulb
    XBULB=ABULB*B9*T9*C3


    Cx3=self.Cx3()    
    LRUN = L9*(1-C9+0.06*C9*(O1*100)/(4*C9-1)) #! Length of run from Holtrop paper
                                               #C -   Remember that lcb in Holtrop paper is percentage!
    K6A= 1+89*np.exp((-(L9/B9)**.80856)*(1-CWAT)**.30484   \
     	*(1-C9-2.25*O1)**.6367*(LRUN/B9)**.34574			\
     	*(100*V0/L9**3)**.16302) # ! iE from Holtrop,half angle of incid.
    VISCO = (12.0433-31.38*S4)*(T6+20)**(-.482+1.72*S4)-1.0312 -5.779*S4
    VT = V5
    IT = 0
    CBAR = 0
    F1 = V5/np.sqrt(L9*9.80665)
    K5 = .5*np.log10(V0)-.1*np.log10(V0)**2
    R1 = V5*L9/VISCO*1000000  
    if  shallow:
      R1 = VT*L9/VISCO*1000000  
    if  R1 < 1000.:
      R1=1000.
    C6 = 1000*.075/(np.log10(R1)-2)**2
    RFRIC = .5*M9*VT**2*C6*S9/1000
    if shallow :
      RFRIC = .5*M9*VT**2*C6*FORMAL*S9*K8/1000
    if not shallow:
      R0 = .5*M9*V5**2*C6*FORMAL*S9*K8/1000 #!This is: Rf(1+k1) including appendages    
    CX7 = .229577*(B9/L9)**.33333
    if B9/L9 > .11:
      CX7 = B9/L9
    if B9/L9 > .25:
      CX7 = .5-.0625*L9/B9
    CX1 = 2223105*CX7**3.78613*(T9/B9)**1.07961*(90-K6A)**(-1.37565)
    CX3 = self.Cx3()
    CX2 = np.exp((-1.89)*np.sqrt(CX3))
    CX5 = 1-.8*ATRANS
    CX16 = 8.07981*C9-13.8673*C9**2+6.984388*C9**3
    if C9 > .8: 
      CX16 =1.73014-.7067*C9
    M1 = .0140407*L9/T9-1.75254/M-4.79323*B9/L9-CX16
    CX15 = -1.69385
    if (M > 8):
       CX15 = -1.69385+(M-8)/2.36
    if (M > 12) :
      CX15 = 0 
    LDA = 1.446*C9-.03*L9/B9
    if (L9/B9 > 12) :
      LDA = 1.446*C9-.36
    FNX = F1
    if (FNX > .4):
     FNX = .4
    RXP = 0
    if (np.abs(FNX) > 1.E-6) :
      rdummy=(-.034)*FNX**(-3.29)
      if (rdummy < -80.0) :
        rdummy=0.0
      else:
	      rdummy=np.exp(rdummy)
	      rdummy=CX15*.4*rdummy
      M2=rdummy
      rdummy=(M1*FNX**(-.9)+M2*np.cos(LDA/FNX**2))
      if (rdummy < -80.0) :
        rdummy=0.0
      else:
	      rdummy=np.exp(rdummy)
	      rdummy=CX1*CX2*CX5*rdummy*V0*M9*GRAV
      RXP=rdummy
    if (F1 > .55):
       RXP = 0
    if (F1 > .4):
       RXP = RXP*(5.5-10*F1)/1.5
    FNX = F1
    if (FNX < .55):
       FNX = .55
    M2 = CX15*.4*np.exp((-.034)*FNX**(-3.29))
    CX17 = 6919.3*C3**(-1.3346)*(V0/L9**3)**2.00977*(L9/B9-2)**1.40692
    M3 = (-7.2035)*(B9/L9)**.326869*(T9/B9)**.605375
    RYP = CX17*CX2*CX5*np.exp(M3*FNX**(-.9)+M2*np.cos(LDA/FNX**2))*V0*M9*GRAV # Rw-b for FN>0.55
    if (F1 < .4):
       RYP = 0
    if (F1 < .55):
       RYP = RYP*(10*F1-4)/1.5 #Interpolation 0.4<FN<.55  
    if (shallow) :
        R0 = RXP+RYP
    else:
        R0 = R0+RXP+RYP
    #CJBP 990222 Holtrop 1984 recommends this statement, i.e. not in I-SHIP
    if ( K6R > 0.6*FMFORE ) :
      K6R = 0.6*FMFORE
    #CJBP Holtrop 1998 bulb correction as 1984 seems too strong for the Lindø VLCC!
    if ( ABULB > 0.12 ) :
	    XBULB = 0.12*B9*T9*C3    #! To ensure reasonable bulb resistance
    FORSINHF = B9*T9*C9*C3*(136.0-316.3*F1)*F1**3/L9
    if ( FORSINHF < -0.01*L9 ):
       FORSINHF = -0.01*L9
    LOCWAVHW = V5**2*K6A/(400*GRAV)
    if ( LOCWAVHW > 0.01*L9 ) :
      LOCWAVHW = 0.01*L9
    FNI = V5/np.sqrt(GRAV*(FMFORE-K6R-0.25*np.sqrt(XBULB)+FORSINHF+LOCWAVHW))
    PB = 0.56*np.sqrt(XBULB)/(FMFORE-1.5*K6R+FORSINHF)
    RBULB = .11*np.exp((-3)/PB**2)*FNI**3*XBULB**1.5*M9*GRAV/(1+FNI**2)    
    if ( ABULB > 0.12 ) :
	    XBULB = ABULB*B9*T9*C3    # To reset bulb area for wave resistance
    #C --- Here we do the transom stern correction
    if (np.abs(ATRANS) > 1.E-6) :
      if (shallow) :
        FNT = VT/np.sqrt(2*GRAV*ATRANS*T9*C3/(1+CWAT))
      else: 
        FNT = V5/np.sqrt(2*GRAV*ATRANS*T9*C3/(1+CWAT))
      CX6 = 0
      if (FNT < 5):
        CX6 = .2*(1-.2*FNT)
      if (shallow) :
        R0 = R0+.5*M9*VT*VT*ATRANS*CX6*B9*T9*C3
      else:
        R0 = R0+.5*M9*V5*V5*ATRANS*CX6*B9*T9*C3 #R0 + Rtr in Holtrop
    CX4 = FMFORE/L9
    if (CX4 >.04) :
      CX4 = .04
    CXA = .00546*(L9+100)**(-.16)-.00205+.003*np.sqrt(L9/7.5)*C8**4*CX2*(.04-CX4) #This is model ship correlation coefficient Ca
    if CXA < 0.0 :
      CXA = 0.0
    if (shallow) :
      R0 = R0+.5*M9*VT**2*CXA*S9*K8+RFRIC
    else:
      R0 = R0+.5*M9*V5**2*CXA*S9*K8         
    return R0
    
  def RDHOLTROP(self,shallow):  
    '''
    This routine calculate the resistance table according to Holtrop
    '''
    MAXSPEED =  1.4 * self.serviceSpeed
    MINSPEED = -1.0 * self.serviceSpeed
    NAHEAINT = 14
    NASTINT = 8
    velnegative = np.linspace(MINSPEED,-1.0E-2,NASTINT,endpoint=True)
    velpositive = np.linspace(0,MAXSPEED,NAHEAINT,endpoint=True)
    #set the first positive speed as in fortran code
    velpositive[1] = 1.0E-2
    velocities = np.concatenate((velnegative, velpositive))
    Frouden = np.zeros(len(velocities))
    Ctotal = np.zeros(len(velocities))
    S9 = self.wettedSurface
    L5 = self.Lpp
    T9 = self.meanDraft
    
    for i in range(len(velocities)):
      v5 = velocities[i]
      if velocities[i] < 0.0:
          v5 = abs(v5)
      Frouden[i] = v5/np.sqrt(9.81*self.Lpp)
      if v5 != 0.0:
          R0 = self.RESIS(shallow,v5)
      else:
          R0 = 0.0
      if np.abs(v5) < 2.0E-3:
          if v5 == 0.0:
              C0 = 0.0
          else:
              C0 = np.abs(R0/(.5*1.025*S9*v5**2))
      else:
          C0 = R0/(.5*1.025*S9*v5**2)

      # When Froude No. is less than .1 the resistance coefficient is a combination of OCIMF and Holtrop
      Cocimf = 0.036*L5*T9/S9  # 0.036 See OCIMF 1994
      Focimf = abs((0.1 - Frouden[i])*1/0.1)*((0.1 - Frouden[i])*1/0.1)
      if Focimf > 0:
          if C0 != 0:
              C0 = Focimf*Cocimf + (1-Focimf)*C0
      if velocities[i] < 0.0:
          C0 = -1.1*C0
          v5 = -v5
          Frouden[i] = -Frouden[i]
      Ctotal[i] = -C0
      #print(v5, ' ', Frouden[i], ' ', Ctotal[i], ' ', R0)
      
    df = pd.DataFrame(index=Frouden,data=Ctotal)
    df.index.name = 'FN'
    df = df.rename(columns={0:'X_HL'})
    return df

  def shallowWaterResistance(self):
    ''' 
    translated from $/SimFlex Classic/src/lib/core/msdat/lib/thmb/hlshrs.f
    the originalname may be translated to hull shallow resistanc schlicting
    The code is kept close to the origanl implementation even though it seems there is a potential for optimizing the code
    but i works so have used time to butify the code BTJ/20240725
    
    '''
    XLPP = self.Lpp
    XBEAM = self.beam
    XDRAFP = self.draftFore
    XDRAAP = self.draftAft
    XWETSU = self.wettedSurface
    XMIDCO = self.midshipSection
    V5 = self.serviceSpeed
    XWAVIS = 1.191e-6
    xwaden = self.rho
    DRAFTM = self.meanDraft
    AX = XBEAM * DRAFTM * XMIDCO    # Mid ship area
    NAHEAINT = 14
    NASTINT = 8
    TOHARR = np.array([0.000, 0.400,0.667,0.833,1.000])
    ITOH = len(TOHARR)
    MAXSPEED = 1.4*V5
    MINSPEED = -1.0*V5
    velnegative = np.linspace(MINSPEED,-1.0E-2,NASTINT,endpoint=True)
    #dont include zero in the velocities
    velpositive = np.linspace(1.0E-2,MAXSPEED,NAHEAINT-1,endpoint=True)
    velocities = np.concatenate((velnegative, velpositive)) 
    MFRU = len(velocities)
    FRUDEP = np.array([ vel/np.sqrt(g*XLPP) for vel in velocities])
    FRUNOH = np.zeros(len(FRUDEP))
    SHAFRU = np.zeros(len(FRUDEP))
    VFINAL = np.zeros(len(FRUDEP))
    CTSHAL = np.zeros(len(FRUDEP))
    CTSCTD = np.zeros(len(FRUDEP))
    CT = np.zeros(ITOH)
    CRITSP = np.zeros(ITOH)
    FRU = np.zeros(MFRU)
    CTS = np.zeros(MFRU)
    CTSH1  =  np.zeros(MFRU)
    CTS1   =  np.zeros(MFRU)
    FRU1    =  np.zeros(MFRU)
    CTSH2  =  np.zeros(MFRU)
    CTS2   =  np.zeros(MFRU)
    FRU2    =  np.zeros(MFRU)
    CTSH3  =  np.zeros(MFRU)
    CTS3   =  np.zeros(MFRU)
    FRU3    =  np.zeros(MFRU)
    CTSH4  =  np.zeros(MFRU)
    CTS4   =  np.zeros(MFRU)
    FRU4    =  np.zeros(MFRU)
    CTSH5  =  np.zeros(MFRU)
    CTS5   =  np.zeros(MFRU)
    FRU5    =  np.zeros(MFRU)              
    IFNO = len(FRUDEP)  
    ICRIT = 0
    FTEMP = np.zeros((IFNO,ITOH))

    SQAXOH =[ 0.0, 0.1, 0.2, 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 ,0.9 , 1.0 , 1.1 ,   1.2 ,   1.3 ,
              1.4 ,1.5 ,1.6 ,1.7 , 1.8 , 1.9 , 4.0 ]
    VHVI = [1.0, 1.0, 0.999, 0.995, 0.987, 0.977, 0.962, 0.946, 0.926, 0.905, 0.883, 0.858, 
        0.830, 0.801, 0.7925, 0.743, 0.710, 0.6795, 0.644, 0.608, 0.155]
        
    RINTP_SQAXOH = interp1d(SQAXOH,VHVI)    
    shallow = False
    CTOTD = self.RDHOLTROP(shallow)
    CTOTD = CTOTD.iloc[[i for i in range(len(CTOTD)) if i != NASTINT]]
    
    RINTP_FRUDEP = interp1d(FRUDEP,CTOTD['X_HL'].values)


#
#---- CALCULATION OF SHALLOW WATER RESISTANCE AND CT-H/CT-INF ----
#
    for L in range(len(CTOTD)):
      for J in range(len(TOHARR)):
        if J == 1:
          idum = 0
        VSPEED = velocities[L]
        if VSPEED  < 0.0 :
          REYNOM = -VSPEED * XLPP/XWAVIS
        else:
          REYNOM = VSPEED * XLPP/XWAVIS
        
        CFDEEP = -np.sign(FRUDEP[L])*0.075/(np.log10(REYNOM) - 2.0)**2
        DELTAC = CTOTD.iloc[L]['X_HL'] - CFDEEP
        if TOHARR[J] == 0.0:
           TOHARR[J] = 1.E-05
        DEPTH  = DRAFTM / TOHARR[J]

        if (L == 0) : 
          CRITSP[J] = 0.75*np.sqrt(9.81*DEPTH)*3600.0/1852.0 

        SQAXH = np.sqrt(AX)/DEPTH
        #C
        #C--- Vi/Vu IS CALCULATED  FOR ACTUAL SPEED ---
        #C
        XVIVU = np.sqrt(np.tanh(9.81*DEPTH/(VSPEED**2)))    
        #
        #--- CFi AND CTi ARE CALCULATED AT SPEED Vi ---
        #
        FRUNOI = FRUDEP[L] * XVIVU
        VSPED = FRUNOI * np.sqrt(9.81 * XLPP)
        if (FRUDEP[L] < 0.0) :
          REYNO1 = -VSPED * XLPP/XWAVIS
        else:
          REYNO1 =  VSPED * XLPP/XWAVIS
        
        CFI = -np.sign(FRUDEP[L])*0.075/(np.log10(REYNO1) - 2.0)**2
        CTI = CFI + DELTAC*(VSPEED**2/VSPED**2)  
        #C
        #C---- LOOKUP IN VH/VI TABLE FOR SQRT(AX)/H ---
        #C
        FRUSHA = RINTP_SQAXOH(SQAXH)

        FRUNOH[L] = FRUNOI * FRUSHA
        VSPED1 = FRUNOH[L] * np.sqrt(9.81 * XLPP)
        V1 = VSPED1 * (3600.0/1852.0)
        #C
        #C--- CT SHALLOW IS CALCULATED ON BASE OF CTi AND Vi/Vh ---
        #C
        CTSHAL[J] = CTI * (VSPED**2/VSPED1**2)
        #C
        #C---- LOOKUP FOR CTOTD IN CTOTD TABLE FOR FN IN SHALLOW WATER
        #C
        CTFRUH = RINTP_FRUDEP(FRUNOH[L])

        #C
        #C--- CT-H/CT-INF ARE CALCULATED ---
        #C
        CTSCTD[J] = CTSHAL[J]/CTFRUH
        if (J == 0) :
          CTSH1[L]  = CTSHAL[J]*1000.0
          CTS1[L]   = 1.0000 #CTSCTD[J]
          FRU1[L]   = FRUNOH[L]

        if (J == 1) :
          CTSH2[L]  = CTSHAL[J]*1000.0
          CTS2[L]   = CTSCTD[J]
          FRU2[L]   = FRUNOH[L]
      
        if (J == 2) :
          CTSH3[L]  = CTSHAL[J]*1000.0
          CTS3[L]   = CTSCTD[J]
          FRU3[L]   = FRUNOH[L]

        if (J == 3) :
          CTSH4[L]  = CTSHAL[J]*1000.0
          CTS4[L]   = CTSCTD[J]
          FRU4[L]   = FRUNOH[L]

        if (J == 4) :
          CTSH5[L]  = CTSHAL[J]*1000.0
          CTS5[L]   = CTSCTD[J]
          FRU5[L]   = FRUNOH[L]
       
        #C
        #C--- FN AND SPEED FOR SHALLOW ARE SAVED FOR FIRST TOH ---
        #C
        if (J == 0) :
          SHAFRU[L] = FRUNOH[L]
          VFINAL[L] = V1
       
        #C
        #C--- FN AND CTS ARE TRANSFERED TO AN WORKING ARRAY TO BE USED IN FINDING
        #C    CTS FOR THE ACTUAL FN.
        #C
    for K1 in range(IFNO):
      for K2 in range(ITOH):
        if (K2 == 0) : 
          for I in range(MFRU):
            FRU[I] = FRU1[I]
            CTS[I] = CTS1[I]

        if (K2 == 1) : 
          for I in range(MFRU): 
            FRU[I] = FRU2[I]
            CTS[I] = CTS2[I]

        if (K2 == 2) : 
          for I in range(MFRU):  
            FRU[I] = FRU3[I]
            CTS[I] = CTS3[I]

        if (K2 == 3) : 
          for I in range(MFRU):
            FRU[I] = FRU4[I]
            CTS[I] = CTS4[I]
        if (K2 == 4) : 
          for I in range(MFRU):  
            FRU[I] = FRU5[I]
            CTS[I] = CTS5[I]
     
        SHAFU_RINTP = interp1d(FRU,CTS)
        if SHAFRU[K1] < FRU[0]:
          CT[K2] = CTS[0]
        elif SHAFRU[K1] > FRU[IFNO-1]:
          print('NOT sure how to do this correct !!!')
          CT[K2] =SHAFU_RINTP(FRU[IFNO-1])
        else:
          CT[K2] =SHAFU_RINTP(SHAFRU[K1])

        #C
        #C--- CHECK FOR SPEED NOT EXCEEDS .75 TIMES CRITICAL SPEED --
        #C
        VS1 = np.abs(SHAFRU[K1]) * np.sqrt(9.81 * XLPP)
        if ((VS1*3600.0/1852.0) > CRITSP[K2]) :
          ICRIT = ICRIT + 1
        #C jbp The values in the critical range are now "adjusted", before they were zero 
          if (CT[K2]  > 5.0 ) :
            CT[K2] = 5.0 + 0.1*(CT[K2]-5.0)

        #C
        #C--- CTH/CTINF ARE WRITTEN TO A TABEL FOR SAME FN IN A ROW ---
        #C  
      #END K2 LOOP fil the row in FTEMP
      for ix in range(ITOH):
        FTEMP[K1,ix]= CT[ix]
      
      

      #if ( ICRIT !=  0 ) :
      #  print(f'<Critical speed exceeded  {ICRIT} times>')
    df = pd.DataFrame(FTEMP,columns=TOHARR)
    df.index = FRUDEP
    return df

    
  def pmmmot(self,ucar,betad,gamma,delta,heel,epsil):
    '''
    BTJ: think the pmmmot is short for pmm motion
        ucar could be speed of carriage in the towingtank
    the routine returns 6 dof dimensional speeds depending on speed and driftangle,turnrate..
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

  def updateMultiplierWithActualSpeed(self,multiPlierDict,speed_dict):
    '''
    if a key from multiplierDict is found in speed_dic then 
    update the actual multiplier dict with the actual speed
    '''
    for key in multiPlierDict.keys():  
      if key in speed_dict.keys():
        multiPlierDict[key] = speed_dict[key]
    return multiPlierDict
    
  def getForceCoefficient(self,multiPlier,force,utot2):
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
    
  def yaw_drift_Tables(self, motion, tableType, absc1, absc2, icoty, uo, pmmcoefs):
    multipliers = self.MultiPliers[tableType]
    force_coefficients = {key: pd.DataFrame(np.zeros((len(absc1), len(absc2))),index=absc1,columns=absc2) for key in multipliers.keys()} 
    correction_table = {key: np.zeros((5,2)) for key in multipliers.keys()}

    for ixa1, betaDEG in enumerate(absc1):
      betad = np.radians(betaDEG)
      for ixa2,gammaDEG in enumerate(absc2):
        gamma = np.radians(gammaDEG)
        
        pmmtyp = 1
        SepPoint = self.SeparationPoint
        coftyp, ucar = self.pmmmsm(uo, betad, gamma, pmmtyp, SepPoint)

        uoo = ucar
        ucar = self.pmmcar(uoo, betad, coftyp)
        tertiary = delta = heel = epsil = 0.0  
        #calulate pure drift force F1 setting gamma 0.0
        udim, vdim, rdim, qdim, pdim, ddim = self.pmmmot(ucar, betad, 0.0, delta, heel, epsil)

        toh = 0
        uo = self.serviceSpeed
        fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)
        F1 = fdimu2
        #calculate pure yaw F2 setting beta = 0
        ucar = self.pmmcar(uoo, betad, coftyp)
        udim, vdim, rdim, qdim, pdim, ddim = self.pmmmot(ucar, 0.0, gamma, delta, heel, epsil)
        fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)
        F2 = fdimu2
        
        #calculate combined
        ucar = self.pmmcar(uoo, betad, coftyp)
        udim, vdim, rdim, qdim, pdim, ddim = self.pmmmot(ucar, betad, gamma, delta, heel, epsil)
        fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)
        FCombined = fdimu2          
        speed_dict = self.hluref(udim, vdim, rdim, pdim)
        # update the multipliers with speed values if they are available
        multipliers = {key: self.updateMultiplierWithActualSpeed(value, speed_dict) for key, value in multipliers.items()}          
        F_Yaw_Drift = FCombined - F1 - F2
        #print(f"BETA,GAMMA {betaDEG},{gammaDEG} FCombined {FCombined[0]} F1 {F1[0]} F2 {F2[0]}")
        
        for key, multiplier in multipliers.items():            
          force_coefficients[key].iloc[ixa1,ixa2] = self.getForceCoefficient(multiplier, F_Yaw_Drift[self.forceIndex[key]], utot2)  
          idum = 0
          if key == 'X_HL' and betaDEG == 10.0 and gammaDEG == -25.0 :
            print(f"X_HL yaw_drift {betaDEG},{gammaDEG} = {force_coefficients[key].iloc[ixa1,ixa2]}")
    # pack to dataFrames      
          
          
    return force_coefficients, correction_table    
    
  def calculate_force_coefficients(self, motion, tableType, absc1, absc2, icoty, uo, pmmcoefs):
      multipliers = self.MultiPliers[tableType]
      force_coefficients = {key: {} for key in multipliers.keys()}
      correction_table = {key: np.zeros((len(absc1), len(absc2))) if absc2 else None for key in multipliers.keys()}
      forceIndex=self.forceIndex
      for ixa1, val in enumerate(absc1):
        dummy = val
        if tableType == 'DRIFT':
          betad = np.radians(val)
          gamma = 0
        elif tableType == 'YAW':
          betad = 0
          gamma = np.radians(val)

          
        pmmtyp = 1
        SepPoint = self.SeparationPoint
        coftyp, ucar = self.pmmmsm(uo, betad, gamma, pmmtyp, SepPoint)

        uoo = ucar
        ucar = self.pmmcar(uoo, betad, coftyp)
        tertiary = delta = heel = epsil = 0.0  

        udim, vdim, rdim, qdim, pdim, ddim = self.pmmmot(ucar, betad, gamma, delta, heel, epsil)
        speed_dict = self.hluref(udim, vdim, rdim, pdim)
        # update the multipliers with speed values if they are available
        multipliers = {key: self.updateMultiplierWithActualSpeed(value, speed_dict) for key, value in multipliers.items()}

        toh = 0
        uo = self.serviceSpeed
        fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)
        
        if icoty == 0: #BaseTable
            for key, multiplier in multipliers.items():
                if key == 'K_HL':
                  continue
                force_coefficients[key][dummy] = self.getForceCoefficient(multiplier, fdimu2[forceIndex[key]], utot2)
        else: ## now we create a correction table
            FBASE = fdimu2
            for ix, toh in enumerate(absc2):
                if ix == 0:
                    for key in correction_table.keys():
                        correction_table[key][ixa1, ix] = 1.0
                    continue
                fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)
                for key in correction_table.keys():
                  if key == 'K_HL':
                    continue
                  correction_table[key][ixa1, ix] = fdimu2[forceIndex[key]] / FBASE[forceIndex[key]]
      return force_coefficients, correction_table

  def pmms(self, pmmcoefs, motion, icoty, uo, ipmms):
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

      absc1, absc2 = self.defval(motion, icoty)

      if motion in [1, 2]: # Drift or Yaw
          tableType = 'DRIFT' if motion == 1 else 'YAW'
          force_coefficients, correction_table = self.calculate_force_coefficients(motion, tableType,  absc1, absc2, icoty, uo, pmmcoefs)
          if icoty == 0:  # construct baseTables
              tables = {}
              for key, value in force_coefficients.items():
                df = pd.DataFrame(index=value.keys(), data=value.values())
                df.index.name = 'BETAD' if motion == 1 else 'GAMMA'
                df = df.rename(columns={0: key})
                if key== 'X_HL' and tableType == 'DRIFT':
                  df.loc[-90] = df.loc[90.0] = 0.0
                  # Change sign of all values where index np.abs(BETAD) > 170
                  df.loc[np.abs(df.index) >= 170.0, 'X_HL'] *= -1
                  # These comments and the following code transfered and translated from fortran code    
                  #cbla    the next if block are made to make the X_HL for drift look right,
                  #cbla    meaning that we multiply the value at 20 degree driftangle with 1.5 
                  #cbla    to get the value at 45 degree and make the table symmetrical. 
                  #cbla    The 1.5 factor is taken from ship3005. 
                  #index135 = np.where(np.abs(df['BETAD'].values) == 135.0)   
                  #index70 =  np.where(np.abs(df['BETAD'].values) ==  70.0) 
                  #index45 =  np.where(np.abs(df['BETAD'].values) ==  45.0)     
                  #index20 =  np.where(np.abs(df['BETAD'].values) ==  20.0) 
                  df.loc[-135.0]  = -1.5* df.loc[-20.0]
                  df.loc[-70.0]   = df.loc[-20.0]
                  df.loc[-45.0]   = 1.5* df.loc[-20.0]
                  df.loc[ 45.0]   = 1.5* df.loc[-20.0]
                  df.loc[ 70.0]   = df.loc[-20.0]
                  df.loc[ 135.0]  = -1.5* df.loc[-20.0]                    
                if key=='X_HL' and tableType == 'YAW' :
                  print('Can not see how this works in the fortan code but Y and X (gamma) seems to be zero so hardcoded here !!')
                  df.loc[np.abs(df.index) == 90.0] = 0.0
                    
                tables[f'{key} = {tableType}'] = df       
          else: # construct correction tables
              tables = {}
              for key in force_coefficients.keys():
                  df = pd.DataFrame(correction_table[key], index=absc1, columns=absc2)
                  tables[f'{key}(BETAD,TOH)' if motion == 1 else f'{key}(GAMMA,TOH)'] = df
      elif motion == 3:
        for ixa1, val in enumerate(absc1):
          
          pmmtyp = 1
          SepPoint = self.SeparationPoint
          coftyp, ucar = self.pmmmsm(uo, betad, gamma, pmmtyp, SepPoint)
          uoo = ucar
          ucar = self.pmmcar(uoo, betad, coftyp)
          tertiary = delta = heel = epsil = 0.0  

          udim, vdim, rdim, qdim, pdim, ddim = self.pmmmot(ucar, betad, gamma, delta, heel, epsil)
          speed_dict = self.hluref(udim, vdim, rdim, pdim)
          # update the multipliers with speed values if they are available
          multipliers = {key: self.updateMultiplierWithActualSpeed(value, speed_dict) for key, value in multipliers.items()}

          toh = 0
          uo = self.serviceSpeed
          fdim, fdimu2, utot2 = self.pmmfor(pmmcoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh)        
          F1 = fdimu2
          idum = 0
        
        
      return tables

  

  def getForceTables(self,pmmcoefs,tableName):
    '''
    returns a dict of baseTables given a set of pmmCoefficients
    
    '''
    baseTable = None
    correctionTable = None
    if tableName == 'DRIFT':
      motion = 1
      icoty = 0
      ipmms = 1
      uo = self.serviceSpeed
      #Get baseTables
      baseTables = self.pmms(pmmcoefs,motion,icoty,uo,ipmms)
   
      
      #Get correction table
      icoty = 1
      correctionTables = self.pmms(pmmcoefs,motion,icoty,uo,ipmms)
    elif tableName == 'YAW':
      motion = 2
      icoty = 0
      ipmms = 1
      uo = self.serviceSpeed
      baseTables = self.pmms(pmmcoefs,motion,icoty,uo,ipmms)
      #Get correction table
      icoty = 1
      correctionTables = self.pmms(pmmcoefs,motion,icoty,uo,ipmms)
    return baseTables,correctionTables
    

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
  def getSwayDamping(self):
    return self.SwayDamping
  def getRollDamping(self):
    return self.RollDamping
  def getPitchDamping(self):
    return self.PitchDamping
  def getHeaveDamping(self):
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


        
    return tables

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
    BTJ  assume the double loop arround the minimize function is to try to ensure that a global minimum is found
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
      ###
      # pmms change -95,95 to -135,135 here BTJ changed this in the defval
      #             -85,85 to   -70,70
      ###
      absc1=[-180.0,-177.0,-175.0,-170.0,-135.0,-90.0,
              -70.0,-45.0,-20.0,-18.0,-15.0,-10.0,-5.0,-2.0,-1.0,
                0.0,1.0,2.0,5.0,10.0,15.0,18.0,20.0,45.0,70.0,90.0,
                135.0,170.0,175.0,177.0,180.0]
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
    if np.abs(np.abs(gamma)-np.radians(90.)) < .01:
        ucar = 0.
    return coftyp,ucar   
      
      
      
      
  def calculate_PMM_force(self, derivatives, pmmCoefs, acoef, bcoef, U, V, R, toh):
    coeff = np.array([])
    speedvector = np.array([])
    TOHCorrection = np.array([])
    for val in derivatives:
      if val in pmmCoefs.keys() and 'DOT' not in val:
        if 'VIVI' in val:
          sfactor = self.speedfactor(val, U, V, R)
          sfactor *= np.sign(U)    
          speedvector = np.append(speedvector, sfactor)
        else:       
          speedvector = np.append(speedvector, self.speedfactor(val, U, V, R))   
           
        coeff = np.append(coeff, pmmCoefs[val])
        a = acoef[val] if val in acoef.keys() else 0
        b = bcoef[val] if val in bcoef.keys() else 1
        TOHCorrection = np.append(TOHCorrection, (1 + a * toh**b))
    return np.sum(coeff * TOHCorrection * speedvector)

  def pmmfor(self, pmmCoefs, coftyp, uo, udim, vdim, rdim, qdim, pdim, ddim, toh):
      acoef = self.acoef
      bcoef = self.bcoef
      fdim = np.array([0., 0., 0.])
      fdimu2 = np.array([0., 0., 0.])
      lpp = self.Lpp
      utot2 = udim*udim + vdim*vdim if coftyp == 1 else 1.
      utot = np.sqrt(utot2)
      U = (udim - uo)/utot if coftyp == 1 else udim
      V = vdim/utot if coftyp == 1 else vdim
      R = rdim/(utot/lpp) if coftyp == 1 else rdim*lpp

      derivatives = ['X0', 'XU', 'XUU', 'XUUU', 'XV', 'XVV', 'XR', 'XRR', 'XVR', 'XUDOT', 'XRDOT', 'XUIUI']
      fdimu2[0] = self.calculate_PMM_force(derivatives, pmmCoefs, acoef, bcoef, U, V, R, toh) * 0.5 * self.rho * lpp**2
      fdim[0] = fdimu2[0]*utot2

      derivatives = ['Y0', 'Y0U', 'YUU', 'YV', 'YVV', 'YVIVI', 'YVVV', 'YVU', 'YR', 'YRR', 'YRIRI', 'YRRR', 'YRIVI', 'YRVV', 'YVIRI', 'YVRR', 'YVDOT', 'YRDOT']
      fdimu2[1] = self.calculate_PMM_force(derivatives, pmmCoefs, acoef, bcoef, U, V, R, toh) * 0.5 * self.rho * lpp**2
      fdim[1] = fdimu2[1]*utot2

      derivatives = ['N0', 'N0U', 'NUU', 'NV', 'NVV', 'NVIVI', 'NVVV', 'NVU', 'NR', 'NRR', 'NRIRI', 'NRRR', 'NRIVI', 'NRVV', 'NVIRI', 'NVRR', 'NVDOT', 'NRDOT ']
      fdimu2[2] = self.calculate_PMM_force(derivatives, pmmCoefs, acoef, bcoef, U, V, R, toh) * 0.5 * self.rho * lpp**3
      fdim[2] = fdimu2[2]*utot2

      return fdim, fdimu2, utot2
   
  
  def getRollHeelTables(self):
    tables = {name: {} for name in ['ROLL', 'DRIFT_BETAD', 'DRIFT_HEEL', 'YAW_DRIFT', 'YAW_HEEL', 'YAW']}
    files = glob.glob(r"h:\GitRepos\ShipYard2\data\hull\*.dat")

    for fil in files:
        filename = os.path.basename(fil)
        key = filename[-8:-4]  
        tableName = filename[:-9]

        if 'ROLL' in fil:
            index_col = 'EPSI'
        elif tableName in ['DRIFT_HEEL', 'DRIFT_BETAD', 'YAW_DRIFT']:
            index_col = 'BETAD'
        elif tableName.startswith('YAW'):
            index_col = 'GAMMA'
        else:
            continue

        df = pd.read_csv(fil, header=0, sep='\s+', index_col=index_col)
        tables[tableName][key] = df

    return tables
  
    

        
        
# %%

if __name__ == '__main__':
  import matplotlib.pyplot as plt
  import matplotlib.ticker as ticker
  from mpl_toolkits.mplot3d import Axes3D
  import numpy as np  
  import plotly.graph_objects as go
  shipDatadict={}
  shipDatadict['shipnr'] = 3949
  shipDatadict['Lpp'] = 340.5
  shipDatadict['Beam'] = 53.5

  shipDatadict['waterPlaneArea'] = 0.893628117 * shipDatadict['Lpp'] * shipDatadict['Beam']
  shipDatadict['propellerType'] ='FP'
  shipDatadict['PropellerPitch'] = 0.71505
  shipDatadict['displacement'] =218220.0
  shipDatadict['propellerDiameter'] = 7.0
  shipDatadict['draftAft']  = 17.
  shipDatadict['draftFore'] = 17.0
  shipDatadict['blockCoefficient'] = 0.704661757
  shipDatadict['meanDraft'] = (shipDatadict['draftAft'] + shipDatadict['draftFore'] )/2.0
  shipDatadict['wettedSurface'] = 0.7801584*(shipDatadict['Lpp'] *(2*shipDatadict['meanDraft']+shipDatadict['Beam']))
  shipDatadict['underWaterLateralArea'] =  0.9843655* shipDatadict['Lpp'] * shipDatadict['meanDraft']
 
  shipDatadict['CenterofGravity'] = np.array([-0.002290749, 0, -0.415058824]) * np.array([shipDatadict['Lpp'] ,shipDatadict['Beam'],shipDatadict['meanDraft']])
  shipDatadict['verticalCenterOfBoyancy'] = 0.453647058824 * shipDatadict['meanDraft'] 
  shipDatadict['GyrationArms'] = np.array([0.4, 0.25, 0.25]) * np.array([shipDatadict['Beam'],shipDatadict['Lpp'],shipDatadict['Lpp']])
  shipDatadict['SeparationPoint'] = 45.0
  hull_Thumbs = HullThumbs(shipDatadict)
# %%
  ### Resistance table
 
  SY1_XL_FN = [
    (-1.878140E-01, 2.104312E-03),    (-1.609834E-01, 1.941680E-03),    (-1.341528E-01, 1.899358E-03),
    (-1.073222E-01, 1.928490E-03),    (-8.049170E-02, 2.293411E-03),    (-5.366113E-02, 3.774366E-03),
    (-2.683057E-02, 6.356310E-03),    (-1.730242E-04, 9.845336E-03),    (0.000000E+00, 0.000000E+00),
    (1.730242E-04, -8.950305E-03),    (2.039914E-02, -6.481410E-03),    (4.062525E-02, -4.450279E-03),
    (6.085137E-02, -2.970214E-03),    (8.107748E-02, -2.067285E-03),    (1.013036E-01, -1.764216E-03),
    (1.215297E-01, -1.733612E-03),    (1.417558E-01, -1.729059E-03),    (1.619819E-01, -1.768446E-03),
    (1.822081E-01, -1.871172E-03),    (2.024342E-01, -2.051965E-03),    (2.226603E-01, -2.318379E-03),
    (2.428864E-01, -2.668396E-03)
    ]
# Convert the data to a pandas DataFrame
  df = pd.DataFrame(SY1_XL_FN,columns=['FN','X_HL']) 
  
  shallow = False
  resistanceTable = hull_Thumbs.RDHOLTROP(shallow)
  fig, ax = plt.subplots()
  # Use scientific notation for y axis
  formatter = ticker.ScalarFormatter(useMathText=True)
  formatter.set_scientific(True) 
  formatter.set_powerlimits((-1,1)) 
  ax.yaxis.set_major_formatter(formatter)
  plt.plot(df['FN'],df['X_HL'],label='SY1')
  plt.plot(resistanceTable.index,resistanceTable['X_HL'],label='pySY2')
  plt.legend()
  plt.xlabel('FN')
  plt.title('X_HL(FN)')
  plt.grid()
  plt.show()
  idum = 0
  pass
# %%
  print("***  Shallow water Schlicting ****")
  shallowCorrection = hull_Thumbs.shallowWaterResistance()
  sy1Data=pd.read_csv(r'H:\GitRepos\ShipYard2\data\hull\FN_TOH_SY1_result.dat',header=0,sep='\s+')  
  sy1Data.set_index('FN',inplace=True)
  
  sy1_surface = go.Surface(z=sy1Data.values,x=sy1Data.index,y=sy1Data.columns, colorscale='Cividis', opacity=0.9, showscale=False)
  sy2_surface = go.Surface(z=shallowCorrection.values, x=shallowCorrection.index, y=shallowCorrection.columns,colorscale='Viridis')
  # Combine the surfaces into a single figure
  fig = go.Figure(data=[sy1_surface, sy2_surface])
  #fig = go.Figure(data=[go.Surface(z=shallowCorrection.values, x=shallowCorrection.index, y=shallowCorrection.columns)])
  # Update layout
  fig.update_layout(title='X_HL(FN,TOH)', autosize=False,
                    width=500, height=500,
                    margin=dict(l=65, r=50, b=65, t=90))

  # Show the plot
  fig.show()  
# %%
  
  # !!! dmimix takes a bit of time !!!
  #a,b =hull_Thumbs.dmimix()  
  with open(r'data\PMMdata\3949NeuPMMDe_005_001.txt','r') as f:
    pmmCoefs = {}
    for line in f:
      splt = line.split()
      pmmCoefs[splt[0]] = float(splt[1])/1.0E5
  #get base and correctionTables for Drift
  tableType =  'DRIFT'
  baseTables,correctionTables = hull_Thumbs.getForceTables(pmmCoefs,tableType)
  # plot results
  print(f"############ DRIFT X;Y;N and corresponding TOH tables are ready")
  dfSY1 = pd.read_csv(r'H:\GitRepos\ShipYard2\data\PMMdata\SY1hull3949XHL_DriftBetaD.dat',header=None,sep='\s+')
  fig, ax = plt.subplots()
  plt.plot(baseTables['X_HL = DRIFT'],'-*',label='SY2:X_HL(BETAD)')
  plt.plot(dfSY1[0],dfSY1[1],'-*',label='SY1')
  plt.legend()
  plt.grid()
  plt.show()
  pass    
# %%

  fig = go.Figure(data=[go.Surface(z=correctionTables['X_HL(BETAD,TOH)'].values, x=correctionTables['X_HL(BETAD,TOH)'].index, y=correctionTables['X_HL(BETAD,TOH)'].columns)])
  # Update layout
  fig.update_layout(title='X_HL(BETAD,TOH)', autosize=False,
                    width=500, height=500,
                    margin=dict(l=65, r=50, b=65, t=90))

  # Show the plot
  fig.show()
  ##
   
           
  # %%
  tableType = 'YAW'
  baseTables,correctionTables = hull_Thumbs.getForceTables(pmmCoefs,tableType)
 
  SY1_gamma = [-9.000000E+01, -2.600000E+01, -2.400000E+01, -2.200000E+01, -2.000000E+01, -1.500000E+01,
                -1.000000E+01, -5.000000E+00, -2.000000E+00, -1.000000E+00,  0.000000E+00,  1.000000E+00,
                2.000000E+00,  5.000000E+00,  1.000000E+01,  1.500000E+01,  2.000000E+01,  2.200000E+01,
                2.400000E+01,  2.600000E+01,  9.000000E+01]
  SY1_X_HL =[  0.000000E+00, -7.491886E-04, -6.449616E-04, -5.470891E-04, -4.560481E-04, -2.611558E-04,
              -1.175568E-04, -2.961415E-05, -4.748381E-06, -1.187457E-06,  0.000000E+00, -1.187457E-06,
              -4.748381E-06, -2.961415E-05, -1.175568E-04, -2.611558E-04, -4.560481E-04, -5.470891E-04,
              -6.449616E-04, -7.491886E-04,  0.000000E+00]

  #print(baseTables['X_HL = YAW'].head())
  fig,ax = plt.subplots()
  plt.plot(baseTables['X_HL = YAW'],'-*',label='X_HL(GAMMA)')
  plt.plot(SY1_gamma,SY1_X_HL,'-*',label='SY1_XHL(Gamma)')
  plt.legend()
  plt.xlabel('GAMMA')
  plt.title('X_HL(GAMMA)')
  plt.grid()
  plt.show()
  pass

# %%
  
  tables = hull_Thumbs.getRollHeelTables()
  #fig,ax = plt.subplots()
  tables['ROLL']['Y_HL'].plot()
  plt.title('ROLL')
  plt.grid()
  idum = 0

  # %%
  with open(r"H:\GitRepos\ShipYard2\data\PMMdata\3949NeuPMMDe_005_001.txt",'r') as f:
    pmmCoefs = {}
    for line in f:
        splt = line.split()
        pmmCoefs[splt[0]] = float(splt[1])/1.0E5   
  motion = 3
  icoty = 0                
  uo = hull_Thumbs.serviceSpeed
  ipmms = 1
  tableType ='YAW_DRIFT'
  absc1,absc2 = hull_Thumbs.defval(motion,icoty)
  tables,correctionTables = hull_Thumbs.yaw_drift_Tables(motion, tableType, absc1, absc2, icoty, uo, pmmCoefs)

  z= tables['X_HL'].values
  x = tables['X_HL'].index
  y = tables['X_HL'].columns
  surface=go.Surface(z=z,x=x,y=y,contours={
        "z": {
            "show": True,
            "start": np.min(z),
            "end": np.max(z),
            "size": 2,
            "color": "white"
        }})
  fig = go.Figure(data=[surface])
  # Update layout
  fig.update_layout(title='YAW_DRIFT(X_HL)', autosize=False,
                    width=500, height=500,
                    margin=dict(l=65, r=50, b=65, t=90))

  # Show the plot
  fig.show()         

# %%
