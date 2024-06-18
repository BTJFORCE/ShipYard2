

# %%
import numpy as np
from scipy.constants import g
from scipy.optimize import fmin_bfgs,minimize

import pandas as pd

import Thumbs as Thumbs


def FUNCAB(XAB,dfTOH,index):
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



    def linearDamping(self):
      RHO=1026.0  
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


    def S_mid(self,CB, T, L, FNH):
        rVal = 0.01 * ((38 * CB * T / L) * FNH**2) / (np.sqrt(1 - FNH**2))
        return rVal

    def trim(self,CB, T, L, FNH):
        rVal = 0.01 * (47.4 * CB * T / L - 1.2) * FNH**2 / (np.sqrt(1 - FNH**2))
        return rVal


    def calculateSquat(self):
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
    
    def getHeaveDamping(self):
      return self.HeaveDamping

    def getPitchDamping(self):
      return self.PitchDamping

    def getRollDamping(self):
      return self.RollDamping

    def getSurgeDamping(self):
      return self.SurgeDamping

    def getSwayDamping(self):
      return self.HeaveDamping

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
      
      
    def pmmfor(self,pmmdata,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,
            		toh,rho,lpp,fdimu2,utot2,fdim):
    
      fdim=np.array[0.,0.,0.]
      fidmu2=np.array[0.,0.,0.]
      ## SKIPPING D,Q,P derivatives as we dont have them
      xderivatives=['X0','XU','XUU','XUUU','XV','XVV','XR','XRR','XVR','XUDOT','XRDOT','XUIUI']
      yderivatives=['Y0','Y0U','YUU','YV','YVV','YVIVI','YVVV','YVU','YR','YRR','YRIRI','YRRR','YRIVI','YRVV','YVIRI','YVRR','YVDOT','YRDOT']
      nderivatives=['N0','N0U','NUU','NV','NVV','NVIVI','NVVV','NVU','NR','NRR','NRIRI','NRRR','NRIVI','NRVV','NVIRI','NVRR','NVDOT','NRDOT ']
      DERPOINT=['IYVDOT','IYRDOT','IYV','IYVL','IYR','IYVIVI','IYRIRI','IYRRR','IYVRR','IYVIRI','IYRIVI','IYRVV','INVDOT','INRDOT','INV',
                'INVL','INR','INVIVI','INVVV','INVRR','INVIRI','INRVV','INRIVI','IXUDOT','IXVR','IXVV','IXRR']
      acoef,bcoef = self.dmimix()


      pmmCoefs={}
      with open(pmmdata) as f:
        for line in f:
          key,value = line.split()
          pmmCoefs[key] = float(value)
          
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
        BETAMSMFQ=np.sign(1.0,udim)*(-np.asin(vdim/np.sqrt(vdim**2+udim**2)))
        usign = np.sign(1.0,U)
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
              continue  # in the fortran PMMFOR not DOT derivatives appears
            xcoeff=np.append(xcoeff,pmmCoefs[val])
            speedvector=np.append(speedvector,self.speedfactor(val,U,V,R))
            TOHCorrection=np.append(TOHCorrection,(1+ acoef[val]*toh**bcoef[val]))
            
        fdimu2[0]= np.sum(xcoeff * TOHCorrection * speedvector) * 0.5 * 1025 * self.Lpp**2
        fdim[0] = fdimu2[0]*utot2
        # Y-force
        ycoeff=np.array([])
        speedvector=np.array([])
        TOHCorrection=np.array([])
        for val in yderivatives:
          if val in pmmCoefs.keys():
            if 'DOT' in val:
              continue  # in the fortran PMMFOR not DOT derivatives appears
            ycoeff=np.append(ycoeff,pmmCoefs[val])
            speedvector=np.append(speedvector,self.speedfactor(val,U,V,R)) 
            TOHCorrection=np.append(TOHCorrection,(1+ acoef[val]*toh**bcoef[val]))       
        
        fdimu2[1]= np.sum(ycoeff * TOHCorrection * speedvector) * 0.5 * 1025 * self.Lpp**2
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
            TOHCorrection=np.append(TOHCorrection,(1+ acoef[val]*toh**bcoef[val]))        
        
        fdimu2[2]= np.sum(ncoeff * TOHCorrection * speedvector) * 0.5 * 1025 * self.Lpp**3
        fdim[2] = fdimu2[2]*utot2      
        
        return fdim,fdimu2      
        
        
# %%
if __name__ == '__main__':
  
  shipDatadict={}
  shipDatadict['shipnr'] = 3949
  shipDatadict['Lpp'] = 312
  shipDatadict['Beam'] = 53.5
  shipDatadict['wettedSurface'] = 23243.8
  shipDatadict['waterPlaneArea'] = 0.893628117 * shipDatadict['Lpp'] * shipDatadict['Beam']
  shipDatadict['propellerType'] ='FP'
  shipDatadict['PropellerPitch'] = 0.71505
  shipDatadict['displacement'] =218220.0
  shipDatadict['propellerDiameter'] = 7.0
  shipDatadict['draftAft']  = 17.
  shipDatadict['draftFore'] = 17.0
  shipDatadict['meanDraft'] = (shipDatadict['draftAft'] + shipDatadict['draftFore'] )/2.0
  shipDatadict['blockCoefficient'] = 0.704661757
  shipDatadict['CenterofGravity'] = np.array([-0.002290749, 0, -0.415058824]) * np.array([shipDatadict['Lpp'] ,shipDatadict['Beam'],shipDatadict['meanDraft']])
  shipDatadict['verticalCenterOfBoyancy'] = 0.453647058824 * shipDatadict['meanDraft'] 
  shipDatadict['GyrationArms'] = np.array([0.4, 0.25, 0.25]) * np.array([shipDatadict['Beam'],shipDatadict['Lpp'],shipDatadict['Lpp']])

  hull_Thumbs = HullThumbs(shipDatadict)
  a,b =hull_Thumbs.dmimix()
  print(a)

           
          

        
        
        
        
        
      
    



# %%
TOH=[.1,.2,.5,.75,.833]
XAB=[2,2]
FUNCV = 0
FNCTOH=[1.01,1.03,1.27,1.8,3.49]
for ix in range(len(TOH)):
  FUNCV = FUNCV+  (XAB[0]*TOH[ix]**XAB[1] -(FNCTOH[ix]-1.0))**2
  print(FUNCV)
# %%
