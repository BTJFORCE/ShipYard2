


import numpy as np


from scipy.constants import g


import Thumbs as Thumbs


class HullThumbs(Thumbs.Thumbs):
    


    def __init__(self,shipnr,propellerDiameter=None):


      super().__init__(shipnr,propellerDiameter=None)


#


# FNH values for squat calculation


#


      self.FNH_values = [-1.90,-1.20,-1.10,-1.05,-1.00, 


                    -0.95,-0.90,-0.80,-0.70,-0.60,


                    -0.50,-0.40,-0.30,-0.20,-0.10, 0.00,


                     0.10, 0.20, 0.30, 0.40, 0.50,


                     0.60, 0.70, 0.80, 0.90, 0.95,


                     1.00, 1.05, 1.10, 1.20, 1.90]




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

    


