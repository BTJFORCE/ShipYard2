# %%
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import propellerThumbs as pThumbs

'''
c
c  single screw ship
c
      Z7=B9/L9
      E7=Err
      IF (Ierr.eq.1) E7=.9922-.05908*Y+.07424*(C9+2.25*(O1-.015))
      W5=Wake
      IF (Iwake.ne.1) goto 1
      A6=.1*Z7+.149
      A7=.05*Z7+.449
      A8=585-5027*Z7+11700*Z7**2
      W1=A6+A7/(A8*(.98-C8)**3+1)+.05*Fmaft/(100*(C8-.7)**2+1)/2 
      W3=-.18+.00756/(Diam/L9+.002)
      IF (W3.gt.1d0) W3=.1d0
      W5=W1+W3
    1 T5=Thru
      IF (Ithru.ne.1) goto 100
      A6=.625*Z7+.08
      A7=.165-.25*Z7
      A8=825-8060*Z7+20300*Z7**2
      T5=A6+A7/(A8*(.98-C8)**3+1)-.01*Fmaft+2*(Diam/L9-.04)
      GOTO 100
'''


class wake(pThumbs.propellerThumbs):
           # Initialize xw1 with the given values
    xw1 = [-180, -165, -135, -90, -45, -20, 0, 
        10, 18.4, 45, 90, 110, 130, 140, 150, 160, 165, 170, 180]

    yw1= [0,0.67,0.8,1.00] # TOH's

    fw1 = [
    0.65, 0.80, 0.90, 1.06,
    0.59, 0.69, 0.76, 0.87,
    0.48, 0.53, 0.58, 0.63,
    0.26, 0.26, 0.26, 0.26,
    0.09, 0.17, 0.17, 0.24,
    0.18, 0.18, 0.30, 0.42,
    0.50, 0.60, 0.90, 1.02,
    0.81, 1.08, 1.33, 1.55,
    1.00, 1.35, 1.50, 1.82,
    1.20, 1.70, 1.87, 2.20,
    1.00, 1.60, 1.75, 2.05,
    0.90, 1.45, 1.57, 1.85,
    0.79, 1.24, 1.38, 1.70,
    0.75, 1.18, 1.28, 1.60,
    0.73, 1.11, 1.20, 1.40,
    0.68, 1.00, 1.09, 1.30,
    0.66, 0.93, 1.05, 1.20,
    0.64, 0.87, 1.00, 1.15,
    0.61, 0.80, 0.90, 1.06
    ]
    xw2=[-180,-150,-130,-110,-90,-75,-65,-50,-45,-35,-15,-5,0,
         5,15,35,45,50,65,75,90,110,130,150,180]
    fw2 = [
    1.00,    1.00,
    0.85,    0.99,
    0.61,    0.95,
    0.17,    0.83,
    -0.26,   0.65,
    -0.50,   0.42,
    -0.40,   0.12,
    -0.07,  -0.48,
    0.03,   -0.50,
    0.25,    0.03,
    0.71,    0.88,
    0.90,    1.00,
    1.00,    1.00,
    1.00,    0.90,
    0.88,    0.71,
    0.03,    0.25,
    -0.50,   0.03,
    -0.48,  -0.07,
    0.12,   -0.40,
    0.42,   -0.50,
    0.65,   -0.26,
    0.83,    0.17,
    0.95,    0.61,
    0.99,    0.85,
    1.00,    1.00
    ]
    xw3 = [-60,-40,-30,-20,-10,0,10,20,30,40,60]
    fw3= [ .16,.09,.06,.03,.01,.0,.01,.03,.06,.09,.16]

    ft1 = [
    1.15, 1.60, 1.60, 2.00,
    2.15, 2.35, 3.00, 4.05,
    2.40, 2.60, 3.45, 4.75,
    2.00, 2.35, 3.10, 4.50,
    1.50, 1.80, 2.50, 3.55,
    0.95, 1.35, 1.85, 2.60,
    0.26, 0.43, 0.64, 0.95,
    0.80, 1.05, 1.20, 1.45,
    1.00, 1.20, 1.40, 1.70,
    1.00, 1.25, 1.47, 1.80,
    1.30, 1.60, 1.83, 2.20,
    1.53, 1.80, 2.00, 2.40,
    1.70, 1.80, 1.80, 1.00,
    1.70, 1.60, 0.30, -2.00,
    1.55, 1.25, -2.05, -3.20,
    0.90, -1.20, -2.10, -2.40,
    -0.40, -1.40, -1.40, -1.40,
    -0.30, -0.50, -0.50, -0.50,
    0.61, 0.80, 0.90, 2.00
]
    def __init__(self,shipnr,propellerDiameter=None):
        super().__init__(shipnr,propellerDiameter=None)

    def R1(self):
        L9 = self.lengthWaterLine 
        Visco = self.viscocity
        V5 = self.serviceSpeed
        if V5 <= 1:
            return L9/Visco*1000000
        else:
            return V5*L9/Visco*1000000
    
    def c6(self):
        R1 = self.R1()
        return 1000*.075/(np.log10(R1)-2)**2

    def Cx8(self):
        S = self. wettedSurface
        B= self.beam
        LWL = self.lengthWaterLine
        D = self.propellerDiameter
        TA = self.draftAft
        if B/TA <= 5:
            return S/(LWL*D)*(B/TA)
        else:
            return S*(7*B/TA - 25)/(LWL*D*(B/TA - 3))

    def Cx9(self):
        if self.Cx8() <= 28 :
            return self.Cx8()
        else:
            return 32 - 16/(self.Cx8() - 24)

    def Cx11(self):
        D = self.propellerDiameter
        TA = self.draftAft
        if TA/D <= 2:
            return TA/D
        else:
            return 0.0833333*(TA/D)**3 + 1.33333

    def Cx19(self):
        C9 = self.prismaticCoefficeint
        C3 = self.midshipSection
        C8 = self.waterLineBlock
        if C9 > 0.7:
            return .18567/(1.3571-C3)-.71276+.38648*C9
        else:
            return .12997/(.95-C8)-.11056/(.95-C9)

    def c20(self):
        Cstern = self.Cstern
        return 1 + 0.015*Cstern

    def Cp1(self) :
        O1 = self.LCB_ratio
        C9 = self.prismaticCoefficeint
        return 1.45*C9 -0.315 -2.25*O1

    def Cxa(self):
        B9 = self.beam
        T9 = self.meanDraft
        C3 = self.midshipSection
        Abulb = self.Abulb
        Fmfore = self.draftFore
        Xbulb=Abulb*B9*T9*C3
        K6r = self.verticalCenterBulb
        L9 = self.lengthWaterLine
        C8 = self.waterLineBlock
        Cx3=.56*Xbulb**1.5/(B9*T9*(.31*np.sqrt(Xbulb)+Fmfore-K6r))
        Cx2=(-1.89)*np.sqrt(Cx3)

        Cx2=np.exp(Cx2)
        Cx4=Fmfore/L9
        if Cx4 > .04:
             Cx4=.04
        Cxa=.006*(L9+100)**(-.16)-.00205 \
           +.003*np.sqrt(L9/7.5)*C8**4*Cx2*(.04-Cx4)
        return Cxa

    def Cv(self):
        Formal = self.Formal
        C6 = self.c6()
        Cxa = self.Cxa()
        return Formal*C6/1000+Cxa

    def wakeCalculation(self):
        L9 = self.lengthWaterLine
        Fmaft = self.draftAft
        C9 = self.prismaticCoefficeint
        T9 = self.meanDraft
        Cv = self.Cv()
        Cx11 = self.Cx11()
        Cp1 = self.Cp1()
        B9 = self.beam
        Cx19 = self.Cx19()
        Cx9 = self.Cx9()
        Fmxaft = self.Fmxaft
        Inrpro = self.nrProp
        if Inrpro == 2 :
            wake = .3095*C8+10*Cv*C8-.23*Diam/SQRT(B9*T9)
        else:
            wake = (1+.05*Fmxaft)* \
                   (Cx9*Cv*L9/Fmaft*(.050776+.93405*Cx11*Cv/(1.0-Cp1)) \
                    +.27915*np.sqrt(B9/L9/(1.0-Cp1))+Cx19)
        self.wake = wake
        return wake

    def thrustDeduction(self):
        B9 = self.beam
        L9 = self.lengthWaterLine
        T9 = self.meanDraft
        C9 = self.prismaticCoefficeint
        Diam = self.propellerDiameter
        O1 = self.LCB_ratio
        Fmxaft = self.Fmxaft
        T5=.25014*(B9/L9)**.28956*(np.sqrt(B9*T9)/Diam)**.2624/   \
                         (1-C9+2.25*O1)**.01762+.005*Fmxaft
        if self.nrProp == 2:
            T5 = .325*C8-.1885*Diam/np.sqrt(B9*T9)
        self.thrustDeduction = T5
        return self.thrustDeduction
    
    def relativeRotativeEfficiency(self):
        Y = self.pRatio
        C9 = self.prismaticCoefficeint
        O1 = self.LCB_ratio
        X = self.PD

        E7=.9922-.05908*Y+.07424*(C9-2.25*O1)
        if self.nrProp == 2:
            E7=.9737+.111*(C9-2.25*O1)-.06325*X
        return E7

    
        
    def wakeTable(self):
        fw1=np.reshape(self.fw1,(19,4))
        wakeTable = fw1[:,0]*self.wake
        print ('Adjusting wake ')
        betap = np.arctan2(self.serviceSpeed,.7*3.141592654*self.propeller_revs_serviceSpeed*self.propellerDiameter)
        betap = np.degrees(betap)
        f = interp1d(self.xw1,wakeTable,kind='linear')
        fval = f(betap)
        fac= self.wake/fval
        for ix in range(len(self.xw1)):
            print (f"{self.xw1[ix]} {wakeTable[ix]*fac}")
    
    def wakeTablePDIR(self):
        fw2 = np.reshape(self.fw2,(len(self.xw2),2))

        print ("WAKE-PROPL-PDIR")
        for ix in range(len(self.xw2)):
            print(f"{self.xw2[ix]}  {fw2[ix][0]} {fw2[ix][1]}")

    def thrustDeductionTable(self):
        ft1 = np.reshape(self.ft1,(19,4))
        thrustDeductionTable = ft1[:,0]* self.thrustDeduction
        for ix in range(len(self.xw1)):
            print (f"{self.xw1[ix]} {thrustDeductionTable[ix]}")

    def wakeDeltaRAS(self):
        print ("WAKEF-DELTA-RAS")
        for ix in range(len(self.xw3)):
            print(f"{self.xw3[ix]} {self.fw3[ix]}")

# %%
        



if __name__ == '__main__':
    propellerDiameter = 7
    nProps = 1
    ship_nr = 3949
    aWake = wake(ship_nr,propellerDiameter)
    wake  = aWake.wakeCalculation()
    thrde = aWake.thrustDeduction()
    print (f" Wake for D ={propellerDiameter} = {wake}")
    aWake.wakeTable()
    #print(f" Thrust deduction for D = {propellerDiameter}  , {thrde}")
    #aWake.thrustDeductionTable()
    #aWake.wakeTablePDIR()
    pass
# %%
    xw1= aWake.xw1
    fw1 = np.reshape(aWake.fw1,(19,4))
    f = interp1d(xw1,fw1[:,0],kind='cubic')
    spx = np.arange(-180,180,0.5)
    splined=f(spx)
    plt.plot(xw1,fw1[:,0],'-*',label='linear')
    plt.plot(spx,splined,label='splined')
    plt.legend()
# %%
