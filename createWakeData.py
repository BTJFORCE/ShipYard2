import numpy as np

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


class wake:
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
    def __init__(self,propellerDiameter,nrprop):
 

        lengthWaterLine = 346.3 
        beam = 53.5
        pRatio = .98
        wettedSurface = 23243.8
        self.nrProp = nrprop
        self.displacement = 218220.0
        self.wettedSurface = wettedSurface
        self.propellerDiameter = propellerDiameter
        self.lengthWaterLine = lengthWaterLine 
        self.beam = beam
        self.draftAft = 17.
        self.draftFore = 17
        self.meanDraft =(self.draftAft + self.draftFore)/2.0
        self.waterLineBlock = self.displacement /(self.lengthWaterLine*self.beam*self.meanDraft) # waterline block coefficient
        print('hardcoded some values here for now')
        self.speed = 10.8548
        self.Formal = 1.21268
        self.viscocity = 1.34066
        self.prismaticCoefficeint = 0.70699
        self.midshipSection = 0.98
        self.Abulb = 0.004
        self.verticalCenterBulb = 5.0 #Vertical position of center of bulb from BL, m 
        self.LCB_ratio = 0.00612186
        self.Fmxaft = 0

        #self.areaRatio = aRatio
       # self.Displacement = 100000
        #self.Cp = self.Displacement / (self.lengthWaterLine * self.crossSectionArea) #
    def R1(self):
        L9 = self.lengthWaterLine 
        Visco = self.viscocity
        V5 = self.speed
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
        
    def wakeTable(self):
        fw1=np.reshape(self.fw1,(19,4))
        wakeTable = fw1[:,0]*self.wake

        for ix in range(len(self.xw1)):
            print (f"{self.xw1[ix]} {wakeTable[ix]}")

        



if __name__ == '__main__':
    propellerDiameter = 10
    nProps = 1
    aWake = wake(propellerDiameter,nProps)
    wake  = aWake.wakeCalculation()
    print (f" Wake for D ={propellerDiameter} = {wake}")

    aWake.wakeTable()
    pass