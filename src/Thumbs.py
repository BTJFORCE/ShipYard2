
import numpy as np

class Thumbs:
    '''
    Base class holding ship particulars used in different thumps
    some data are hard coded as SY2 should have these variables already
    methods needed to set these values from SY2
    a new extra data fuile introduced with data normaly inserted by user in SY1
    '''


    def __init__(self,shipDatadict):
 
        self.read_SY1Data(shipDatadict['shipnr'])
        print('hardcoded some values here for now')
        self.lengthWaterLine = self.ship_data['Length_of_Waterline'] 
        self.Lpp = shipDatadict['lpp'] ##312 # From ShipData
        self.beam = shipDatadict['Beam'] ##53.5  ## get from shipDataLib
        self.pRatio = self.ship_data['Propeller_Area_Ratio']
        self.wettedSurface = shipDatadict['wettedSurface'] ##23243.8  ## get from shipDataLib
        self.waterPlaneArea = shipDatadict['waterPlaneArea'] ##0.893628117 * self.Lpp * self.beam
        self.underWaterLateralArea = shipDatadict['underWaterLateralArea']
        self.propellerType = shipDatadict['propellerType'] #'FP'
        self.PD = shipDatadict['PropellerPitch'] ## 0.71505  ## get from shipDataLib
        self.nrProp = self.ship_data['number_of_propellers']
        self.nrBlades = self.ship_data['number_of_blades']
        self.displacement = shipDatadict['displacement'] ##218220.0 ## get from shipDataLib
        self.serviceSpeed = self.ship_data['service_speed']
        self.propeller_revs_serviceSpeed = self.ship_data['propeller_revs_serviceSpeed']
        self.propellerDiameter = shipDatadict['propellerDiameter'] #7 ## get from shipData
        self.draftAft = shipDatadict['draftAft'] #17. ## get from shipData
        self.draftFore = shipDatadict['draftFore'] #17. ## get from shipData
        self.nozzle = 0
        self.meanDraft =(self.draftAft + self.draftFore)/2.0
        self.waterLineBlock = self.displacement /(self.lengthWaterLine*self.beam*self.meanDraft) # waterline block coefficient
        self.blockCoefficient = shipDatadict['blockCoefficient'] #0.704661757
        self.COG=shipDatadict['CenterofGravity'] #np.array([-0.002290749, 0, -0.415058824]) * np.array([self.Lpp ,self.beam,self.meanDraft])
        self.verticalCenterOfBoyancy = shipDatadict['verticalCenterOfBoyancy'] #'0.453647058824 * self.meanDraft
        self.gyration = shipDatadict['GyrationArms'] ##np.array([0.4, 0.25, 0.25]) * np.array([self.beam,self.Lpp,self.Lpp])
        self.SeparationPoint = shipDatadict['SeparationPoint']
       
        self.Formal = self.ship_data['Total_FormFactor'] ## 1.21268
        self.viscocity = self.ship_data['viscocity'] ##1.34066
        self.prismaticCoefficeint = self.ship_data['Prismatic_Coefficient'] ##0.70699
        self.midshipSection = self.ship_data['midshipSection'] ##0.98
        self.Abulb = self.ship_data['BulbArea'] ##0.004
        self.verticalCenterBulb = self.ship_data['verticalCenterBulb']## 5.0 #Vertical position of center of bulb from BL, m 
        self.LCB_ratio = self.ship_data['LCB_ratio'] ##0.00612186
        self.Fmxaft = self.ship_data['FormFactorStern']

        self.GMT = self.ship_data['GMT']
        self.GML = self.ship_data['GML']
        self.rho = 1025.0


        #self.areaRatio = aRatio
       # self.Displacement = 100000
        #self.Cp = self.Displacement / (self.lengthWaterLine * self.crossSectionArea) #
    def read_SY1Data(self, shipnr):
        self.ship_data = {}
        # Assuming the data is in a file called 'data.txt'
        with open(f'U:/ships/ship{shipnr}/ShipYardData{shipnr}.dat', 'r') as file:
            for line in file:
                line = line.strip()
                if line.find("##") == 0 or len(line)==0:
                    continue
                # Split the line into key and value
                key, value = line.split()
                # Convert the value to a float and add it to the dictionary
                self.ship_data[key] = float(value)