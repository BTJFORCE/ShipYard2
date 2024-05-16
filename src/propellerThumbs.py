

class propellerThumbs:


    def __init__(self,shipnr,propellerDiameter=None):
 
        self.read_SY1Data(shipnr)
        print('hardcoded some values here for now')
        self.lengthWaterLine = self.ship_data['Length_of_Waterline'] 
        self.beam = 53.5  ## get from shipDataLib
        self.pRatio = self.ship_data['Propeller_Area_Ratio']
        self.wettedSurface = 23243.8  ## get from shipDataLib
        self.PD = 0.71505  ## get from shipDataLib
        self.nrProp = self.ship_data['number_of_propellers']
        self.nrBlades = self.ship_data['number_of_blades']
        self.displacement = 218220.0 ## get from shipDataLib
        self.serviceSpeed = self.ship_data['service_speed']
        self.propeller_revs_serviceSpeed = self.ship_data['propeller_revs_serviceSpeed']
        if propellerDiameter:
            self.propellerDiameter = propellerDiameter ## get from shipDataLib
        else:
            self.propellerDiameter = 7 ## get from shipData
        self.draftAft = 17. ## get from shipData
        self.draftFore = 17. ## get from shipData
        self.nozzle = 0
        self.meanDraft =(self.draftAft + self.draftFore)/2.0
        self.waterLineBlock = self.displacement /(self.lengthWaterLine*self.beam*self.meanDraft) # waterline block coefficient
    
       
        self.Formal = self.ship_data['Total_FormFactor'] ## 1.21268
        self.viscocity = self.ship_data['viscocity'] ##1.34066
        self.prismaticCoefficeint = self.ship_data['Prismatic_Coefficient'] ##0.70699
        self.midshipSection = self.ship_data['midshipSection'] ##0.98
        self.Abulb = self.ship_data['BulbArea'] ##0.004
        self.verticalCenterBulb = self.ship_data['verticalCenterBulb']## 5.0 #Vertical position of center of bulb from BL, m 
        self.LCB_ratio = self.ship_data['LCB_ratio'] ##0.00612186
        self.Fmxaft = self.ship_data['FormFactorStern']

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