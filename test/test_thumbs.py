# %%
import unittest
import sys
import numpy as np
from scipy.constants import g
from scipy.optimize import fmin_bfgs,minimize
from scipy.interpolate import interp1d
import json
import os
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))


import KtKq
from HullThumbs import HullThumbs
from propellerSteering import propellerStering

# %%

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
shipDatadict['meanDraft'] = (shipDatadict['draftAft'] + shipDatadict['draftFore'] )/2.0
shipDatadict['wettedSurface'] = .794059157 *(shipDatadict['Lpp']*(2*shipDatadict['meanDraft']+shipDatadict['Beam']))
shipDatadict['underWaterLateralArea'] = shipDatadict['Lpp']*shipDatadict['meanDraft']*0.9843655                      
shipDatadict['blockCoefficient'] = 0.704661757
shipDatadict['CenterofGravity'] = np.array([-0.002290749, 0, -0.415058824]) * np.array([shipDatadict['Lpp'] ,shipDatadict['Beam'],shipDatadict['meanDraft']])
shipDatadict['verticalCenterOfBoyancy'] = 0.453647058824 * shipDatadict['meanDraft'] 
shipDatadict['GyrationArms'] = np.array([0.4, 0.25, 0.25]) * np.array([shipDatadict['Beam'],shipDatadict['Lpp'],shipDatadict['Lpp']])
shipDatadict['SeparationPoint'] = 45.0

class KtKqMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up any class-level fixtures, if necessary
        pass
    
    def setUp(self):
        self._KtKqObj = KtKq.KtKq(shipDatadict)

    @classmethod    
    def tearDownClass(cls):
        pass

    def test_openwaterprop(self):
      
        print(self._KtKqObj)
        betas_arr=[ 
                    -10.0,
                    0.,
                    5.2000 ,
                    10.3000 ,
                    15.3000 ,
                    20.0000 ,
                    24.5000 ,
                    29.6000 
                    
                                ]
        BETASWAG = 30.0
        
        CtArr=[]
        CqArr=[]
        for BETASWAG in betas_arr:
            FCTWAG, FCQWAG = self._KtKqObj.openwaterprop(BETASWAG)
            if BETASWAG >=0 and FCQWAG >= 0: # on in 1. Quadrant we use wageningen
                CtArr.append((BETASWAG,FCTWAG))
                CqArr.append((BETASWAG,FCQWAG))
            elif FCQWAG < 0:
                break
            else:
                continue
        CtArr = np.array(CtArr)
        CqArr = np.array(CqArr)
        #fig = go.Figure(data=go.Scatter(x=CtArr[:,0], y=CtArr[:,1], mode='lines+markers', name='Ct'))
        #fig.add_trace(go.Scatter(x=CqArr[:,0], y=CqArr[:,1], mode='lines+markers', name='Cq'))
        #fig.show()
        self.assertTrue(FCQWAG <0) # we know we are done with 1. quadrant when FCQWAG becomes negatives
        
        print(f"FCTWAG: {FCTWAG}")
        print(f"FCQWAG: {FCQWAG}")



class TestHullThumbs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Set up any class-level fixtures, if necessary
        pass

    def setUp(self):
        # Initialize the HullThumbs object before each test method
        self.hull_thumbs = HullThumbs(shipDatadict)

    def tearDown(self):
        # Clean up any resources after each test method
        pass

    @classmethod
    def tearDownClass(cls):
        # Clean up any class-level fixtures, if necessary
        pass
    
    def test_linear_damping_initialization(self):
        # Test the linearDamping method
        self.hull_thumbs.linearDamping()
        
        # Verify the damping tables are initialized correctly
       
        # ... (set the expected values for the rest of the table)
        
        np.testing.assert_allclose(self.hull_thumbs.getHeaveDamping()[0,1], 2.18409471e+08,rtol=10)
        print(f" Heave damping {self.hull_thumbs.getHeaveDamping()}")
        
        # Repeat for PitchDamping, RollDamping, SurgeDamping, and SwayDamping
        # ...
    def test_squat(self):
        squaFP = self.hull_thumbs.getSquatAP()
        print (f" SquatFP {squaFP}")
        np.testing.assert_allclose(squaFP[0,1],-1.22197368e-02,rtol=0.1)
        pass
    
    def test_addedMassTOH(self):
        XDCOR = self.hull_thumbs.addedMassTOH('XDCOR')
        print('***** ADDED MASS XDCOR(TOH)')
        print (f"{XDCOR}")
        np.testing.assert_allclose(XDCOR[4,1],2.024,rtol=0.1)
        pass
    
    def test_dmimix(self):
        a,b = self.hull_thumbs.dmimix()
        np.testing.assert_allclose(a['XUDOT'],5.872,rtol=0.001)
        #print(a)
        pass
    def test_pmmmot(self):
        ucar = 10 * 0.5144 # 10 knots
        betad = np.radians(45.)
        gamma = delta = heel = epsil = 0.0
        udim,vdim,rdim,qdim,pdim,ddim = self.hull_thumbs.pmmmot(ucar,betad,gamma,delta,heel,epsil)
        np.testing.assert_allclose(np.sqrt(udim**2 + vdim**2),5.144,rtol=0.001)
        pass
    def test_speedfactor(self):
        U=5.144*np.cos(np.radians(45.0))
        V=5.144*np.sin(np.radians(45.0))
        R = 0.0
        coefs = {'XVR':497.3150,'XVV':46.72300}
        factors =[]
        for coef in coefs.keys():
            factor = self.hull_thumbs.speedfactor(coef,U,V,R)
            factors.append(factor)
        pass
        np.testing.assert_allclose(factors[0],0.0,rtol=0.001) #Since R = 0, factor should be 0
        np.testing.assert_allclose(factors[1],V*V,rtol=0.001) #the XVV should be multiplied with V*V
    def test_defval(self):
        motion = 1
        icoty = 0 # Base model Beta
        absc1,absc2 = self.hull_thumbs.defval(1,0)
        assert absc2 is None, "absc2 should be none"
        assert len(absc1)==31,"absc1 should have length 31"
        pass
    def test_calculate_force_coefficients(self):
        motion = 2 # GAMMA table
        tableType = 'YAW'
        icoty = 0 # Base model Beta
        uo = self.hull_thumbs.serviceSpeed
        with open(r"H:\GitRepos\ShipYard2\data\PMMdata\3949NeuPMMDe_005_001.txt",'r') as f:
            pmmCoefs = {}
            for line in f:
                splt = line.split()
                pmmCoefs[splt[0]] = float(splt[1])/1.0E5           
        absc1, absc2 = self.hull_thumbs.defval(motion, icoty)
        force_coefficients, correction_table = self.hull_thumbs.calculate_force_coefficients(motion, tableType, absc1, absc2, icoty, uo, pmmCoefs)
        print ("test_calculate_force_coefficients MISSING assert")

    def test_yaw_drift_Tables(self):
        with open(r"H:\GitRepos\ShipYard2\data\PMMdata\3949NeuPMMDe_005_001.txt",'r') as f:
            pmmCoefs = {}
            for line in f:
                splt = line.split()
                pmmCoefs[splt[0]] = float(splt[1])/1.0E5   
        motion = 3
        icoty = 0                
        uo = self.hull_thumbs.serviceSpeed
        ipmms = 1
        tableType ='YAW_DRIFT'
        absc1,absc2 = self.hull_thumbs.defval(motion,icoty)
        tables,correctionTables = self.hull_thumbs.yaw_drift_Tables(motion, tableType, absc1, absc2, icoty, uo, pmmCoefs)
        dfsy1 = pd.read_csv(r"H:\GitRepos\ShipYard2\test\testData\SY13949YAW_DRIFT_X_HL.dat",header=0,sep='\s+')
        dfsy1.set_index('BETAD',inplace=True)        
        np.testing.assert_allclose(tables['X_HL'].loc[5.0][25.0]-dfsy1.loc[5.0]['25'],0.00,rtol=0.01,atol=1e-9)

        pass
    def test_getRollHeelTables(self):
        tables = self.hull_thumbs.getRollHeelTables()
        idum = 0
    def test_getForceTables(self):
        with open(r"H:\GitRepos\ShipYard2\data\PMMdata\3949NeuPMMDe_005_001.txt",'r') as f:
            pmmCoefs = {}
            for line in f:
                splt = line.split()
                pmmCoefs[splt[0]] = float(splt[1])/1.0E5       
        # test getting drift tables 
        baseTables,correctionTables = self.hull_thumbs.getForceTables(pmmCoefs,'DRIFT')
        np.testing.assert_allclose(baseTables['X_HL = DRIFT'].loc[-135.0],-0.0016682163766285373,rtol=0.01)
        # test getting YAW tables
        baseTables,correctionTables = self.hull_thumbs.getForceTables(pmmCoefs,'YAW')
        np.testing.assert_allclose(baseTables['X_HL = YAW'].loc[26.0],-7.491886E-04,rtol=0.001)
        idum = 0
        pass
    def test_resistance(self):
        shallow = False
        resistanceTable = self.hull_thumbs.RDHOLTROP(shallow)
        np.testing.assert_allclose(resistanceTable.iloc[0][0],0.0020987559933257744,rtol=0.001)
        pass
    def test_shallowWaterResistance(self):
        shallowCorrection = self.hull_thumbs.shallowWaterResistance()
        np.testing.assert_allclose(shallowCorrection.iloc[-1,-1],5.043079,rtol=0.001)
        pass
    


class TestPropellerSteering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up any class-level fixtures, if necessary
        pass

    def setUp(self):
        # Initialize the HullThumbs object before each test method
        x_values = np.array([-180, -165, -150, -135, -120, -110, -100, -90, -70, -60, -45, -30, -20, -10, 0, 5.2, 10.31, 15.26, 19.99, 25.86, 30, 45, 60, 70, 90, 100, 110, 120, 135, 150, 165, 180])
        y_values = np.array([-0.25, -0.1, 0.0875, 0.2841, 0.4684, 0.5021, 0.4973, 0.51, 0.5738, 0.6166, 0.5484, 0.3444, 0.2472, 0.19, 0.26449, 0.21944, 0.17034, 0.12139, 0.07175, 0, -0.0606, -0.2801, -0.5668, -0.4527, -0.4627, -0.4833, -0.5157, -0.5668, -0.5169, -0.4, -0.3, -0.25])

# Create a NumPy array with shape (n, 2) where n is the number of data points
        self.CTdata = np.column_stack((x_values, y_values))
        
    def test_propellerSteering(self):
        
        self.propellerSteer = propellerStering(shipDatadict, CTdata = self.CTdata)
        FY = self.propellerSteer.propellerSteeringForce(1)
        np.testing.assert_allclose(FY[0,0],0.025,rtol=0.01)
        print ('NOTE *** Propeller steering calculation only testet for FP prop')
        pass


    def tearDown(self):
        # Clean up any resources after each test method
        pass

    @classmethod
    def tearDownClass(cls):
        # Clean up any class-level fixtures, if necessary
        pass    
# %%


if __name__ == '__main__':
    unittest.main()

