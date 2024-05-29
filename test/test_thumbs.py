# %%
import unittest

import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import os
import plotly.graph_objects as go

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))


import KtKq
from HullThumbs import HullThumbs
from propellerSteering import propellerStering

# %%
class KtKqMethods(unittest.TestCase):
    KtKqObj = KtKq.KtKq(3949)
    @classmethod
    def setUpClass(cls):
        cls._KtKqObj = KtKq.KtKq(3949)

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
            FCTWAG, FCQWAG = self.KtKqObj.openwaterprop(BETASWAG)
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
        self.hull_thumbs = HullThumbs(shipnr=3949, propellerDiameter=5.0)

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
        
        self.propellerSteer = propellerStering(shipnr=3949, CTdata = self.CTdata)
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
# %%
