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
shipDatadict['blockCoefficient'] = 0.704661757
shipDatadict['meanDraft'] = (shipDatadict['draftAft'] + shipDatadict['draftFore'] )/2.0
shipDatadict['wettedSurface'] = 0.7801584*(shipDatadict['Lpp'] *(2*shipDatadict['meanDraft']+shipDatadict['Beam']))
shipDatadict['underWaterLateralArea'] = 0.9124 * shipDatadict['Lpp'] * shipDatadict['meanDraft']

shipDatadict['CenterofGravity'] = np.array([-0.002290749, 0, -0.415058824]) * np.array([shipDatadict['Lpp'] ,shipDatadict['Beam'],shipDatadict['meanDraft']])
shipDatadict['verticalCenterOfBoyancy'] = 0.453647058824 * shipDatadict['meanDraft'] 
shipDatadict['GyrationArms'] = np.array([0.4, 0.25, 0.25]) * np.array([shipDatadict['Beam'],shipDatadict['Lpp'],shipDatadict['Lpp']])
shipDatadict['SeparationPoint'] = 45.0

class KtKqMethods(unittest.TestCase):

    KtKqObj = KtKq.KtKq(shipDatadict)
    @classmethod
    def setUpClass(cls):
        cls._KtKqObj = KtKq.KtKq(shipDatadict)

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
        self.hull_thumbs = HullThumbs(shipDatadict)

    def test_linear_damping_initialization(self):
        # Test the linearDamping method
        self.hull_thumbs.linearDamping()
        
        # Verify the damping tables are initialized correctly
       
        # ... (set the expected values for the rest of the table)
        
        np.testing.assert_allclose(self.hull_thumbs.getHeaveDamping()[0,1], 2.18409471e+08,rtol=10)
        print(f" Heave damping {self.hull_thumbs.getHeaveDamping()}")
        
        # Repeat for PitchDamping, RollDamping, SurgeDamping, and SwayDamping
        # ...

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