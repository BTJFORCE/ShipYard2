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

# %%
class KtKqMethods(unittest.TestCase):

    def test_openwaterprop(self):

        KtKqObj = KtKq.KtKq(3949)

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
            FCTWAG, FCQWAG = KtKqObj.openwaterprop(BETASWAG)
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


# %%


if __name__ == '__main__':
    unittest.main()