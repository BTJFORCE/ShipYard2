import math
import Thumbs as pThumbs

import numpy as np
import scipy.io



class propellerStering(pThumbs.Thumbs):
    
    def __init__(self,shipnr,CTdata,propellerDiameter=None):
        super().__init__(shipnr,propellerDiameter=None)    
        
        self.CTdata = CTdata
        self.numBETAs = len(CTdata)
        self.numPDs = np.shape(np.array(self.CTdata))[1]-1
        
    def findZeroCrossings(self,CTdata):
        x = CTdata[:,0]
        y = CTdata[:,1]
        # Find the indices where y changes sign
        zero_crossings_indices = np.where(np.diff(np.sign(y)))[0]
        # The actual zero crossings are between these indices in the x array
        zero_crossings_x = x[zero_crossings_indices]
        return zero_crossings

    def propellerSteeringForce(self,propellerDirection):
        betas_arr= self.CTdata[:,0]
        fsway = np.zeros((self.numBETAs,self.numPDs))
     
        for ixPD in range(len(self.numPDs)):
            for ixBETA in range(len(self.numBETAs)):
                FCOR = self.propStee(betas_arr[ixBETA],self.CTdata[:,ixPD+1])
                fsway[ixBETA,ixPD] = self.CTdata[ixBETA,ixPD]*FCOR
                if self.propellerType =='FP':
                    fsway[ixBETA,ixPD] = np.abs(fsway[ixBETA,ixPD] )
                fsway[ixBETA,ixPD] *= propellerDirection
        return fsway
            
        

    def propStee(self,xeval,PD):
        '''
        rewrite lib/core/msdat/lib/thmb/propstee.f by BTJ
        '''
        
        if self.propellerType == 'FP':
            xvals =  [-180.0,self.zero_crossings[0],-90.0,0.0,self.zero_crossings[1] ,90.0,135.0,180.0]
            yvals = [0.1,0.0,0.0,0.02,0.0,0.0,0.01,0.1]
        elif self.propellerType == 'CP' and PD > 0:
            xvals = [-180.0,-90.0,-45.0,0.0,self.zero_crossings[0],45.0,90.0,180.0]
            yvals = [0.0,0.0,0.03,0.02,0.0,0.01,0.0,0.0]
        elif self.propellerType == 'CP' and PD == 0:
            xvals = [-180.0,-90.0,-45.0,self.zero_crossings[0],45.0,90.0,180.0]
            yvals = [0.0,0.0,0.01,0.00,0.01,0.0,0.0]
        elif self.propellerType == 'CP' and PD < 0:
            xvals = [-180.0,-90.0,self.zero_crossings[0]/2.0,self.zero_crossings[0],0.0,45.0,90.0,180.0]
            yvals = [0.0,0.0,0.01,0.00,0.1,0.01,0.0,0.0]
            if self.zero_crossings[0] == 0:
                xvals[2] = -45
        else:
            assert('some thing is wrong should never end here !!!')
        ## calculate the correction factor
        
        for ix in range(len(xvals)) :
            if xeval < xvals[ix]:
                slope = yvals[ix] - yvals[ix-1]
                slope = slope / (xvals[xi] - xvals[ix-1])
                yfunc = (xeval-xvals[ix-1])*slope+yvals[ix-1]
                break
        return yfunc


    
        
        
    