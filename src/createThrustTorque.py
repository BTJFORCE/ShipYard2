# %%

import numpy as np
import pandas as pd
from scipy import interpolate
from io import StringIO

from pathlib import Path


# Create an interpolation function based on the depth (z) values
def interpolate_matrices(a1, a2, interval, z):
    m, n = a1.shape
    a3 = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            a3[i, j] = interpolate.interp1d(interval, [a1[i, j], a2[i, j]])(z)
    return a3


class ThrustTorque:


    fileNames={}

    fileNames['fp80Thrust'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\fp80_Thrust.dat")

    fileNames['fp50Thrust'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\fp50_Thrust.dat")

    fileNames['fp80Torque'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\fp80_Torque.dat")

    fileNames['fp50Torque'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\fp50_Torque.dat")


    fileNames['cp50thrust'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp50_thrust.dat")

    fileNames['cp65thrust'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp65_thrust.dat")

    fileNames['cp80thrust'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp80_thrust.dat")

    fileNames['cp50torque'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp50_torque.dat")

    fileNames['cp65torque'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp65_torque.dat")

    fileNames['cp80torque'] = Path(r"H:\GitRepos\ShipYard2\data\propeller\cp80_torque.dat")

    cpAeA0Data=[0.50,0.65,0.80]


    def ThrustTorque(self,ThrustBool,Aratio,PD,propellerType):

        FP = True

        file1 = None

        file2 = None

        if propellerType.upper() == 'CP':

            FP = False

        if FP:   # we only have data for Aratio between 50 and 80

            if Aratio >= 0.8:

                if ThrustBool:

                    file1 = self.fileNames['fp80Thrust']

                else:

                    file1 = self.fileNames['fp80Torque']

            elif Aratio <= 0.5:

                if ThrustBool:

                    file1 = self.fileNames['fp50Thrust']

                else:

                    file1 = self.fileNames['fp50Torque']

            else:

                if ThrustBool:

                    file1 = self.fileNames['fp50Thrust']

                    file2 = self.fileNames['fp80Thrust']

                else:

                    file1 = self.fileNames['fp50Torque']

                    file2 = self.fileNames['fp80Torque']


            if not file2: # no need for aratio interpolation

                # perform PD interpolation

                df = pd.read_csv(file1, sep='\s+', header = 0)

                # find the columns to interpolate between

                cols = [float(val) for val in df.columns]

                index = np.searchsorted(cols, PD)

                xp= cols[index-1:index+1]

                df['Thrust'] = df.apply(lambda row: np.interp(PD,xp,row[index-1:index+1]), axis=1)

                return df['Thrust']

            else:

                print('we need to do AeAo interpolation')

                df50 = pd.read_csv(file1, sep='\s+', header = 0)

                df80 = pd.read_csv(file2, sep='\s+', header = 0)

                # we know the PD colums are the same so find the PD columns 

                # using df50 

                cols = [float(val) for val in df50.columns]

                index = np.searchsorted(cols, PD)

                xp= cols[index-1:index+1]

                # Now create the thrust curve for the 50% and 80% Arearatio

                df=df50.apply(lambda row: np.interp(PD,xp,row[index-1:index+1]), axis=1).to_frame(name='50')

                df['80']=df80.apply(lambda row: np.interp(PD,xp,row[index-1:index+1]), axis=1)

                # now interpolate to the actaual arearatio

                df['Thrust'] = df.apply(lambda row: np.interp(Aratio,[50,80],row[0:2]),axis =1)

                return df['Thrust']

        else: #" CP prop"
            # one or two files and which files to use
            index = np.searchsorted(self.cpAeA0Data, Aratio)
            if index == 3 :
                if ThrustBool:
                    file1 = self.fileNames['cp80thrust'] 
                else:
                    file1 = self.fileNames['cp80torque']
            elif index == 0:
                if ThrustBool:
                    file1 = self.fileNames['cp50thrust'] 
                else:
                    file1 = self.fileNames['cp50torque'] 
            elif index == 1:
                AeA0Range=[.50,.65]
                if ThrustBool:
                    file1 = self.fileNames['cp50thrust']
                    file2 = self.fileNames['cp65thrust']
                    
                else:
                    file1 = self.fileNames['cp50torque']
                    file2 = self.fileNames['cp65torque']
            elif index == 2:
                AeA0Range=[.65,.80]
                if ThrustBool:
                    file1 = self.fileNames['cp65thrust']
                    file2 = self.fileNames['cp80thrust']
                else:
                    file1 = self.fileNames['cp65torque']
                    file2 = self.fileNames['cp80torque']
            else:
                print('SOMETHING IS WRONG we only have 3 different CP data')

            if not file2: #No interpolation just return the df
                df = pd.read_csv(file1, sep='\s+', header = 0)
                return df
            else:
                df1 = pd.read_csv(file1, sep='\s+', header = 0)
                df2 = pd.read_csv(file2, sep='\s+', header = 0)
                a1 = df1.values
                a2 = df2.values
                a3=interpolate_matrices(a1,a2,AeA0Range,Aratio)
                df = pd.DataFrame(a3,columns=df1.columns)
                df.index=df1.index
                
                
       

                # Then, we interpolate the salinity at the target depth
                return df


    

    
    

# %%

if __name__ == '__main__':

    AeA0 = 0.98

    PD = 0.71505  # 0.9040 
    ThrustData = True

    aThrustTorqueData = ThrustTorque()
    df = aThrustTorqueData.ThrustTorque(ThrustData,AeA0,PD,'FP')
    #print(df)
    #df = aThrustTorqueData.ThrustTorque(ThrustData,AeA0,PD,'CP')
    df.head()
    pass

    
    
   




# %%



# %%

