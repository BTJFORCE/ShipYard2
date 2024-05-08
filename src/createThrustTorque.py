# %%
import numpy as np
import pandas as pd
from io import StringIO
from pathlib import Path

def thrust(Aratio,PD):

    df = pd.read_csv(Path("U:/ships/lib/propeller/fp50_Thrust.dat"), delim_whitespace=True, header = 0)
    return df
# %%
if __name__ == '__main__':
    df = thrust(0.98,77)



# %%
