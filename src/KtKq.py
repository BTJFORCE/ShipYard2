import math
import Thumbs as pThumbs

import numpy as np
import scipy.io


class KtKq(pThumbs.Thumbs):

    def __init__(self,shipnr,propellerDiameter=None):
        super().__init__(shipnr,propellerDiameter=None)



    def wageningen(self,Ja, PD, AEAO, z):
        """
        [KT, KQ] = wageningen(Ja,PD,AEAO,z) computes the thrust and torque
        coefficients of the Wageningen B-series propellers for:

        Ja = Va/(n*D)  advance number
        PD             pitch/diameter ratio (typically 0.5-2.5)
        AEAO           blade area ratio (ratio of expanded blade area to propeller disk area)
        z              number of propeller blades

        The B-series propellers were designed and tested at the Netherlands Ship
        Model Basin in Wageningen. The open_water characteristics of 120
        propeller models of the B_series were tested at N.S.M.B. and analyzed
        with multiple polynomial regression analysis. The derived polynomials
        express the thrust and torque coefficients in terms of the number of
        blades, the blade area ratio, the pitch_diameter ratio and the advance
        coefficient.

        Reference: Barnitsas, M.M., Ray, D. and Kinley, P. (1981).
        KT, KQ and Efficiency Curves for the Wageningen B-Series Propellers
        http://deepblue.lib.umich.edu/handle/2027.42/3557
        """
        data = scipy.io.loadmat(r"H:\GitRepos\ShipYard2\data\propeller\WageningData.mat")

        WagCThrust_stuv = data['WagCThrust_stuv']
        WagThrust_s = data['WagThrust_s']
        WagThrust_t = data['WagThrust_t']
        WagThrust_u = data['WagThrust_u']
        WagThrust_v = data['WagThrust_v']

        WagCTorque_stuv = data['WagCTorque_stuv']
        WagTorque_s = data['WagTorque_s']
        WagTorque_t = data['WagTorque_t']
        WagTorque_u = data['WagTorque_u']
        WagTorque_v = data['WagTorque_v']

        KT = np.sum(WagCThrust_stuv * ((Ja)**WagThrust_s) * (PD**WagThrust_t) * (AEAO**WagThrust_u) * (z**WagThrust_v))
        KQ = np.sum(WagCTorque_stuv * ((Ja)**WagTorque_s) * (PD**WagTorque_t) * (AEAO**WagTorque_u) * (z**WagTorque_v))

        return KT, KQ

        


    def openwaterprop(self, BETASWAG):
        
       

        Ja = (np.tan((BETASWAG*np.pi)/180.))*(0.7*np.pi)    

        FCTWAG,FCQWAG = self.wageningen(Ja,self.PD,self.pRatio,self.nrBlades)

        return FCTWAG, FCQWAG
    
    def CtCq1Qaudrant(self,betas_arr):
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
        return CtArr,CqArr




    def Ktkq(self,ajinp, X, Y, Z, Diam, N5, Tcs, Tc, Rnx, K1, K3, K1noz, V5, V6, Ic, Jc,L9,Txbase):

        assert False, "Warning: 'old_method' is deprecated and just a copy from FORTAN t"
        """
        This function calculates Kt and Kq for given J.
        Uncorrected Kt and Kq are calculated by call of Ktkqx.
        These values are then corrected, if necessary and wanted, for:
            emergence
            partial cavitation
                                                    23. March 1994
        ajinp	Ja avanceringstal (number of advance)
        X		Ae/A0 area ratio
        Y		P/D at 0.7*R 
        Z		Number of propeller blades
        Diam	Propeller diameter
        N5		Rouns per second of propeller revolutions
        Tcs		Standard thickness-chord length ratio at 0.75*R
        Tc		Userdefined thickness-chord length ratio at 0.75*R
        Rnx		Something to do with Tc or Tcs??????
        K1		kt
        K3		kq
        K1noz	nozzle kt???????
        V5		Ship speed
        V6		Va water inflow speed (wake corrected)
        Ic		emergence parameter
        Jc		cavitation parameter

                Other input prametersers via common areas!!!!!
        L9		is ship length, required as input in common for some calculations


        """
        # Initialize variables
        C8 = self.waterLineBlock
        B9 = self.beam
        aj = ajinp if ajinp > 1e-10 else 1e-10
        Ic = 0
        Jc = 0
        Aji = aj
        Ans = 46.4301 * Y**(-1.746) * (10-Z)**(-2.223)
        Aps = 15.1845 * Y**(-2.2514) * (10-Z)**(-1.4478)
        Urx = .09
        Txaft = Txbase
        Fns = V5 / math.sqrt(L9 * 9.8066)
        if Fns < .3:
            Urx = Fns**2
        Aux = Txbase + Diam * .5 - Txaft - .6 * C8 * B9 * Urx
        Urele = 0.0

        Ktbi = .06218 + .1194 * Y - .00249 * Z
        self.ktkqx(Aji, X, Y, Z, Diam, N5, Tcs, Tc, Rnx, K1, K3, K1noz)
        if Isubco == 0:
            if Urele > 1e-3:
                Ic = 1
                Gf = K1 / Aji / Aji
                Aux = (1 + .2 / .4) / 3 / Urele - .4 / Aji
                if Gf > Aux:
                    K1 = .00001
                    K3 = .00001
                    Ic = -1
                    return
                Gf = 1. + 3. * Urele * Gf
                Aj1 = Aji * Gf
                self.ktkqx(Aj1, X, Y, Z, Diam, N5, Tcs, Tc, Rnx, K1, K3, K1noz)
        if Icavco == 1:
            Ktac = K1 / Aji / Aji
            if Ktac > 200.:
                Ktac = 200.
            Ktac = Ktac / Sigma0
            Fns = 1.
            Fps = 1.
            Jc = 0
            Aux = Ktac - Ktbi
            if Aux > 0.:
                Jc = 1
                Fns = 1 + Ans * Aux**1.2
            Ajx = Aji / Fns
            Auy = Aj - Ajx
            if abs(Auy) > 1e-6:
                Aji = Auy * .5 + Aji
                return Ktkq(ajinp, X, Y, Z, Diam, N5, Tcs, Tc, Rnx, K1, K3, K1noz, V5, V6, Ic, Jc)
            Aux = Aux - .01
            if Aux > 0.:
                Fps = 1. + Aps * Aux**1.2
            K1 = K1 / Fns / Fns
            K3 = K3 * Fps / Fns**3
