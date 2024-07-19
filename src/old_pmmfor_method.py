def pmmfor(self,pmmCoefs,coftyp,uo,udim,vdim,rdim,qdim,pdim,ddim,toh):
    '''
    this function calculates the forces on the hull using the PMM coefficient/hydrodynamic derivatives
    the function is translate to python by BTJ originally part of shipYard1 see tfs
    $/SimFlex Classic/src/lib/core/msdat/lib/pmm/pmmfor.f
    '''
    acoef = self.acoef
    bcoef = self.bcoef
    fdim=np.array([0.,0.,0.])
    fdimu2=np.array([0.,0.,0.])
    ## SKIPPING D,Q,P derivatives as we dont have them
    xderivatives=['X0','XU','XUU','XUUU','XV','XVV','XR','XRR','XVR','XUDOT','XRDOT','XUIUI']
    yderivatives=['Y0','Y0U','YUU','YV','YVV','YVIVI','YVVV','YVU','YR','YRR','YRIRI','YRRR','YRIVI','YRVV','YVIRI','YVRR','YVDOT','YRDOT']
    nderivatives=['N0','N0U','NUU','NV','NVV','NVIVI','NVVV','NVU','NR','NRR','NRIRI','NRRR','NRIVI','NRVV','NVIRI','NVRR','NVDOT','NRDOT ']
    lpp = self.Lpp
    
    if coftyp == 1:
      utot2 = udim*udim + vdim*vdim
      utot  = np.sqrt(utot2)
      U   = (udim - uo)/utot
      V   = vdim/utot
      R   = rdim/(utot/lpp)
      Q  = qdim
      P   = pdim/(utot/lpp)
      D   = ddim
      usign = 1.0
    else:
      utot2 = 1.
      utot  = 1.0
      U   = udim
      V   = vdim
      R   = rdim*lpp
      Q   = qdim
      P   = pdim*lpp
      D  = ddim
#     if vdim and udim is zero we jump over the next because when gamma is -90 or 90
#     this would lead to not a number (NAN) for BETAMSMFQ, which are used in the XVVSPEC coef.
      if (np.sqrt(vdim**2+udim**2) < 0.000001) :
        usign = np.sign(U)
        '''
        if U == 0:
          usign = 0
        elif U > 0:
          usign = 1
        elif U < 0:
          usign = -1
        '''
      else:
        BETAMSMFQ=np.sign(udim)*(-np.arcsin(vdim/np.sqrt(vdim**2+udim**2)))
        usign = np.sign(U)
        pass
    # X --- force
    xcoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    #inserting coefs and speed in arrays so we can do simple array multiplication to calculate force
    # check for key in dict needs to be added
    for val in xderivatives:
      if val in pmmCoefs.keys():
        if 'DOT' in val:
          continue  # in the fortran PMMFOR no DOT derivatives appears
        xcoeff=np.append(xcoeff,pmmCoefs[val])
        speedvector=np.append(speedvector,self.speedfactor(val,U,V,R))
        if val in acoef.keys():
          a = acoef[val]
          b = bcoef[val]
        else:
          a = 0
          b = 1
        TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))
        
    fdimu2[0]= np.sum(xcoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**2
    fdim[0] = fdimu2[0]*utot2
    # Y-force
    ycoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    for val in yderivatives:
      if val in pmmCoefs.keys():
        if 'DOT' in val:
          continue  # in the fortran PMMFOR no DOT derivatives appears
        ycoeff=np.append(ycoeff,pmmCoefs[val])
        speedvector=np.append(speedvector,self.speedfactor(val,U,V,R)) 
        if val in acoef.keys():
          a = acoef[val]
          b = bcoef[val]
        else:
          a = 0
          b = 1
        TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b))       
    
    fdimu2[1]= np.sum(ycoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**2
    fdim[1] = fdimu2[1]*utot2    
    # N-moment
    ncoeff=np.array([])
    speedvector=np.array([])
    TOHCorrection=np.array([])
    for val in nderivatives:
      if val in pmmCoefs.keys():
        if 'DOT' in val:
          continue  # in the fortran PMMFOR not DOT derivatives appears
        ncoeff=np.append(ncoeff,pmmCoefs[val])
        speedfactor = self.speedfactor(val,U,V,R)
        speedvector=np.append(speedvector,speedfactor) 
        # comment from fortran code
        #It is believed that N must change its sign for ship going astern
        #Therefore, USIGN is multiplied to the NvIvI
        if 'VIVI' in val:
          speedfactor *= usign
        if val in acoef.keys():
          a = acoef[val]
          b = bcoef[val]
        else:
          a = 0
          b = 1            
        TOHCorrection=np.append(TOHCorrection,(1+ a * toh**b)*speedfactor)        
    
    fdimu2[2]= np.sum(ncoeff * TOHCorrection * speedvector) * 0.5 * self.rho * lpp**3
    fdim[2] = fdimu2[2]*utot2      
    
    return fdim,fdimu2,utot2  