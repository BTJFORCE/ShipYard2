
import Thumbs as Thumbs
class HullThumbs(Thumbs.Thumbs):
    
    def __init__(self,shipnr,propellerDiameter=None):
        super().__init__(shipnr,propellerDiameter=None)

    def linearDamping(self):
        	pi=3.14159265359
      raddeg = 180./pi
	RHO=1026.0  
	grav = 9.81

C Below, the necessary inputs to the damping method are collected. 

	ires = 0
      irds = 0

      call getvanc('<LPP>',2,IDUM,L5,ires)
      call getvanc('<LPP>',2,IDUM,L5,ires)
c	write(6,*)'lpp ',L5
	call getvanc('<BEAM>',2,IDUM,B9,ires)
C --- Dimensionalize
	B9 = B9*L5
c	write(6,*)'beam ',B9
	call getvanc('<BDRAU(2)>',2,IDUM,TXAFT,ires)
	TXAFT = TXAFT*L5
c	write(6,*)'Draught aft ',TXAFT
	call getvanc('<BDRAU>',2,IDUM,TXFORE,ires)
	TXFORE = TXFORE*L5
c	write(6,*)'Draught fore ',TXFORE
	FMFORE = TXFORE
	FMAFT = TXAFT
	T9 = .5*(TXAFT+TXFORE)  ! Mean Draught
	call getvanc('<BBLCK>'    ,2,IDUM,V0,ires)
	C2 = V0   ! Block coefficient (Lpp)
	V0 = V0*(L5*B9*T9)    !Displacement volume
	call getvanc('<BCGG>'    ,2,IDUM,COG,ires)
	COG(1) = COG(1)*L5
	COG(3) = COG(3)*T9
	write(6,*)'COG meters, ',COG(1),COG(3)
C	calculate water plane area
	call getvanc('<BALW>'    ,2,IDUM,ALW,ires)
	ALW=ALW*L5*T9   ! water plane area
	call getvanc('<BROG>'    ,2,IDUM,GYRATION,ires)
	GYRATION(1) = GYRATION(1)*B9
	GYRATION(2) = GYRATION(2)*L5
	write(6,*)'Radii of gyration ',GYRATION(1), GYRATION(2)
	call getvanc('<BCBG>'    ,2,IDUM,VCOB,ires)
	VCOB =VCOB*T9
	IXX = GYRATION(1)**2*V0*RHO
	IYY = GYRATION(2)**2*V0*RHO
C	call getvanc('<BAFW>'    ,2,IDUM,AFRONT,ires)
C Now input additional DEN-Mark1 variables
	IRDS=0
	CALL RDS('<Enter Transverse BM in m>','>',2,1,IDUM,BMT,IRDS)
	IRDS=0
	CALL RDS('<Enter Longitudinal BM in m>','>',2,1,IDUM,BML,IRDS)
C	
C	calculate the GMT and GML: GM = KB + BM - KG
C
      IGMTNEG = 0
	GMT = (T9-VCOB) + BMT - (T9 - COG(3))
      IF ( GMT .LT. 0.0 ) THEN
	  WRITE(6,*) 'Negative GM, the ship will capsize'
	  IGMTNEG = 1
	ENDIF
	GML = (T9-VCOB) + BML - (T9 - COG(3))
C
c	calculate the resonance frequencies
C
	IF (IGMTNEG .EQ. 0 ) THEN
	  RDUMMY = grav*V0*rho*GMT/IXX
	  omerol = SQRT(RDUMMY)
	ELSE
	  omerol = 0.000001
      ENDIF
	RDUMMY = grav*V0*rho*GML/IYY
	IF (RDUMMY .GT. 0.0 ) THEN
	  omepit = SQRT(RDUMMY)
	ELSE
	  omepit = 0.000001
      ENDIF
	RDUMMY = ALW*grav/V0
	IF (RDUMMY .GT. 0.0 ) THEN
	  omehea = SQRT(RDUMMY)
	ELSE
	  omehea = 0.000001
	ENDIF
	write(*,*)'resonans periods, roll , pitch, heave ',
     &             1/omerol*2*PI, 1/omepit*2*PI, 1/omehea*2*PI

C	Very rough assumptions are applied here. one section shape, 
C	pitch equal to roll
C	calculate the Ap coefficients (dp assumed = 0.5)
	DPROLL = 0.5
C	DPPITCH = 0.05
	DPPITCH = 0.1*B9/L5
	APROLL = DPROLL*(omerol**2*B9/(2.0*grav))**2
	APPITC = DPPITCH*(omepit**2*L5/(2.0*grav))**2
	APHEAV = 2*SIN(omehea**3*B9/(2*grav))*EXP(-omehea**2*T9/grav)

	BSECROLL = RHO*grav**2/omerol**3*(B9/2)**2*APROLL**2
	ROLLDAMP = BSECROLL*L5

	BSECPITC = RHO*grav**2/omepit**3*(L5/2)**2*APPITC**2
	PITCDAMP = BSECPITC*B9

	BSECHEAV = RHO*grav**2/omehea**3*APHEAV**2
	HEAVDAMP = BSECHEAV*L5
	WRITE(6,*)'roll , pitch, heave, ',ROLLDAMP, PITCDAMP, HEAVDAMP
C	write(*,*)fcombn
C	NOW MAKE THE TABLES
C------------------------------------------------------------------------
C	HEAVE TABLES
C------------------------------------------------------------------------
      NDEF1 = 3
	ABSC1(1) = -10.0
	ABSC1(2) =   0.0
	ABSC1(3) =  10.0
	YVAL(1) = ABSC1(3)*HEAVDAMP
	YVAL(2) = 0.0
	YVAL(3) = -ABSC1(3)*HEAVDAMP
	
	DO 104 I=1,3
	  TABSC1bla(i)=0
 104	CONTINUE

C --- Generate table name
        
	  NOTAB = 1
C
        CALL coptxti('<LINEAR-HEAVE>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<Z_HL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<W_HL>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

        CALL coptxti('<HF-HEAVE-DAMPING>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<ZGHL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<DWHBW>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

C------------------------------------------------------------------------
C	PITCH TABLES
C------------------------------------------------------------------------
      NDEF1 = 3
	ABSC1(1) = -1.0
	ABSC1(2) =   0.0
	ABSC1(3) =  1.0
	YVAL(1) = PITCDAMP
	YVAL(2) = 0.0
	YVAL(3) = -PITCDAMP
	
C --- Generate table name
        
	  NOTAB = 1
C
        CALL coptxti('<LINEAR-PITCH>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<M_HL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<Q_HL>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

        CALL coptxti('<HF-PITCH-DAMPING>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<MGHL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<DQHBW>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

C------------------------------------------------------------------------
C	ROLL TABLES
C------------------------------------------------------------------------
      NDEF1 = 3
	ABSC1(1) = -1.0
	ABSC1(2) =   0.0
	ABSC1(3) =  1.0
	YVAL(1) = ROLLDAMP
	YVAL(2) = 0.0
	YVAL(3) = -ROLLDAMP
	
C --- Generate table name
        
	  NOTAB = 1
C
        CALL coptxti('<LINEAR-ROLL>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<K_HL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<P_HL>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

        CALL coptxti('<HF-ROLL-DAMPING>',TABNAM)
        CALL coptxti('<     >',TABUNA)
        CALL coptxti('<KGHL>',CONNAM(1,NOTAB))
        CALL coptxti('<     >',CONUNA(1,NOTAB))
        CALL coptxti('<DPHBW>',INDNAM(1,1))
        CALL coptxti('<>',INDUNA(1,1))

        NVAR = 1
        IBTCM = 1
        ICTCM = 0
        J0 = 1
        IDUM = 1
        FAC = 1.0
        CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &              YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	  CALL UPDLM8

C     ----------------------------------------------------------------
C	Now calculate the high frequency dampings in sway, yaw and surge
C     ----------------------------------------------------------------
      NDEF1 = 3
	ABSC1(1) = -10.0
	ABSC1(2) =   0.0
	ABSC1(3) =  10.0
	YVAL(1) = 20.0
	YVAL(2) = 0.0
	YVAL(3) = -20.0

      CALL coptxti('<HF-SURGE-DAMPING>',TABNAM)
      CALL coptxti('<     >',TABUNA)
      CALL coptxti('<XGHL>',CONNAM(1,NOTAB))
      CALL coptxti('<     >',CONUNA(1,NOTAB))
      CALL coptxti('<DUHBW>',INDNAM(1,1))
      CALL coptxti('<>',INDUNA(1,1))

      NVAR = 1
      IBTCM = 1
      ICTCM = 0
      J0 = 1
      IDUM = 1
      FAC = 1.0
      CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &            YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	CALL UPDLM8

      CALL coptxti('<HF-SWAY-DAMPING>',TABNAM)
      CALL coptxti('<     >',TABUNA)
      CALL coptxti('<YGHL>',CONNAM(1,NOTAB))
      CALL coptxti('<     >',CONUNA(1,NOTAB))
      CALL coptxti('<DVHBW>',INDNAM(1,1))
      CALL coptxti('<>',INDUNA(1,1))

      NVAR = 1
      IBTCM = 1
      ICTCM = 0
      J0 = 1
      IDUM = 1
      FAC = 1.0
      CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &            YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	CALL UPDLM8

      NDEF1 = 3
	ABSC1(1) = -10.0
	ABSC1(2) =   0.0
	ABSC1(3) =  10.0
	YVAL(1) = 10.0
	YVAL(2) = 0.0
	YVAL(3) = -10.0

      CALL coptxti('<HF-YAW-DAMPING>',TABNAM)
      CALL coptxti('<     >',TABUNA)
      CALL coptxti('<NGHL>',CONNAM(1,NOTAB))
      CALL coptxti('<     >',CONUNA(1,NOTAB))
      CALL coptxti('<DRHBW>',INDNAM(1,1))
      CALL coptxti('<>',INDUNA(1,1))

      NVAR = 1
      IBTCM = 1
      ICTCM = 0
      J0 = 1
      IDUM = 1
      FAC = 1.0
      CALL MKTABL(IBTCM,ICTCM,J0,NDEF1,IDUM,absc1,RDUM,tabsc1bla,IDUM,
     &            YVAL,FAC)
C --- UPDLM8 updates the lump with the modified table
	CALL UPDLM8

C	Update the hydrstatic common area for later printout of hydrostatic properties
      
	KEELGRAV = (T9 - COG(3))
	KEELBUOY = (T9-VCOB)
	BUOMETTR = BMT
	BUOMETLO = BML
	METCENTR = GMT
	METCENLO = GML
	INERIAXX = IXX
	INERIAYY = IYY
	PERIODRO = 1/omerol*2*PI 
	PERIODPI = 1/omepit*2*PI
      PERIODHE = 1/omehea*2*PI
