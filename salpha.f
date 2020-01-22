      SUBROUTINE SALPHA(istout,ison,vdat,sdat,adat,edat,tcof,ntot,
     .                  pres,temp,NGAS2IMOL,sfin,ALPHAL,ALPHAD1)

********************************************************************
* PURPOSE:
*    In this subroutine the line intensity sfin is calculated at 
*    temperature temp, and the air-broadened halfwidth ALPHAL is 
*    calculated at temperature temp, and pressure pres. Also the
*    Doppler half-width is calculated and stored without the wave-
*    number in array ALPHAD1.
*
* temp0 : reference temperature of the HITR92 data base
*         temp0= 296 K
* pres0 : reference pressure of the HITR92 data base
*         pres0= 1 atm
* C2    : second radiation constant = hc/k
*         h : Planck's constant (6.6262E-34 J s)
*         c : speed of light (2.9979E10 cm s^-1)
*         k : Boltzmann's constant (1.3807E-23 J K^-1)
*         C2= 1.4388E0 cm K
* C3    : constant = c^{-1} (2k/mu)^{0.5}   (dimensionless)
*         mu : mass of 1 atomic mass unit (1.66E-27 kg)
* LOSCH : Loschmidt's number, the number of molecules per unit
*         volume at STP
*         LOSCH= 2.6868E19 cm^-3
*         
* PARAMETER: 
*   IWOUT  (1) write to output file kd.out
*          (0) don't
*
* DATE:
*    June , 1995
*
* Added: SFB = alpha_L(self)/alpha_L(foreign)
*
* AUTHOR:
*    D.M.Stam
********************************************************************
      IMPLICIT REAL*8 (a-h,o-z)

      INCLUDE 'max.incl'

      INTEGER istout,IWOUT,ntot,NGAS2IMOL,NMOL,NSPEC
      PARAMETER (IWOUT=0,NMOL=32,NSPECI=75)

      INTEGER ison(nvMAX)

      REAL*8  pres,temp,temp0,pres0,C2,C3,SFB,xch4,AD,xo2
      REAL*8  vdat(nvMAX),sdat(nvMAX),adat(nvMAX),edat(nvMAX),
     .        tcof(nvMAX),sfin(nvMAX),QRAT(nvMAX),ALPHAL(nvMAX),
     .        ALPHAD1(nvMAX)

      PARAMETER (temp0=296.0,pres0=1.0,C2=1.4388,C3=4.3022E-7,
     .           SFB=1.1,xch4=0.0022,xo2=0.209)

      COMMON / ISVEC / NISO(NMOL), MISO(NSPECI), isonm(NMOL)

Cf2py intent(in) istout,ison,vdat,sdat,adat,edat,tcof,ntot,pres,temp,NGAS2IMOL
Cf2py intent(out) sfin,ALPHAL,ALPHAD1
      
*-------------------------------------------------------------------
*     Initialize arrays, and calculate the constant TS1:
*-------------------------------------------------------------------
      DO i=1,ntot
         sfin(i)= 0.D0
         QRAT(i)= 0.D0
         ALPHAL(i)= 0.D0
         ALPHAD1(i)= 0.D0
      ENDDO

      TS1= C2*(temp-temp0)/(temp*temp0)

*-------------------------------------------------------------------
*     Loop over the wavenumbers:
*-------------------------------------------------------------------
      DO i=1,ntot
         TS2= (1.0-DEXP(-C2*vdat(i)/temp)) /
     .        (1.0-DEXP(-C2*vdat(i)/temp0))

*-------------------------------------------------------------------
*     Call the function QOFT to calculate the total internal
*     partition sum:
*-------------------------------------------------------------------
         CALL QOFT(NGAS2IMOL,ison(i),temp0,QTEMP0)
         CALL QOFT(NGAS2IMOL,ison(i),temp, QTEMP)
         QRAT(i)= QTEMP0 / QTEMP

*-------------------------------------------------------------------
*     Calculate the line strength (in cm/molecule):
*     (to obtain the line strength as used in Lacis & Oinas:
*     multiply with LOSCH (Loschmidt's number))
*-------------------------------------------------------------------
         sfin(i)= sdat(i) * QRAT(i) * TS2 * DEXP(edat(i)*TS1)

*-------------------------------------------------------------------
*     Calculate the Lorentz line width ALPHAL, and the Doppler
*     line width ALPHAD1 (without the wavenumber!):
*-------------------------------------------------------------------
*         AD= adat(i)*(xch4 + (1.0-xch4)/SFB)
         AD= adat(i)*(xo2 + (1.0-xo2)/SFB)
         ALPHAL(i)= AD * (pres/pres0) * (temp0/temp)**tcof(i)
         Ip= isonm(NGAS2IMOL) + ison(i)
         Mp= MISO(Ip)
         ALPHAD1(i)= C3 * (temp/DFLOAT(Mp))**0.5

      ENDDO

*-------------------------------------------------------------------
*     In case IWOUT=1, write the data to standard outputfile:
*-------------------------------------------------------------------
      IF (IWOUT.EQ.1) THEN
         WRITE(istout,*)
         WRITE(istout,300)
         WRITE(istout,305)
         DO 310 i=1,ntot
            WRITE(istout,315) ison(i),vdat(i),QRAT(i),sdat(i),
     .                   sfin(i),adat(i),ALPHAL(i),ALPHAD1(i)
310      CONTINUE
      ENDIF

300   FORMAT('SALPHA OUTPUT')
305   FORMAT('ison   vdat           QRAT    sdat        sfin',
     .       '       adat        ALPHAL      ALPHAD1')
315   FORMAT(I3,2X,F10.3,2X,F8.5,2X,E10.5,2X,E10.5,2X,E10.5,
     .       2X,E10.5,2X,E10.5)

********************************************************************
      RETURN
      END


********************************************************************
      BLOCK DATA ISOTOPES

*-------------------------------------------------------------------
* This data block contains the common block ISVEC
* ISVEC contains:
*     NISO(NMOL) : the number of isotopes of each molecule specie
*     MISO(NISO) : the molecular weight (in amu) of each isotope,
*                  arranged according to molecule specie
*     isonm(NMOL): the total number of preceding isotopes of array
*                  NISO, to find the place in array WISO having a 
*                  molecule number and isotope number
*-------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER NMOL,NISO
      PARAMETER (NMOL=32,NSPECI=75)

      COMMON / ISVEC / NISO(NMOL), MISO(NSPECI), isonm(NMOL)

* The number of isotopes for a particular molecule specie:
      DATA (NISO(i), i=1,NMOL) /
*     H2O,    CO2,     O3,    N2O,     CO,    CH4,     O2,
     +  4,	8,	3,	5,	5,	3,	3,
*      NO,    SO2,    NO2,    NH3,   HNO3,     OH,     HF,
     +  3,	2,	1,	2,	1,	3,	1,
*     HCL,    HBR,     HI,    CLO,    OCS,   H2CO,   HOCL,    
     +  2,	2,	1,	2,	4,	3,	2,
*      N2,    HCN,  CH3CL,   H2O2,   C2H2,   C2H6,    PH3,
     +  1,	3,	2,	1,	2,	1,	1,
*    COF2,    SF6,    H2S,  HCOOH
     +  1,	1,	1,	1/

* The molecular weights (in amu) of each molecule and its isotopes:
      DATA (MISO(i), i=1,NSPECI) /
*      H2O,
     + 	18,	20,	19,	19,	 
*      CO2,
     +	44,	45,	46,	45,	47,	46,	48,	47,
*       O3,                    N2O,
     +  48,	50,	50,	44,     45,     45,     46,     45,
*       CO,                                    CH4,
     +  28,	29,	30,	29,	31,	16,	17,	17,
*       O2,                     NO,
     +	32,	34,	35,	30,	31,	32,	
*      SO2,            NO2,    NH3,           HNO3,
     +  64,	66,	46,	17,	18,	63,
*       OH,                     HF,    HCL,            HBR,
     +  17,	19,	18,	20,	36,	38,	80,	81,
*       HI,    CLO,            OCS,
     + 128,	51,	53,	60,	62,	61,	62,
*     H2CO,                   HOCL,             N2,
     +  30,	31,	32,	52,	54,	28,
*      HCN,                  CH3CL,           H2O2,   C2H2,
     +  27,	28,	28,	50,	52,	34,	26,	27,
*     C2H6,    PH3,   COF2,     SF6,   H2S,  HCOOH
     +  30,	34,	66,	146,	34,	46/

* The total number of preceding isotopes per molecule specie:
      DATA (isonm(i), i=1,NMOL) /
*     H2O,    CO2,     O3,    N2O,     CO,    CH4,     O2,
     +  0,      4,     12,     15,     20,     25,     28,
*      NO,    SO2,    NO2,    NH3,   HNO3,     OH,     HF,
     + 31,     34,     36,     37,     39,     40,     43,
*     HCL,    HBR,     HI,    CLO,    OCS,   H2CO,   HOCL, 
     + 44,     46,     48,     49,     51,     55,     58,
*      N2,    HCN,  CH3CL,   H2O2,   C2H2,   C2H6,    PH3,
     + 60,     61,     64,     66,     67,     69,     70,
*    COF2,    SF6,    H2S,  HCOOH
     + 71,     72,     73,     74/

      END
