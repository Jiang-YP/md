C     User customized AspenPlus unit model to simulate the membrane
C     distillation in a hollow fiber membrane module
C
C     Developed by Guan Guoqiang @ SCUT, Guangzhou, China
C
C     Revision 2.1(16022), code alias: boiler
C
C     User Unit Operation Model (or Report) Subroutine for USER2
C
      SUBROUTINE DCMD  (NMATI,  SIN0,    NINFI,   SINFI,  NMATO,
     2                   SOUT,   NINFO,  SINFO,   IDSMI,  IDSII,
     3                   IDSMO,  IDSIO,  NTOT,    NSUBS,  IDXSUB,
     4                   ITYPE,  NINT,   INT,     NREAL,  REAL,
     5                   IDS,    NPO,    NBOPST,  NIWORK, IWORK,
     6                   NWORK,  WORK,   NSIZE,   SIZE,   INTSIZ,
     7                   LD   )
C
C     Invoke the common module in commod.f to define global variables
      use CommonDef
C      
      IMPLICIT NONE
C	Use Aspen build-in Terminal File Writer Utility to show some message
C	on the control panel
#include "dms_maxwrt.cmn"
C     Include files to pass additional varibles via COMMONs
C     Pass USER_NHSRY
#include "ppexec_user.cmn"
C     Pass arrays containing component data such as molecular weight
#include "dms_plex.cmn"
      real*8 B(1)
      equivalence (B(1), IB(1))
C     PROPERTY PARAMETERS OFFSETS, E.G., FOR MW
#include "dms_ipoff1.cmn"
C     Pass NCOMP_NCC
#include "dms_ncomp.cmn"
C
#include "shs_stwork.cmn"
C
C     Declare arguments
      INTEGER NMATI, NINFI, NMATO, NINFO, NTOT,
     +        NSUBS, NINT,  NPO,   NIWORK, NWORK,
     +        NSIZE
      integer IDSMI(2,NMATI), IDSII(2,NINFI),
     +        IDSMO(2,NMATO), IDSIO(2,NINFO),
     +        IDXSUB(NSUBS), ITYPE(NSUBS), INT(NINT),
     +        IDS(2,3), NBOPST(6,NPO),
     +        IWORK(NIWORK), INTSIZ(NSIZE), NREAL, LD 
      real*8 SIN0(NTOT,NMATI), SINFI(NINFI),
     +       SOUT(NTOT,NMATO), SINFO(NINFO),
     +       WORK(NWORK), SIZE(NSIZE), REAL(NREAL)
C     Declare local variables
      integer OFFSET, IERR, LDATA, KDIAG, IDX(10), NCP, I, J, 
     +        INDEX, LMW, IFAIL
      real*8 SIN(NTOT,NMATI), X(10), Y(10), PHI(10), DPHI(10), FLOW
      integer XMW(NCOMP_NCC), LXMW
      real*8 CMW(NCOMP_NCC)
C     Declare dummies used only for invoking subroutines
      real*8 DUMMY, DUMMY2(2)
C     Declare Aspen I/O functions
      INTEGER USRUTL_GET_INT_PARAM, USRUTL_GET_REAL_PARAM,
     +        USRUTL_SET_INT_PARAM, USRUTL_SET_REAL_PARAM 
C     Side index of hollow fiber tubes, 1 - lumen side; 2 - shell side
      integer ISIDE     
      
C     Initiate membrane module
      call InitMOD(COM_MOD)
C     Initiate streams
      do ISIDE = 1, 2
        call InitStream(COM_SIN(ISIDE))
        call InitStream(COM_SOUT(ISIDE))
      end do
C
C     Get configured integer variables from user2 unit in Aspen
      IFAIL = 0
      INDEX = 0
C      
C     Option for calculating the membrane area
C     OPT1 = 0 - membrane area calculation according to ID1
C            1 -                                        OD1
C            2 -                                        averged diameter
      IERR = USRUTL_GET_INT_PARAM('OPT1', INDEX, COM_OPT(1))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 1'
        IFAIL = 1
      END IF
C     Option for exporting temperature profile
C     OPT2 = 1 - export the temperature profile into the data files
      IERR = USRUTL_GET_INT_PARAM('OPT2', INDEX, COM_OPT(2))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 2'
        IFAIL = 1
      END IF
C     Option for exporting the data files for debugging in vs project
C     OPT3 = 1 - export temporary files
      IERR = USRUTL_GET_INT_PARAM('OPT3', INDEX, COM_OPT(3))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 3'
        IFAIL = 1
      END IF 
C     Option for setting the feeding side in lumen side (OPT4 = 1)
C                                      or in shell side (OPT4 = 2)           
      IERR = USRUTL_GET_INT_PARAM('OPT4', INDEX, COM_OPT(4))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 4'
        IFAIL = 1
      END IF      
      IERR = USRUTL_GET_INT_PARAM('OPT5', INDEX, COM_OPT(5))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 5'
        IFAIL = 1
      END IF      
      IERR = USRUTL_GET_INT_PARAM('OPT6', INDEX, COM_OPT(6))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 6'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_INT_PARAM('OPT7', INDEX, COM_OPT(7))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 7'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_INT_PARAM('OPT8', INDEX, COM_OPT(8))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 8'
        IFAIL = 1
      END IF      
      IERR = USRUTL_GET_INT_PARAM('OPT9', INDEX, COM_OPT(9))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING RUNNING OPTION 9'
        IFAIL = 1
      END IF
C     Get configured real variables from user2 unit in Aspen
      IERR = USRUTL_GET_REAL_PARAM('HFNUMBER', INDEX, COM_MOD%NUM)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING TUBES NUMBER'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('HFLENGTH', INDEX, COM_MOD%LEN)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING TUBE LENGTH'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('HFINNERD', INDEX, COM_MOD%ID1)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING INNER DIAMETER OF TUBES'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('HFOUTERD', INDEX, COM_MOD%OD1)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING OUTER DIAMETER OF TUBES'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('HOUSEID', INDEX, COM_MOD%ID2)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY,*) 'ERROR FETCHING INNER DIAMETER OF SHELL'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('HOUSEOD', INDEX, COM_MOD%OD2)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING OUTER DIAMETER OF SHELL'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('THERMCND', INDEX, 
     +                             COM_MOD%Membrane%KM)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING THERMAL CONDUCTIVITY 
     +                          OF MEMBRANE'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('POROSITY', INDEX, 
     +                             COM_MOD%Membrane%porosity)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING POROSITY'
        IFAIL = 1
      END IF
      IERR = USRUTL_GET_REAL_PARAM('PORESIZE', INDEX, 
     +                             COM_MOD%Membrane%PoreRadius)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR FETCHING PORE SIZE'
        IFAIL = 1
      END IF

C     Get the streams from Aspen Plus
      KDIAG = 4
C     Check whether are there two inlet streams
      IF (NMATI .NE. 2) THEN
        WRITE(USER_NHSTRY, *) 'ERROR: TWO INLET STREAMS NEEDED'
        RETURN
      END IF
C     Check whether are there two outlet streams
      IF (NMATO .NE. 2) THEN
        WRITE(USER_NHSTRY, *) 'ERROR: TWO OUTLET STREAMS NEEDED'
        RETURN
      END IF

C     Set the feeding stream index as declaration of COM_OPT(4)
C     ISIDE = 1 means feeding in lumen side, 
C         and 2 means feeding in shell side
      call SetStreamIndex(COM_OPT(4), SIN0, SIN)

C     Get streams from Aspen Plus
      DO ISIDE = 1, 2
C       Get mass flowrates = (molar flowrate)*(molecular weight), [kg/s]
        COM_SIN(ISIDE)%W = SIN(NCOMP_NCC+1,ISIDE)*SIN(NCOMP_NCC+9,ISIDE)
C       Get temperature, [K]
        COM_SIN(ISIDE)%T = SIN(NCOMP_NCC+2,ISIDE)
C       Get pressure, [Pa]
        COM_SIN(ISIDE)%P = SIN(NCOMP_NCC+3,ISIDE)
C       Get density of liquid mixture
        COM_SIN(ISIDE)%PhysProp%rho = SIN(NCOMP_NCC+8,ISIDE) 
C       Get viscosity of liquid mixture, [Pa-s]
        CALL SHS_CPACK(SIN(1,ISIDE), NCP, IDX, X, FLOW)
        CALL PPMON_VISCL(SIN(NCOMP_NCC+2,ISIDE), SIN(NCOMP_NCC+3,ISIDE),
     +                   X, NCP, IDX, NBOPST, KDIAG, 
     +                   COM_SIN(ISIDE)%PhysProp%mu, IERR)
        IF (IERR .NE. 0) THEN
          WRITE(USER_NHSTRY, *) 'ERROR EVALUATING VISCOSITY'
          IFAIL = 1
        END IF
C       Get thermal conductivity of liquid mixture, [W/m-K]
        CALL PPMON_TCONL(SIN(NCOMP_NCC+2,ISIDE), SIN(NCOMP_NCC+3,ISIDE),
     +                   X, NCP, IDX, NBOPST, KDIAG, 
     +                   COM_SIN(ISIDE)%PhysProp%KM, IERR)
        IF (IERR .NE. 0) THEN
          WRITE(USER_NHSTRY, *) 'ERROR EVALUATING THERMAL CONDUCTIVITY'
          IFAIL = 1
        END IF
C       Get specific heat of liquid mixture, [J/kg-K]
C       Enthalpy monitor is called with KH=2 to compute the specific heat
C       Refer to "Aspen Properties: toolkit manual" P57
C       The cp output is in [J/kmol-K],
C       ref to "Aspen Properties: toolkit manual" P24
        CALL PPMON_ENTHL(SIN(NCOMP_NCC+2,ISIDE), SIN(NCOMP_NCC+3,ISIDE),
     +                   X, NCP, IDX, NBOPST, KDIAG, 0, 2, DUMMY,
     +                   COM_SIN(ISIDE)%PhysProp%cp, IERR)
        IF (IERR .NE. 0) THEN
          WRITE(USER_NHSTRY, *) 'ERROR EVALUATING SPECIFIC HEAT'
          IFAIL = 1
        END IF
C       Get the average molecular weight
!        COM_SIN(ISIDE)%PhysProp%AMW = PPUTL_AVEMW(NCP, IDX, X)
        LXMW = IPOFF1_IPOFF1(306)
        do i = 1, NCOMP_NCC
          XMW(i) = LXMW+i
          CMW(i) = B(XMW(i))
        end do
        COM_SIN(ISIDE)%PhysProp%AMW = 
     +                  AvgMolWeight(NCP, CMW(IDX), X(1:NCP))
!        WRITE(MAXWRT_MAXBUF, "(2E12.4)") (CMW(IDX(i)), X(i), i=1, NCP)
C       Convert the unit from [J/kmol-K] to [J/kg-K]
        COM_SIN(ISIDE)%PhysProp%cp = COM_SIN(ISIDE)%PhysProp%cp/
     +                               COM_SIN(ISIDE)%PhysProp%AMW
      END DO

C     Run the simulation
      call CalcSOUT(NTOT, NMATI, SIN)

C     Fill user2 unit with resulted parameters
C     Set integer variables
      IERR = USRUTL_SET_INT_PARAM('OPT10', INDEX, COM_OPT(10))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 10'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_INT_PARAM('OPT11', INDEX, COM_OPT(11))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 11'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_INT_PARAM('OPT12', INDEX, COM_OPT(12))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 12'
        IFAIL = 1
      END IF 
      IERR = USRUTL_SET_INT_PARAM('OPT13', INDEX, COM_OPT(13))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 13'
        IFAIL = 1
      END IF                 
      IERR = USRUTL_SET_INT_PARAM('OPT14', INDEX, COM_OPT(14))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 14'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_INT_PARAM('OPT15', INDEX, COM_OPT(15))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 15'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_INT_PARAM('OPT16', INDEX, COM_OPT(16))
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING OPTION 10'
        IFAIL = 1
      END IF                  
C     Set real variables
      IERR = USRUTL_SET_REAL_PARAM('RE1', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING FEEDING SIDE REYNOLDS
     +                         NUMBER'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_REAL_PARAM('NU1', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING FEEDING SIDE NUSSELT
     +                         NUMBER'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_REAL_PARAM('RE2', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING PERMEATING SIDE REYNOLDS
     +                          NUMBER'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_REAL_PARAM('NU2', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING PERMEATING SIDE NUSSELT
     +                          NUMBER'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_REAL_PARAM('JM', INDEX, COM_MOD%Performance%JM)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING AVG. PERMEATION FLUX'
        IFAIL = 1
      END IF      
      IERR = USRUTL_SET_REAL_PARAM('VAR1', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING MONITOR VARIABLE 1'
        IFAIL = 1
      END IF
      IERR = USRUTL_SET_REAL_PARAM('VAR2', INDEX, 0.0)
      IF (IERR .NE. 0) THEN
        WRITE(USER_NHSTRY, *) 'ERROR STORING MONITOR VARIABLE 2'
        IFAIL = 1
      END IF

C     Set outlet streams
!	SOUT = SIN
      DO ISIDE = 1, 2
C       Water molar flow, [kmol/s]      
        SOUT(1,ISIDE) = COM_SOUT(ISIDE)%MolarFlow%H2O
C       Other components' molar flow, [kmol/s]
        DO I = 2, NCOMP_NCC
          SOUT(I,ISIDE) = SIN(I,ISIDE)
        END DO
C       Total molar flow, [kmol/s]
        SOUT(NCOMP_NCC+1,ISIDE) = SUM(SOUT(1:NCOMP_NCC,ISIDE))
C       Temperature, [K]
        SOUT(NCOMP_NCC+2,ISIDE) = COM_SOUT(ISIDE)%T
C       Pressure, [Pa]
        SOUT(NCOMP_NCC+3,ISIDE) = COM_SOUT(ISIDE)%P
      END DO

!      WRITE(MAXWRT_MAXBUF, "(2E12.4)") (COM_SOUT(ISIDE)%MolarFlow%H2O, 
!     +                                   ISIDE = 1, 2)
      CALL DMS_WRTTRM(11)

      RETURN
      END


