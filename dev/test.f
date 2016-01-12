      subroutine test_CalcProfile
        use CommonDef
        use VaporConc
        type(StreamProp) :: S1, S2, S3, S4
!        ! Construct the membrane module
!        COM_MOD%NUM = 8.4d2
!        COM_MOD%LEN = 0.3867
!        COM_MOD%ID1 = 0.98d-3
!        COM_MOD%OD1 = 1.46d-3
!        COM_MOD%PHI = 0.5
!        COM_MOD%Membrane%porosity = 0.7
!        COM_MOD%Membrane%PoreRadius = 0.15d-6
!        call CalcModule(COM_MOD)
!        ! Construct the stream
!        S1%W = 1.0d-1 ! lumen-side mass flowrate is 0.1 kg/s
!        S1%T = UnitConvert(8.0d1, "C", "K")
!        S1%P = 2.0d+5 ! lumen-side inlet pressure is 200 kPa
!        S1%PhysProp%rho = 1.d3
!        S1%PhysProp%mu = 1.d-3
!        S1%PhysProp%cp = 4.12d3 ! specific heat is 4.12e3 J/kg-K
!        S1%PhysProp%KM = WaterConductivity(S1%T)
!        call CalcStream(S1, COM_MOD%CSA1)
!        S2 = S1
!        S2%T = UnitConvert(4.0d1, "C", "K")
!        call CalcStream(S2, COM_MOD%CSA2)

        ! get MD module parameters
        call ReadMOD(mod_def_file)
        ! get feeding streams
        call ReadStream(stm_def_file)
       
        S1 = COM_SIN(1)
        S2 = COM_SIN(2)
        call CalcProfile(S1, S2, S3, S4)
        COM_SOUT(1) = S3
        COM_SOUT(2) = S4
        
      end subroutine
      
      subroutine test_DAVG
        use CommonDef
        real*8, dimension(5) :: y
        y = (/2.1, 4.2, 9.5, 6.3, 8.7/)
        write(*, *) DAVG(y)
      end subroutine
      
      subroutine test_MolarLatentHeat
      use CommonDef
      use VaporConc
      real*8 :: T, MLH
!     Initiate the temperature
!      T = UnitConvert(6.d1, "C", "K") ! absolute temperature
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      MLH = MolarLatentHeat(T)
      write (*, *) "Known molar latent heat of water in 60 C"
      write (*, *) "42.624 kJ/mol"
      write (*, *) "The calculated value is"
      write (*, *) MLH
      pause
      end subroutine
      
      subroutine test_PressureDerivation_dP(T)
      use CommonDef
      use vaporConc
      real*8 :: T,PreDe  
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      PreDe  = dP(T)
      write (*, *) "The calculated value is"
      write (*, *) PreDe 
      pause
      end subroutine
      
      subroutine test_powerfunction_L(T)
      use CommonDef
      use vaporConc
      real*8 :: T,POWERFUNC  
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      POWERFUNC = L(T)
      write (*, *) "The calculated value is"
      write (*, *) POWERFUNC
      pause
      end subroutine
      
      subroutine test_VaporPressure
      use CommonDef
      use vaporConc
      real*8 :: T, VP
!     Initiate the temperature
!      T = UnitConvert(6.d1, "C", "K") ! absolute temperature
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      VP = VaporPressure(T)
      write (*, *) "Known Vapor Pressure of saturated water in 50 C"
      write (*, *) "12334 Pa"
      write (*, *) "The calculated value is"
      write (*, *) VP
      pause
      end subroutine
      
      subroutine test_tortuosity(porosity)
      use CommonDef
      use vaporConc
      real*8 :: porosity,tortuo
      write(*, *) "Input the porosity for evaluation"
      read(*, *) porosity
      tortuo = tortuosity(porosity)
      write (*, *) "The calculated value is"
      write (*, *) tortuo
      pause
      end subroutine
            

      
      subroutine test_EffThermalConductivity
      use CommonDef
      use vaporConc
      real*8 :: T,porosity, solidConductivity,ETherCon
      write(*, *) "Input the posrosity for evaluation"
      read(*, *) porosity
      write(*, *) "Input the SolidConductivity for evaluation"
      read(*, *)  solidConductivity
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      ETherCon = EffThermalCond(T)
      write (*, *) "The calculated value is"
      write (*, *) ETherCon
      pause
      end subroutine
      
      subroutine test_EffectiveMolarEnthalpy
      use CommonDef
      use vaporConc
      real*8 :: T,EMolarEn
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      EMolarEn = EffMolEnthalpy(T)
      write (*, *) "The calculated value is"
      write (*, *) EMolarEn
      pause
      end subroutine      
      
      subroutine test_ EffectiveDiffusivity
      use CommonDef
      use vaporConc
      real*8 :: T,EfDiffu
      write(*, *) "Input the PoreRadius for evaluation"
      read(*, *) PoreRadius
      write(*, *) "Input the temperture for evaluation"
      read(*, *) T
      EfDiffu = EffDiffusivity(T)
      write (*, *) "The calculated value is"
      write (*, *) EfDiffu
      pause
      end subroutine         
      
!      subroutine test_qag
!      use CommonDef
!!     define the arguments for the solver
!      real, external :: f, g1
!      real :: a, b, epsabs, epsrel
!      real :: result, abserr
!      integer :: key
!      integer :: neval, ier
!!     set the arguments
!      epsabs = 1.d-6
!      epsrel = 5.d-4
!      key = 1
!      a = 40
!      b = 80
!!     invoke the solver      
!      call qag  (f, a, b, epsabs, epsrel, key, 
!     &           result, abserr, neval, ier )
!!     output the results
!      write(*, *) "Integral of f(x)=SIN(x) for x=0..pi/2"
!      write(*, '(A40, F6.4)') "Known result of above integral: ", 1.0
!      write(*, '(A40, F6.4)') "Calculated result of above integral: ", 
!     &                          result
!      write(*, '(A40, F6.4)') "Relative error: ", (1.0-result)/1.0
!
!
!      pause
!                 
!      end subroutine test_qag

      subroutine test_dqk15
      use CommonDef
!     define the arguments for the solver
      real*8, external :: g
      real*8 :: a, b
      real*8 :: result, abserr, resabs, resasc
      a = zero
      b = pi/two
!     invoke the solver      
      call dqk15(g, a, b, result, abserr, resabs, resasc)
!     output the results
      write(*, *) "Integral of f(x)=SIN(x) for x=0..pi/2"
      write(*, '(A40, F7.4)') "Known result of above integral: ", 1.0
      write(*, '(A40, F7.4)') "Calculated result of above integral: ", 
     &                          result
      write(*, '(A40, F7.4)') "Relative error: ", (1.0-result)/1.0


      pause
                 
      end subroutine test_dqk15
      
      function f(x)
      real x
      f = sin(x)
      end function
      
      real*8 function g(x)
      real*8, intent(in) :: x
      g = DCOS(x)
      end function

      
      subroutine test_moduleVaporConc
      use VaporConc
      real*8 estT, estP
      estT = 6.0d1
      estP = 7.6d2
      write(*, *) "Check the result of vapor pressure:"
      write(*, *) VaporPressure(estT)
      pause
      end subroutine test_moduleVaporConc     
      
!      subroutine test_conductivity
!      use VaporConc
!      real*8 T,porosity,solidConductivity
!      T = 2.96d2
!      porosity = 0.8
!      solidConductivity = 0.2 
!      write(*,*) "check the result of conductivity: "
!      write(*,*) ThermalConductivity(T,porosity,solidConductivity)
!      pause
!      end subroutine test_conductivity
      
!      subroutine test_heatTransferRate
!      use HeatTransferRate
!      real*8 :: a,b,porosity,solidConductivity,poreRadius
!      real :: Ta,Tb
!      a=0.5
!      b=1.0
!      Ta=2.96d2
!      Tb=3.30d2
!      porosity = 0.8
!      solidConductivity = 0.2
!      poreRadius = 0.001
!      write(*,*) "check the result of heat transfer rate: "
!      write(*,*) Sq(a,b,Ta,Tb,porosity,solidConductivity,poreRadius)
!      pause
!      end subroutine test_heatTransferRate

    
      
      subroutine ReadMOD(DataFileName)
        use CommonDef
        use toolkits
        
        character*(*), intent(in) :: DataFileName
        integer :: i
        real, dimension(21) :: mod_data
        character*31 :: dummy_str
        
        open(11, file = DataFileName, action = 'read')
        do i = 1, 21
          read(11, mod_def_file_fmt) dummy_str, mod_data(i)
        end do
        close(11)
        
        COM_MOD%NUM = mod_data(1) 
        COM_MOD%LEN = mod_data(2)
        COM_MOD%ID1 = mod_data(3)
        COM_MOD%OD1 = mod_data(4)
        COM_MOD%ID2 = mod_data(5)
        COM_MOD%OD2 = mod_data(6)
        COM_MOD%THK = mod_data(7)
        COM_MOD%PHI = mod_data(8)
        COM_MOD%AM = mod_data(9)
        COM_MOD%CSA1 = mod_data(10)
        COM_MOD%CSA2 = mod_data(11)
        COM_MOD%Membrane%rho = mod_data(12)
        COM_MOD%Membrane%mu = mod_data(13)
        COM_MOD%Membrane%cp = mod_data(14)
        COM_MOD%Membrane%KM = mod_data(15)
        COM_MOD%Membrane%CM = mod_data(16)
        COM_MOD%Membrane%porosity = mod_data(17)
        COM_MOD%Membrane%tortuosity = mod_data(18)
        COM_MOD%Membrane%PoreRadius = mod_data(19)
        COM_MOD%Performance%JM = mod_data(20)
        COM_MOD%Performance%eta = mod_data(21)
        
      end subroutine
      
      subroutine ReadStream(DataFileName)
        use CommonDef
        use toolkits
        
        character*(*), intent(in) :: DataFileName
        real, dimension(15,2) :: stream_data
        integer :: ISIDE
        character*31 :: dummy_str
        
        open(11, file = DataFileName, action = 'read')
        do i = 1, 15
          read(11, stm_def_file_fmt) dummy_str, stream_data(i,1), 
     +                                    stream_data(i,2)
        end do
        close(11)
        
        do ISIDE = 1, 2
          COM_SIN(ISIDE)%W = stream_data(1,ISIDE)
          COM_SIN(ISIDE)%T = stream_data(2,ISIDE)
          COM_SIN(ISIDE)%P = stream_data(3,ISIDE)
          COM_SIN(ISIDE)%G = stream_data(4,ISIDE)
          COM_SIN(ISIDE)%u = stream_data(5,ISIDE)
          COM_SIN(ISIDE)%Re = stream_data(6,ISIDE)
          COM_SIN(ISIDE)%Nu = stream_data(7,ISIDE)
          COM_SIN(ISIDE)%PhysProp%rho = stream_data(8,ISIDE)
          COM_SIN(ISIDE)%PhysProp%mu = stream_data(9,ISIDE)
          COM_SIN(ISIDE)%PhysProp%cp = stream_data(10,ISIDE)
          COM_SIN(ISIDE)%PhysProp%KM = stream_data(11,ISIDE)
          COM_SIN(ISIDE)%MassFrac%H2O = stream_data(12,ISIDE)
          COM_SIN(ISIDE)%MassFrac%NaCl = stream_data(13,ISIDE)
          COM_SIN(ISIDE)%MolarFlow%H2O = stream_data(14,ISIDE)
          COM_SIN(ISIDE)%MolarFlow%NaCl = stream_data(15,ISIDE)  
        end do
               
      end subroutine