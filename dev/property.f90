    module VaporConc
    ! Calculate the variation of vapor molar concentration with respect to membrane temperature, d(nw)/dT
    ! Ideal gas is applied to model the P-V-T of vapor, i.e., nw = p/R/T
    ! Thus, d(nw)/dT = (1/T dp/dT-p/T^2)/R
    ! Clausius-Clapeyron equation is used to correlate the derivation with respect to temperature
    ! dp/dT = p\frac{l(T)}{RT^2}, where molar latent heat of vapor, l(T), is obtained from data fitting
    ! Calculate the molar latent heat of saturated water based on data fitting 
    ! at temperature range of 0-100 C
    ! Kim A.S., Journal of membrane science, 428(2013), 410-424
    ! [2016/01/18] The saturation vapor pressure (SVP) of water, correlated by Kim, produces poor outputs, which are far away from the experimental data. The 5-parameter model was finally used to obtain the SVP and deriviation of SVP.
    
        use CommonDef
        use ToolKits

        real*8, parameter :: L0 = 57.075, L1 = 4.3856d-2, P0 = 1.0684d4, R = 8.314

        contains
!       
!        real(8) function dnw(T) ! derivation of molar concentration
!            real(8), intent(in) :: T ! Temperature [K] 
!            real(8) :: p ! vapor pressure corresponding to T [Pa]
!            ! For ideal gas
!            call SatVapPres(2, T, p)
!            rho = SteamDensity(T) 
!            dnw = 1/(R)*(dP(T)/T-p/(T*T))
!        end function

        real(8) function diffnw(temp) ! derivation of molar concentration
          real(8), intent(in) :: temp ! Temperature [K] 
          real(8) :: p ! vapor pressure corresponding to T [Pa]
          real(8) :: diffP ! deriviation of saturated vapor pressure [Pa]
          ! For ideal gas
!         Define the arguments required for subroutine DIFF()
          real :: T, Tmin, Tmax, eps, acc, deriv, err
          integer :: IORD, IFAIL
          Tmin = 0.
          Tmax = 100.
          eps = 1.e-5
          acc = 0.
          IORD = 1
!         Check the input temperature's unit
          if (temp .LE. 273.15) then
            T = temp
          else
            T = temp-273.15
          end if
!         Invoke the subroutine diff() to calculate the deriviation of molar concentration with respect to temperature
!         Due to the default REAL used in the subroutine diff(), a interface function r4MolConc() is required.
          call diff(IORD, T, Tmin, Tmax, r4MolConc, eps, acc, deriv, err, IFAIL)
          diffnw = deriv*2d1
        end function
        
        real(8) function MolarConc(T)
          real(8), intent(in) :: T
          real(8) :: SVV
          integer :: opt = 2
          call SpecVolV(opt, T, SVV)
          MolarConc = one/SVV/18.*1d3 ! molar concentration of steam [mol/m3]
        end function

        real function r4MolConc(r4T)
          real, intent(in) :: r4T
          real(8) :: SVV, T
          integer :: opt = 2
          T = r4T
          call SpecVolV(opt, T, SVV)
          r4MolConc = real(one/SVV/18.*1d3)
        end function

!       Specific volume of saturated vapor
        subroutine SpecVolV(opt, temp, SVV, IFAIL)
          integer, intent(in) :: opt ! options for using specific interpolation method
          real*8, intent(in) :: temp ! input temperature [K]
          real*8, intent(out) :: SVV ! output specific volume of vapor [m3/kg]
          integer, intent(out), optional :: IFAIL ! options for running
          real*8 :: T
          real*8, pointer :: xd(:), yd(:), xi(:), yi(:), cd(:)
          integer :: nd, ni
          ! Initiation
          if (present(IFAIL)) IFAIL = 0
          SVV = zero
          call def_PropWaterSteam()
          xd => COM_PWS%t
          yd => COM_PWS%vv
          nd = size(COM_PWS)
          ni = 1
          allocate(xi(ni), yi(ni), cd(nd))
          ! Check the input temperature's unit [C]
          call CheckTempUnit('C', temp, T)
          ! Check the input temperature's range
          if (T .GE. 373.15) then
            if (present(IFAIL)) IFAIL = 2
            return
          end if
          xi(ni) = T

          select case(opt)
            case(1)
              ! Use 1d piecewise linear interpolation              
              call pwl_value_1d(nd, xd, yd, ni, xi, yi)
            case(2)
              ! Use 1d Newton interpolation
              call newton_coef_1d(nd, xd, yd, cd)
              call newton_value_1d(nd, xd, cd, ni, xi, yi)
          end select
          SVV = yi(ni)
          deallocate(xi, yi, cd)
        end subroutine

        subroutine def_PropWaterSteam()
          type(Steam) :: PWS(10)
          data PWS( 1)%t, PWS( 1)%p, PWS( 1)%vv / 10,   1228, 106.31/
          data PWS( 2)%t, PWS( 2)%p, PWS( 2)%vv / 20,   2339, 57.761/
          data PWS( 3)%t, PWS( 3)%p, PWS( 3)%vv / 30,   4247, 32.882/
          data PWS( 4)%t, PWS( 4)%p, PWS( 4)%vv / 40,   7384, 18.517/
          data PWS( 5)%t, PWS( 5)%p, PWS( 5)%vv / 50,  12351, 12.028/
          data PWS( 6)%t, PWS( 6)%p, PWS( 6)%vv / 60,  19946, 7.6677/
          data PWS( 7)%t, PWS( 7)%p, PWS( 7)%vv / 70,  31201, 5.0397/
          data PWS( 8)%t, PWS( 8)%p, PWS( 8)%vv / 80,  47415, 3.4053/
          data PWS( 9)%t, PWS( 9)%p, PWS( 9)%vv / 90,  70182, 2.3591/
          data PWS(10)%t, PWS(10)%p, PWS(10)%vv /100, 101420, 1.6719/
          COM_PWS = PWS
        end subroutine

!       Vapor pressure of saturated steam
        subroutine SatVapPres(opt, temp, SVP, IFAIL)
          integer, intent(in) :: opt ! options for using specific correlation
          real*8, intent(in) :: temp ! input temperature [K]
          real*8, intent(out) :: SVP ! output saturation vapor pressure [Pa]
          integer, intent(out), optional :: IFAIL ! options for running result, if no error, IFAIL = 0
          real*8 :: T
          real*8 :: A(3)
          real*8 :: C(5)
          data A /23.238, 3841.0, -45.0/
          data C /73.649, -7258.2, -7.3037, 4.1653E-6, 2.0/
          ! Initiation
          SVP = zero
          if (present(IFAIL)) IFAIL = 0
          ! Check the input temperature's unit [K]
          IF (temp .LE. 273.15) THEN
            T = 273.15+temp
          END IF
          ! Check the input temperature's range
          if (T .GE. 373.15) then
            if (present(IFAIL)) IFAIL = 2
            return
          end if
          select case(opt)
            case(1)
              ! Antonie correlation
              SVP = DEXP(A(1)-A(2)/(A(3)+T))
            case(2)
              SVP = DEXP(C(1)+C(2)/T+C(3)*DLOG(T)+C(4)*T**C(5))
            case(3)
              ! Kim correlation
              SVP = T**(-L1/R)*dexp(dlog(p0)-L0/R/T)
          end select
        end subroutine

        real function SVP1(temp)
          real, intent(in) :: temp
          real*8 :: T, Psat
          T = Temp
          call SatVapPres(1, T, Psat)
          SVP1 = Psat
        end function

        real function SVP2(temp)
          real, intent(in) :: temp
          real*8 :: T, Psat
          T = Temp
          call SatVapPres(2, T, Psat)
          SVP2 = Psat
        end function

        real function SVP3(temp)
          real, intent(in) :: temp
          real*8 :: T, Psat
          T = Temp
          call SatVapPres(3, T, Psat)
          SVP3 = Psat
        end function

        real(8) function tortuosity(porosity)
        ! Calculate the tortuosity for highly interconnected pores
        ! The analytical approach is developed in Beekman's work
        ! Cited in A.S. Kim, Journal of membrane science, 455(2014) 168-186
            real*8, intent(in) :: porosity
            real*8 :: epsilon, tau
            epsilon = porosity
            if ((porosity .ge. zero).and.(porosity .le. one))  then
            tau = epsilon/(one-(one-epsilon)**(one/three))
            tortuosity = tau
            else
            tortuosity = zero
            end if
        end function
        
        real(8) function EffThermalCond(T)
        !Calculate the conductivity of the membrane
        !Since the conductivity is related with the solid conductivity and gas conductivity ,the gas conductivity can be calculated with equations
        !solid conductivity is assumed approximately constant over the temperature range of DCMD operation
        !porostiy is determined by the membrane material
            real(8), intent(in) :: T
            real(8) :: porosity, solidConductivity
            real(8) :: beta, gasConductivity, adjustT
            !K0,K1,K2,K3,K4 are constant suggested by Kim
            real(8) :: K0 = 2.40073953d-2, K1 = 7.278410162d-5, K2 = -1.788037411d-20, K3 = -1.351703529d-9, K4 = -3.322412767d-11 
            
            if (CheckRange(T, 2.7315d2, 3.7315d2)) then ! the input temperature is absolute temperature
            !   Check the input argument in the specific range
                if (CheckRange(UnitConvert(T, "K", "C"), 0.d0, 1.d2)) then
                ! Convert temperature unit to Celcius degree
                     adjustT = UnitConvert(T, "K", "C")
                else
                     adjustT = zero
                end if                
            else ! the input temperature is in Celcius degree
            !   Check the input argument in the specific range after unit conversion
                if (CheckRange(T, 0.d0, 1.d2)) then
                    adjustT = T
                else
                    adjustT = zero
                end if                  
            end if
            ! Get the solid thermal conductivity of membrane and porosity
            porosity = COM_MOD%Membrane%Porosity
            solidConductivity = COM_MOD%Membrane%KM
            ! Correlate the vapor thermal conductivity
            gasConductivity = K0 + K1*adjustT + K2*adjustT**TWO + K3*adjustT**THREE + K4*adjustT**FOUR
            ! Calculate the effective conductivity for porosity higher than 0.6
            beta = (solidConductivity - gasConductivity)/(solidConductivity + TWO*gasConductivity)
            EffThermalCond = gasConductivity*(ONE+TWO*beta*(ONE-porosity))/(ONE-beta*(ONE-porosity))
        end function
        
        real(8) function EffMolEnthalpy(T)
        !Calculate the effective molar enthalpy    
            real*8,intent(in) :: T
            real*8 :: H0 = 2.4527d1 , Cp = 5.2434d-2
            ! temperature can be in absolute degree or Celcius degree
            !   Check the input temperture to be absolute temperature
            if (CheckRange(T, 2.7315d2, 3.7315d2)) then ! the input temperature is absolute temperature
            !   Check the input argument in the specific range
                if (CheckRange(UnitConvert(T, "K", "C"), 0.d0, 1.d2)) then
                    EffMolEnthalpy = H0 - Cp*T
                else
                    EffMolEnthalpy = zero
                end if                
            else ! the input temperature is in Celcius degree
            !   Check the input argument in the specific range after unit conversion
                if (CheckRange(T, 0.d0, 1.d2)) then
                    EffMolEnthalpy = H0 - Cp*UnitConvert(T, "C", "K")
                else
                    EffMolEnthalpy = zero
                end if                  
            end if
        end function
        
        real(8) function EffDiffusivity(T)
        !Calculate the effective diffusion coefficient
            real*8,intent(in) :: T
            real*8 :: PoreRadius
            real*8 :: a = 1.072, PT = 1.01325d5,R = 8.314, Mw = 1.801524d1, B, Db , Dk, K
            PoreRadius = COM_MOD%Membrane%PoreRadius
                ! temperature can be in absolute degree or Celcius degree
        !   Check the input temperture to be absolute temperature
            if (CheckRange(T, 2.7315d2, 3.7315d2)) then ! the input temperature is absolute temperature
            !   Check the input argument in the specific range
                if (CheckRange(UnitConvert(T, "K", "C"), 0.d0, 1.d2)) then
                    B = 1.895d-5/PT
                    Db = B * T**(1+a)
                    K = ONE/THREE*poreRadius*DSQRT(8.0*R/(pi*Mw))
                    Dk = K * DSQRT(T)
                    EffDiffusivity = ONE/(ONE/Db+ONE/Dk)
                else
                    EffDiffusivity = zero
                end if                
            else ! the input temperature is in Celcius degree
            !   Check the input argument in the specific range after unit conversion
                if (CheckRange(T, 0.d0, 1.d2)) then
                    B = 1.895d-5/PT
                    Db = B * UnitConvert(T, "C", "K")**(1+a)
                    K = ONE/THREE*poreRadius*DSQRT(8.0*R/(3.1415926*Mw))
                    Dk = K * DSQRT(UnitConvert(T, "C", "K"))
                    EffDiffusivity = ONE/(ONE/Db+ONE/Dk)
                else
                    EffDiffusivity = zero
                end if                  
            end if
        end function

        ! The thermal conducticity of water was developed by Ramires et al. [Kim 2014 JMS]
        real(8) function WaterConductivity(T) ! [W/m-K]
            real*8,intent(in) :: T ! [K]
            WaterConductivity = 6.065d-1 * (-1.48445 + 4.12292*T/298.15 -1.63866*(T/298.15)**TWO)
        end function
                
    end module