subroutine CalcProfile(LumenIN, ShellIN, LumenOUT, ShellOUT)
  use CommonDef
  use toolkits

  type(StreamProp), intent(in) :: LumenIN ! lumen-side inlet stream
  type(StreamProp), intent(in) :: ShellIN ! shell-side inlet stream
  type(StreamProp), intent(out) :: LumenOUT ! lumen-side outlet stream
  type(StreamProp), intent(out) :: ShellOUT ! shell-side outlet stream
  
  type(StreamProp) :: LumenStream(COM_RadialGridNum,COM_AxialGridNum)
  type(StreamProp) :: ShellStream(COM_RadialGridNum,COM_AxialGridNum)

  real*8 :: LumenRadialLoc(COM_RadialGridNum), LumenAxialLoc(COM_AxialGridNum)
  real*8 :: ShellRadialLoc(COM_RadialGridNum), ShellAxialLoc(COM_AxialGridNum)
  integer :: opt
  integer :: i, j
  real*8 :: Ta, Tb, Sq(COM_AxialGridNum), Fw(COM_AxialGridNum), AvgSq, PrevIterAvgSq
  real*8 :: dT_lmn, dT_shl
  real*8 :: a, b, c, L
  real*8 :: sigma, eta, zeta, xi
  real*8 :: F_eta, H_xi
  real*8 :: kappa_lmn, rho_lmn, cp_lmn, v_lmn, alpha_lmn, beta_lmn, Pe_lmn
  real*8 :: kappa_shl, rho_shl, cp_shl, v_shl, alpha_shl, beta_shl, Pe_shl
  real*8 :: JM

!  open(12, file = "tmp_CalcProf.log")
!  write(12, *) "Invoking CalcProfile()" 
  
  ! Retrieve the geometrics of hollow fibers
  a = COM_MOD%ID1/two
  b = COM_MOD%OD1/two
  L = COM_MOD%LEN
  sigma = one/dsqrt(COM_MOD%PHI)
      
  ! Grid the location arrays
!  call Grid1D(zero, one, LumenRadialLoc, opt)
  LumenRadialLoc(1) = zero
  do i=2, COM_RadialGridNum
    LumenRadialLoc(i) = LumenRadialLoc(i-1)+(one-zero)/(COM_RadialGridNum-1)
  end do
!  call Grid1D(zero, one, LumenAxialLoc, opt)
  LumenAxialLoc(1) = zero
  do i=2, COM_AxialGridNum
    LumenAxialLoc(i) = LumenAxialLoc(i-1)+(one-zero)/(COM_AxialGridNum-1)
  end do
!  call Grid1D(one, sigma, ShellRadialLoc, opt)
  ShellRadialLoc(1) = one
  do i=2, COM_RadialGridNum
    ShellRadialLoc(i) = ShellRadialLoc(i-1)+(sigma-one)/(COM_RadialGridNum-1)
  end do
!  call Grid1D(zero, one, ShellAxialLoc, opt)
  ShellAxialLoc(1) = zero
  do i=2, COM_AxialGridNum
    ShellAxialLoc(i) = ShellAxialLoc(i-1)+(one-zero)/(COM_AxialGridNum-1)
  end do
  
  ! Initiate the flow field
  do j = 1, COM_AxialGridNum
    do i = 1, COM_RadialGridNum
      LumenStream(i,j) = LumenIn
      ShellStream(i,j) = ShellIn
    end do
  end do

  PrevIterAvgSq = zero
  do IterNum = 1, MaxIterNum
    ! Calculate the local heat and mass transfer rate per unit length
    do j = 1, COM_AxialGridNum
      Ta = LumenStream(1,j)%T
      Tb = ShellStream(1,j)%T
      call LocalRate(Ta, Tb, Fw(j), Sq(j))
    end do
    
    ! Calculate the average local heat transfer rate per unit length
    AvgSq = DAVG(Sq)
    
    ! Check for converge
    if (IsRelDiff(AvgSq, PrevIterAvgSq, 1.d-6)) then
!      write(12, *) "Averaged heat flux is convergent."
!      write(12, "(e12.4)") AvgSq
      exit
    else
!      write(12, "(A, I5, A, E12.4)" ) "Iteration # ", IterNum, "  Averaged heat flux: ", AvgSq
      PrevIterAvgSq = AvgSq
    end if

    ! Calculate the temperature changes in both lumen and shell sides
    ! Retrive the averaged thermal conductivity, density, specific heat and velocity of the lumen-side stream along the central line
    kappa_lmn = DAVG(LumenStream(1,:)%PhysProp%KM)
    rho_lmn = DAVG(LumenStream(1,:)%PhysProp%rho)
    cp_lmn = DAVG(LumenStream(1,:)%PhysProp%cp)
    v_lmn = DAVG(LumenStream(1,:)%u)
    ! Retrive the averaged thermal conductivity, density, specific heat and velocity of the shell-side stream along the boundary line
    kappa_shl = DAVG(ShellStream(COM_RadialGridNum,:)%PhysProp%KM)
    rho_shl = DAVG(ShellStream(COM_RadialGridNum,:)%PhysProp%rho)
    cp_shl = DAVG(ShellStream(COM_RadialGridNum,:)%PhysProp%cp)
    v_shl = DAVG(ShellStream(COM_RadialGridNum,:)%u) ! the average bulk speed
    v_shl = v_shl/(sigma**two-one) ! the average outer speed Eq.(21)
!    write(12, "(A, 2E12.4)") "Velocity ", v_lmn, v_shl
    ! Calculate the thermal diffusivity in both lumen and shell sides
    alpha_lmn = alpha(kappa_lmn, rho_lmn, cp_lmn)
    alpha_shl = alpha(kappa_shl, rho_shl, cp_shl)
!    write(12, "(A, 2E12.4)") "Thermal diffusivities: ", alpha_lmn, alpha_shl
    ! Calculate the Peclet numbers in both lumen and shell sides
    Pe_lmn = Pe(v_lmn, L, alpha_lmn)
    Pe_shl = Pe(v_shl, L, alpha_shl)
!    write(12, "(A, 2E12.4)") "Pelect numbers: ", Pe_lmn, Pe_shl
    ! Calculate the beta values in both lumen and shell sides
    beta_lmn = Pe_lmn*(a/L)**two ! Eq.(12)
    PackEffect = one/(sigma**four*dlog(sigma**four)+four*sigma**two-three*sigma**four-one) ! Eq.(20)
    beta_shl = Pe_shl*(b/L)**two*PackEffect*(sigma**two-one) ! Eq.(24)
!    write(12, "(A, 2E12.4)") "Eeta values: ", beta_lmn, beta_shl
    ! Calculate the temperature differences in both lumen and shell sides
    dT_lmn = four*AvgSq/kappa_lmn/beta_lmn 
    dT_shl = four*AvgSq/kappa_shl/beta_shl
!    write(12, "(A, 2E12.4)") "Temperature differences: ", dT_lmn, dT_shl 
    
    ! Calculate the temperature profile
    do j = 1, COM_AxialGridNum
      do i = 1, COM_RadialGridNum
        eta = LumenRadialLoc(i)
        zeta = LumenAxialLoc(j)
        F_eta = eta**two-one/four*eta**four ! Eq.(A.10)
        LumenStream(i,j)%T = LumenIn%T-dT_lmn*(LumenAxialLoc(j)+one/four*beta_lmn*F_eta) ! Eq.(61)
        xi = ShellRadialLoc(i)
        zeta = ShellAxialLoc(j)
        G_xi = (sigma**two*(two*sigma**two*dlog(sigma**two)-three*sigma**two+four)-one)
        H_xi = (G(sigma)-G(xi))/(sigma**two*(two*sigma**two*dlog(sigma**two)-three*sigma**two+four)-one) ! Eq.(65) and (A.27)
        ShellStream(i,j)%T = ShellIn%T+dT_shl*(1-zeta+one/four*beta_shl*H_xi) ! Eq.(62)
      end do
    end do
  end do
    
! call Export2DArray(LumenStream%T, 'Tlmn.txt')
! call Export2DArray(ShellStream%T, 'TShl.txt')
! Export data file for COM_OPT(2) = 1
  if (COM_OPT(2) .eq. 1) then
    call ProfileExport(LumenRadialLoc, LumenAxialLoc, LumenStream%T, 'T_inr.dat')
    call ProfileExport(ShellRadialLoc, ShellAxialLoc, ShellStream%T, 'T_otr.dat')
  end if
    
! Calculate outlet averaged temperature in lumen side
  LumenOUT%T = DAVG(LumenStream(:,1)%T)
! Calculate outlet averaged temperature in shell side
  ShellOUT%T = DAVG(ShellStream(:,1)%T)
  
! Calculate permeation flux [kg/m2-s]
  COM_MOD%Performance%JM = DAVG(Fw)/b*18.0
  
! Calculate the outlet streams
  LumenOUT%MolarFlow%H2O = LumenIN%MolarFlow%H2O+COM_MOD%AM*DAVG(Fw)/b
  ShellOUT%MolarFlow%H2O = ShellIN%MolarFlow%H2O-COM_MOD%AM*DAVG(Fw)/b
  LumenOUT%MolarFlow%NaCl = LumenIN%MolarFlow%NaCl
  ShellOUT%MolarFlow%NaCl = ShellIN%MolarFlow%NaCl
  LumenOUT%W = LumenIN%W+COM_MOD%AM*COM_MOD%Performance%JM
  ShellOUT%W = ShellIN%W-COM_MOD%AM*COM_MOD%Performance%JM
  LumenOUT%P = LumenIN%P
  ShellOUT%P = ShellIN%P  

!  write(12, "(A, 2F7.2)") "Outlet temperatures: ", LumenOUT%T, ShellOUT%T
!  write(12, "(A, E12.4, A)") "Permeation molar flux: ", COM_MOD%AM*DAVG(Fw)/b, " [kmol/m2-s]"
!  write(12, "(A, E12.4)") "Water molar flow of lumen-side influent: ", LumenIN%MolarFlow%H2O
!  write(12, "(A, E12.4)") "Water molar flow of shell-side influent: ", ShellIN%MolarFlow%H2O
!  write(12, "(A, 4E12.4)") "NaCl molar flow of lumen-side influent: ", LumenIN%MolarFlow%NaCl, LumenIN%MolarFlow%Na, LumenIN%MolarFlow%NaClS, LumenIN%MolarFlow%Cl
!  write(12, "(A, 4E12.4)") "NaCl molar flow of shell-side influent: ", ShellIN%MolarFlow%NaCl, ShellIN%MolarFlow%Na, ShellIN%MolarFlow%NaClS, ShellIN%MolarFlow%Cl
!  close(12)

  contains
  
  real(8) function G(xi) ! Eq.(A.26)
    real*8, intent(in) :: xi
    G = xi**four/four+(two*sigma**two-one)*xi**two+sigma**two*(four*sigma**two*dlog(sigma)-three*sigma**two-two*(xi**two-one))*dlog(xi)
  end function
  
end subroutine