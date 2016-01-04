subroutine LocalRate(Ta, Tb, Fw, Sq)
! Calculate local mass and heat transfer rates per unit membrane length, 
! referred to Eq.(45) and (46) in A S Kim, JMS, 455(2014) 168-186
    use VaporConc
    use CommonDef
    
    real*8,intent(in) :: Ta, Tb
    real*8, intent(out) :: Fw, Sq
    real*8, external :: f1, f2
    real*8 :: IntegralValue, AbsErr, ResAbs, ResAsc
    real*8 :: a, b

    ! Get the values from global variables
    a = COM_MOD%ID1/two
    b = COM_MOD%OD1/two

    ! Get the integral value of f1
    call dqk15(f1, Tb, Ta, IntegralValue, AbsErr, ResAbs, ResAsc)
    ! Calculate the local mass flux
    Fw = ONE/DLOG(b/a)*IntegralValue
    ! Get the integral value of f2
    call dqk15(f2, Tb, Ta, IntegralValue, AbsErr, ResAbs, ResAsc)
    ! Calculate the local mass flux
    Sq = ONE/DLOG(b/a)*IntegralValue
    
end subroutine
         
real*8 function f1(T)
    use VaporConc
    use CommonDef    
    real*8, intent(in) :: T
    real*8 :: epsilon, tau
    epsilon = COM_MOD%Membrane%porosity
    tau = tortuosity(epsilon)
    Deff = EffDiffusivity(T)
    f1 = Deff*epsilon/tau*dnw(T)
end function

real*8 function f2(T)
    use VaporConc
    use CommonDef    
    real*8, intent(in) :: T
    real*8 :: kappa, H
    real*8, external :: f1
    kappa = EffThermalCond(T)
    H = EffMolEnthalpy(T)
    f2 = kappa+H*f1(T)
end function
