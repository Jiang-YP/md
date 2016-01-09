module TemperatureDifference
    use CommonDef
    contains
    
    !calculate the axial temperature difference along the lumen center,see equation (59)
    real(8) function LumenTempDiff(Ta,Tb,Fw,Sq)
        real*8,intent(in) :: Ta,Tb,Sq,Fw
        call LocalFlux(Ta, Tb, Fw, Sq)
        LumenTempDiff = FOUR * Sq/(WaterConductivity(Ta) * LumenPeclet(Ta))
    end function
    
    
    !calculate the axial temperature difference along the cell surface,see equation (60)
    real(8) function ShellTempDiff(Ta,Tb,Fw,Sq)
        real*8,intent(in) ::Ta,Tb,Sq,Fw
        call LocalFlux(Ta, Tb, Fw, Sq)
        ShellTempDiff = FOUR * Sq/(WaterConductivity(Tb) * ShellPeclet(Tb))
    end function
    
    !see equation (12)
    real(8) function LumenPeclet(T)
        !thermalDiffusivity is the thermal diffusivity of the lumen stream water
        !a is the inner radius      
        !T is the water temperature
        real*8,intent(in) ::  T
        real*8 :: a
        a = COM_MOD%ID1/two
        LumenPeclet = TWO * COM_STREAM%u * a * a / (ThermalDiffusivity(T) * COM_MOD%LEN)
    end function
    
    !see equation (24)
    real(8) function ShellPeclet(T)
        !b is outer radius of a hollow fiber membrane
        !v is axial velocity in the outer channel
        !f represents the effect of packing fraction on the flow field, see equation(20)                                                              
        real*8,intent(in) :: T
        real*8 :: b,v,f
        b = COM_MOD%OD1/two
        v = COM_STREAM%u * COM_MOD%EPSILON /(ONE - COM_MOD%EPSILON)
        f = COM_MOD%EPSILON**TWO / (FOUR * COM_MOD%EPSILON - THREE - TWO * DLOG(COM_MOD%EPSILON) - COM_MOD%EPSILON**TWO)
        ShellPeclet = TWO * v  * b *b * f * (ONE/COM_MOD%EPSILON - ONE)/(ThermalDiffusivity(T) * COM_MOD%LEN)
    end function
    
    real(8) function ThermalDiffusivity(T)
        real*8,intent(in) :: T
        ThermalDiffusivity = WaterConductivity(T) / (COM_STREAM%PhysProp%rho * COM_STREAM%PhysProp%cp)
    end function 
    
    !see Appendix
    real(8) function WaterConductivity(T)
        real*8,intent(in) :: T
        WaterConductivity = 6.065d-1 * (-1.48445 + 4.12292*T/298.15 -1.63866*(T/298.15)**TWO)
    end function
end module