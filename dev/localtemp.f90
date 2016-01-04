module LocalTemp
    use CommonDef
    use TemperatureDifference
    contains
    
    !calculate the inner temperature,see equation (61)
    real(8) function InnerTemp(Ta,Tb,Sq,Fw,Tlmn,r,y,z)
    !y is the dimensionless radial distance in lumen space y = r/a ( 0<y<1)
    !z is the dimensionless axial distance z = z/L
    !F see A.10
    real*8,intent(in) :: Ta,Tb,Tlmn,Sq,Fw,r,z
    real*8 :: y,F
    y = r/COM_MOD%ID1/two
    F = y*y - ONE/FOUR * y**FOUR
    InnerTemp = Tlmn - LumenTempDiff(Ta,Tb,Fw,Sq) * (z + ONE/FOUR * LumenPeclet(Tlmn) * F)
    end function
    
    !calculate the outer temperature,see equation (62)
    real(8) function OuterTemp(Ta,Tb,Sq,Fw,Tshl,x,r,z)
    !r is radial distance in shell space
    !x is the dimensionless radial distance in shell space x = r/b(1<x<¦Ò£©
    !z is the dimensionless axial distance z = z/L
    !o is the dimensionless radial distance of the cell surface o = c/b
    real*8,intent(in) :: Ta,Tb, Tshl,Sq,Fw,r,z
    real*8 :: x,o,H,Gx,Go1,go2
    x = r/COM_MOD%OD1/TWO
    o = (COM_MOD%ID2/TWO)/(COM_MOD%OD1/TWO)
    !see A.26
    Gx = ONE/FOUR*x**FOUR + (TWO*o*o-ONE)*x*x + o*o*(FOUR*o*o*DLOG(o)-THREE*o*o-TWO*(x*x-ONE))*DLOG(x)
    Go1 = ONE/FOUR*o**FOUR + (TWO*o*o-ONE)*o*o + o*o*(FOUR*o*o*DLOG(o)-THREE*o*o-TWO*(o*o-ONE))*DLOG(o)
    !see A.27
    go2 = o*o*(TWO*o*o*DLOG(o*o)-THREE*o*o+FOUR)-ONE
    H = (Go1 - Gx)/go2 !correct the previous equation,see equation (65)
    OuterTemp = Tshl + ShellTempDiff(Ta,Tb,Fw,Sq) * (ONE - z + ONE/FOUR * ShellPeclet(Tshl) * H)
    end function
end module