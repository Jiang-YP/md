subroutine UpdateTemp(Tlmn,Tshl,Sq,Fw)
    use LocalTemp
    use CommonDef
    
    real*8,intent(in) : : Tlmn,Tshl
    real*8,intent(out) :: Sq,Fw
    real*8 :: Ta,Tb,z,SqSum,FwSum,SqNew,FwNew,SqPre,FwPre
    Ta = Tlmn
    Tb = Tshl
    SqSum = ZERO
    FwSum = ZERO
    call LocalFlux(Ta, Tb, Fw, Sq)
    SqSum = Sq
    FwPre = Fw
    do 
        z = ZERO
        SqPre = SqSum
        FwPre = FwSum
        SqSum = ZERO
        FwSum = ZERO
        do 
            Ta = InnerTemp(Ta,Tb,SqPre,FwPre,Tlmn,COM_MOD%ID1/two,z)
            Tb = OuterTemp(Ta,Tb,SqPre,FwPre,Tshl,COM_MOD%OD1/two,z)
            call LocalFlux(Ta, Tb, Fw, Sq)
            SqSum = SqSum + Sq
            FwSum = FwSum + Fw
            if (ABS(z - ONE)<1.0e-4) exit
            z = z + 0.05
        end do
        if (ABS(SqSum - SqPre)<1.0e-6 .and. ABS(FwSum - FwPre)<1.0e-6) exit
    end do
    Sq = SqSum
    Fw = FwSum
end subroutine UpdateTemp