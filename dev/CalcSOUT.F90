subroutine CalcSOUT(ncomp, nmate, influents)
  use CommonDef
  use VaporConc
  integer, intent(in) :: ncomp, nmate
  real*8, dimension(ncomp,nmate), intent(in) :: influents
  integer :: i, j

!  open(13, file="tmp_CalcSOUT.log")
!  write(13, *) "Invoking CalcSOUT()"
!  do i = 1, ncomp
!    write(13, *) (influents(i,j), j = 1, nmate)
!  end do

! Set the COM_SIN  
  do i = 1, 2
!   Molar flow
    COM_SIN(i)%MolarFlow%H2O   = influents(1,i)
    COM_SIN(i)%MolarFlow%NaCl  = influents(2,i)
    COM_SIN(i)%MolarFlow%Na    = influents(3,i)
    COM_SIN(i)%MolarFlow%NaClS = influents(4,i)
    COM_SIN(i)%MolarFlow%Cl    = influents(5,i)
  end do

! Complete the rest parameters in global variable of COM_MOD 
  call CalcModule(COM_MOD)
! Complete the rest parameters in global stream variables 
  call CalcStream(COM_SIN(1), COM_MOD%CSA1)
  call CalcStream(COM_SIN(2), COM_MOD%CSA2)

  if (COM_OPT(3) .eq. 1) then
  ! Write MD module parameters into data file for debugging
    call WriteMOD(mod_def_file)
  ! Write streams into data file for debugging
    call WriteStream(stm_def_file)
  end if
  
! Calculate the temperature profiles in both lumen and shell sides  
  call CalcProfile(COM_SIN(1), COM_SIN(2), COM_SOUT(1), COM_SOUT(2))

!  write(13, "(A5, A12, A7)") "ISIDE", "Mass flow", "Temp."
!  write(13, "(I5, E12.4, F7.2)") ((I, COM_SIN(i)%W, COM_SIN(i)%T), i = 1, 2)
!  write(13, "(I5, E12.4, F7.2)") ((I, COM_SOUT(i)%W, COM_SOUT(i)%T), i = 1, 2)
!  close(13)  

end subroutine

subroutine WriteMOD(DataFileName)
  use CommonDef
  use toolkits
  
  character*(*), intent(in) :: DataFileName
  character*30, dimension(21) :: mod_data_name
  real, dimension(21) :: mod_data
  integer :: i
  
  mod_data(1) = COM_MOD%NUM
  mod_data(2) = COM_MOD%LEN
  mod_data(3) = COM_MOD%ID1
  mod_data(4) = COM_MOD%OD1
  mod_data(5) = COM_MOD%ID2
  mod_data(6) = COM_MOD%OD2
  mod_data(7) = COM_MOD%THK
  mod_data(8) = COM_MOD%PHI
  mod_data(9) = COM_MOD%AM
  mod_data(10) = COM_MOD%CSA1
  mod_data(11) = COM_MOD%CSA2
  mod_data(12) = COM_MOD%Membrane%rho
  mod_data(13) = COM_MOD%Membrane%mu
  mod_data(14) = COM_MOD%Membrane%cp
  mod_data(15) = COM_MOD%Membrane%KM
  mod_data(16) = COM_MOD%Membrane%CM
  mod_data(17) = COM_MOD%Membrane%porosity
  mod_data(18) = COM_MOD%Membrane%tortuosity
  mod_data(19) = COM_MOD%Membrane%PoreRadius
  mod_data(20) = COM_MOD%Performance%JM
  mod_data(21) = COM_MOD%Performance%eta

  mod_data_name(1) = 'COM_MOD%NUM'
  mod_data_name(2) = 'COM_MOD%LEN'
  mod_data_name(3) = 'COM_MOD%ID1'
  mod_data_name(4) = 'COM_MOD%OD1'
  mod_data_name(5) = 'COM_MOD%ID2'
  mod_data_name(6) = 'COM_MOD%OD2'
  mod_data_name(7) = 'COM_MOD%THK'
  mod_data_name(8) = 'COM_MOD%PHI'
  mod_data_name(9) = 'COM_MOD%AM'
  mod_data_name(10) = 'COM_MOD%CSA1'
  mod_data_name(11) = 'COM_MOD%CSA2'
  mod_data_name(12) = 'COM_MOD%Membrane%rho'
  mod_data_name(13) = 'COM_MOD%Membrane%mu'
  mod_data_name(14) = 'COM_MOD%Membrane%cp'
  mod_data_name(15) = 'COM_MOD%Membrane%KM'
  mod_data_name(16) = 'COM_MOD%Membrane%CM'
  mod_data_name(17) = 'COM_MOD%Membrane%porosity'
  mod_data_name(18) = 'COM_MOD%Membrane%tortuosity'
  mod_data_name(19) = 'COM_MOD%Membrane%PoreRadius'
  mod_data_name(20) = 'COM_MOD%Performance%JM'
  mod_data_name(21) = 'COM_MOD%Performance%eta'
  
  open(11, file = DataFileName, status = 'replace')
  write(11, mod_def_file_fmt) (mod_data_name(i), mod_data(i), i = 1, 21)
  write(11, *)
  write(11, *) 'Successfully generate data file at '//clock()//' '//date()
  close(11)
  
end subroutine

subroutine WriteStream(DataFileName)
  use CommonDef
  use toolkits
  
  character*(*), intent(in) :: DataFileName
  real, dimension(15,2) :: stream_data
  character*30, dimension(15) :: stream_data_name
  integer :: ISIDE
  
  do ISIDE = 1, 2
    stream_data(1,ISIDE) = COM_SIN(ISIDE)%W
    stream_data(2,ISIDE) = COM_SIN(ISIDE)%T
    stream_data(3,ISIDE) = COM_SIN(ISIDE)%P
    stream_data(4,ISIDE) = COM_SIN(ISIDE)%G
    stream_data(5,ISIDE) = COM_SIN(ISIDE)%u
    stream_data(6,ISIDE) = COM_SIN(ISIDE)%Re
    stream_data(7,ISIDE) = COM_SIN(ISIDE)%Nu
    stream_data(8,ISIDE) = COM_SIN(ISIDE)%PhysProp%rho
    stream_data(9,ISIDE) = COM_SIN(ISIDE)%PhysProp%mu
    stream_data(10,ISIDE) = COM_SIN(ISIDE)%PhysProp%cp
    stream_data(11,ISIDE) = COM_SIN(ISIDE)%PhysProp%KM
    stream_data(12,ISIDE) = COM_SIN(ISIDE)%MassFrac%H2O
    stream_data(13,ISIDE) = COM_SIN(ISIDE)%MassFrac%NaCl
    stream_data(14,ISIDE) = COM_SIN(ISIDE)%MolarFlow%H2O
    stream_data(15,ISIDE) = COM_SIN(ISIDE)%MolarFlow%NaCl    
  end do
    stream_data_name(1) = 'COM_SIN(ISIDE)%W'
    stream_data_name(2) = 'COM_SIN(ISIDE)%T'
    stream_data_name(3) = 'COM_SIN(ISIDE)%P'
    stream_data_name(4) = 'COM_SIN(ISIDE)%G'
    stream_data_name(5) = 'COM_SIN(ISIDE)%u'
    stream_data_name(6) = 'COM_SIN(ISIDE)%Re'
    stream_data_name(7) = 'COM_SIN(ISIDE)%Nu'
    stream_data_name(8) = 'COM_SIN(ISIDE)%PhysProp%rho'
    stream_data_name(9) = 'COM_SIN(ISIDE)%PhysProp%mu'
    stream_data_name(10) = 'COM_SIN(ISIDE)%PhysProp%cp'
    stream_data_name(11) = 'COM_SIN(ISIDE)%PhysProp%KM'
    stream_data_name(12) = 'COM_SIN(ISIDE)%MassFrac%H2O'
    stream_data_name(13) = 'COM_SIN(ISIDE)%MassFrac%NaCl'
    stream_data_name(14) = 'COM_SIN(ISIDE)%MolarFlow%H2O'
    stream_data_name(15) = 'COM_SIN(ISIDE)%MolarFlow%NaCl'
  
  open(11, file = DataFileName, status = 'replace')
  do i = 1, 15
    write(11,  stm_def_file_fmt) stream_data_name(i), stream_data(i,1), stream_data(i,2)
  end do
  write(11, *)
  write(11, *) 'Successfully generate data file at '//clock()//' '//date()  
  close(11)
  
end subroutine