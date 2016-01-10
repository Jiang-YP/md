module toolkits
use ifport

contains
      
! Export 2d double array in text file
  subroutine Export2DArray(DArray, DataFileName)
  real*8, dimension(:,:), intent(in) :: DArray
  character*(*), intent(in) :: DataFileName
  integer, dimension(2) :: DArrayShape
  integer :: i, j
! The DArray must have no more than 999 columns, else an error mark writes in data file
  character*3 :: ColNum_str
  character*10 :: fmt_str
  open(11, file = DataFileName, status='replace')  
  DArrayShape = shape(DArray)
  if (DArrayShape(2) .LE. 999) then
    write(ColNum_str, '(i3)') DArrayShape(2)
    fmt_str = '('//trim(ColNum_str)//'E12.5)'
    do i = 1, DArrayShape(1)
      write(11, fmt_str) (DArray(i,j), j = 1, DArrayShape(2))
    end do
  else
    write(1, *) 'The given array has more than 999 columns.'
  end if
  write(11, *)
  write(11, *) 'Successfully generate data file at '//clock()//' '//date()
  close(11)
  end subroutine
  
! Profile exporter
  subroutine ProfileExport(XLoc, YLoc, ZVal, DataFileName)
  real*8, dimension(:), intent(in) :: XLoc, YLoc
  real*8, dimension(:,:), intent(in) :: ZVal
  character*(*), intent(in) :: DataFileName
  integer, dimension(2) :: ZVal_shape
  integer :: NumPtX, NumPtY, i, j
  open(11, file = DataFileName, status = 'replace')
  NumPtX = size(XLoc)
  NumPtY = size(YLoc)
  ZVal_shape = shape(ZVal)
! Check the size of location arrays
  if ((NumPtX .eq. ZVal_shape(1)) .and. (NumPtY .eq. ZVal_shape(2))) then
    do i = 1, NumPtX
      do j = 1, NumPtY
        write(11, '(3e12.5)') XLoc(i), YLoc(j), ZVal(i,j)
      end do
    end do
    write(11, *)
    write(11, *) 'Successfully generate data file at '//clock()//' '//date()
  else
    write(11, *) 'The specified location arrays are not consistent with the given profile data!'
  end if
  close(11) 
  end subroutine
  
end module