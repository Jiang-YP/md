!=========================================================================
!
! Author: Christian W. Straka
! e-mail: cstraka@ita.uni-heidelberg.de
!
! Institut fuer Theoretische Astrophysik    
! Tiergartenstr. 15
! 69121 Heidelberg
!
!==========================================================================
!
! DESCRIPTON: AUTOMATIC DIFFERENTIATION
!
!--------------------------------------------------------------------------

!============== MAY BE ALTERED BY USER ====================================

!definitions for precicion
module mod_precision
  !note: spr must not equal dpr
  integer, parameter :: spr = KIND(1.0)
  integer, parameter :: dpr = KIND(1.D0)
  integer, parameter :: ipr = KIND(1)
end module mod_precision

!definitions used by ADF95
module mod_ADF95defs
  use mod_precision
  integer(ipr), parameter :: LDsize = 20
end module mod_ADF95defs


!============== DO NOT ALTER BELOW HERE ===================================

!data structures and overloaded operators
module mod_ADF95types

  use mod_precision
  use mod_ADF95defs, only : LDsize
  
  !core data structure: double precision
  type ADF95_dpr
     real   (dpr)                      :: value  = 0.0_dpr
     real   (dpr), dimension(1:LDsize) :: deriv  = 0.0_dpr
     integer(ipr), dimension(0:LDsize) :: index  = 0_ipr
  end type ADF95_dpr

  !operators: assignment
  interface assignment(=)
     module procedure assignment1, assignment2, assignment3, assignment4 
  end interface

  !operators: elementary math
  interface operator(**)
     module procedure power1, power2, power3, power4, power5, &
                    & power6, power7
  end interface
  interface operator(*)
     module procedure multiply1, multiply2, multiply3, multiply4, &
                    & multiply5, multiply6, multiply7
  end interface
  interface operator(/)
     module procedure divide1, divide2, divide3, divide4, divide5, &
                    & divide6, divide7
  end interface
  interface operator(+)
     module procedure add1, add2, add3, add4, add5, add6, add7, unary_plus1
  end interface
  interface operator(-)
     module procedure minus1, minus2, minus3, minus4, minus5,  &
                    & minus6, minus7, unary_minus1
  end interface

  !operators: comparison
  interface operator(.gt.)
     module procedure gt1, gt2, gt3, gt4, gt5, gt6, gt7
  end interface
  interface operator(.ge.)
     module procedure ge1, ge2, ge3, ge4, ge5, ge6, ge7
  end interface
  interface operator(.eq.)
     module procedure eq1, eq2, eq3, eq4, eq5, eq6, eq7
  end interface
  interface operator(.ne.)
     module procedure ne1, ne2, ne3, ne4, ne5, ne6, ne7
  end interface
  interface operator(.le.)
     module procedure le1, le2, le3, le4, le5, le6, le7
  end interface
  interface operator(.lt.)
     module procedure lt1, lt2, lt3, lt4, lt5, lt6, lt7
  end interface

  !operators: elemental mathematical functions
  interface log
     module procedure log_e
  end interface
  interface log10
     module procedure log_10
  end interface
  interface exp
     module procedure exponent1
  end interface
  interface sqrt
     module procedure square_root
  end interface
  interface sin
     module procedure sinus
  end interface
  interface cos
     module procedure cosinus
  end interface
  interface tan
     module procedure tangens
  end interface
  interface sinh
     module procedure sinus_hyperbolicus
  end interface
  interface cosh
     module procedure cosinus_hyperbolicus
  end interface
  interface tanh
     module procedure tangens_hyperbolicus
  end interface
  interface asin
     module procedure arcus_sinus
  end interface
  interface acos
     module procedure arcus_cosinus
  end interface
  interface atan
     module procedure arcus_tangens
  end interface
  interface atan2
     module procedure arcus2_tangens1, arcus2_tangens2, arcus2_tangens3
  end interface

  !operators: inquiry functions for any type
  interface kind
     module procedure kind1
  end interface

  !operators: transformal functions that reduce arrays
  interface sum
     module procedure sum_vec
  end interface
  interface product
     module procedure product_vec
  end interface
  interface minval
     module procedure minimum_value
  end interface
  interface minloc
     module procedure minimum_loc1
  end interface
  interface maxval
     module procedure maximum_value
  end interface
  interface maxloc
     module procedure maximal_loc1
  end interface

  !operators: elemental functions that do not convert
  interface dim
     module procedure dim1
  end interface
  interface max
     module procedure max_value1, max_value2, max_value3
  end interface
  interface min
     module procedure min_value1, min_value2, min_value3
  end interface
  interface mod
     module procedure mod1, mod2, mod3
  end interface
  interface modulo
     module procedure modulo1, modulo2, modulo3
  end interface
  interface sign
     module procedure signum1, signum2, signum3
  end interface

  !operator: elemental functions that may convert
  interface abs
     module procedure absolute_value1
  end interface
  interface aint
     module procedure aint1
  end interface
  interface anint
     module procedure anint1
  end interface
  interface ceiling
     module procedure ceiling1
  end interface
  interface floor
     module procedure floor1
  end interface
  interface int
     module procedure int1
  end interface
  interface nint
     module procedure nint1
  end interface


  !operator: numerical inquiry functions
  interface digits
     module procedure digits1
  end interface
  interface epsilon
     module procedure epsilon1
  end interface
  interface huge
     module procedure huge1
  end interface
  interface maxexponent
     module procedure maxexponent1
  end interface
  interface minexponent
     module procedure minexponent1
  end interface
  interface precision
     module procedure precision1
  end interface
  interface radix
     module procedure radix1
  end interface
  interface range
     module procedure range1
  end interface
  interface tiny
     module procedure tiny1
  end interface

  !operator: elemental functions to manipulate reals
  interface exponent
     module procedure exponent_iqr1
  end interface
  interface fraction
     module procedure fraction1
  end interface
  interface nearest
     module procedure nearest1, nearest2, nearest3
  end interface
  interface rrspacing
     module procedure rrspacing1
  end interface
  interface scale
     module procedure scale1
  end interface
  interface set_exponent
     module procedure set_exponent1
  end interface
  interface spacing
     module procedure spacing1
  end interface

  !operator: vector and matrix multipliacation functions
  !matmul does not work with Absoft (Version 8.2a): BEGIN comment out
  interface dot_product
     module procedure dot_product1
  end interface
  interface matmul
     module procedure matmul1, matmul2, matmul3
  end interface
  !matmul does not work with Absoft (Version 8.2a): END comment out

  contains

    !operators: assignment

    !assignment1: (ADF95_dpr, int(ipr))
    elemental subroutine assignment1(a,b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(out) :: a
      integer(ipr)   , intent(in)  :: b

      a%value    = real(b,dpr)
      a%index(0) = 0
    end subroutine assignment1

    !assignment2: (ADF95_dpr, dpr)
    elemental subroutine assignment2(a,b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(out) :: a
      real(dpr)      , intent(in)  :: b

      a%value    = b
      a%index(0) = 0
    end subroutine assignment2

    !assignment3: (ADF95_dpr, spr)
    elemental subroutine assignment3(a,b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(out) :: a
      real(spr)      , intent(in)  :: b

      a%value    = real(b,dpr)
      a%index(0) = 0
    end subroutine assignment3

    !assignment4: (ADF95_dpr, ADF95_dpr)
    elemental subroutine assignment4(a,b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(out) :: a
      type(ADF95_dpr), intent(in)  :: b
      integer(ipr)                 :: lenb

      lenb = b%index(0)

      a%value = b%value

      a%deriv(1:lenb) = b%deriv(1:lenb)
      a%index(0:lenb) = b%index(0:lenb)
    end subroutine assignment4

    !operators: elementary math

    !OPERATOR(**)

    !power1: (int(ipr),ADF95_dpr)
    elemental function power1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in)  :: a
      type(ADF95_dpr), intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lenb

      f%value = a ** b%value

      if(a .ne. 0) then
         lenb = b%index(0)
         f%deriv(1:lenb) = f%value * log(real(a,dpr)) &
              & * b%deriv(1:lenb)
         f%index(0:lenb) = b%index(0:lenb)
      endif
    end function power1

    !power2: (ADF95_dpr,int(ipr))
    elemental function power2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in)  :: a
      integer(ipr)   , intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lena

      lena = a%index(0)

      f%value = a%value ** b

      f%deriv(1:lena) = b * a%value ** (b-1) * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function power2

    !power3: (real(dpr),ADF95_dpr)
    elemental function power3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in)  :: a
      type(ADF95_dpr), intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lenb
      
      f%value = a ** b%value

      if(a .ne. 0.0_dpr) then
         lenb = b%index(0)
         f%deriv(1:lenb) = f%value * log(a) * b%deriv(1:lenb)
         f%index(0:lenb) = b%index(0:lenb)
      endif
    end function power3

    !power4: (ADF95_dpr, real(dpr))
    elemental function power4(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in)  :: a
      real(dpr)      , intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lena

      lena = a%index(0)
      
      f%value = a%value ** b

      f%deriv(1:lena) = b * a%value ** (b-1.0_dpr)*a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function power4

    !power5: (ADF95_dpr,ADF95_dpr)
    elemental function power5(a, b) result(f)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), intent(in)       :: a, b
      type(ADF95_dpr)                   :: f
      integer(ipr), dimension(1:LDsize) :: fidxa, fidxb
      integer(ipr)                      :: j, k, l, lena, lenb

      f%value = a%value ** b%value

      if(a%value .ne. 0.0_dpr) then      

         lena = a%index(0)
         lenb = b%index(0)

         j=1; k=1; l=0
         do
            if((a%index(j) .le. b%index(k) .and. &
                 & a%index(j) .ne. 0) .or. (b%index(k) .eq. 0)) then
               l = l + 1
               f%index(l) = a%index(j)
               fidxa(j) = l
               if(a%index(j) .eq. b%index(k)) then
                  fidxb(k) = l
                  k = k + 1
                  if(k .gt. LDsize) exit 
               endif
               j = j + 1
               if(j .gt. LDsize) exit 
            else
               l = l + 1
               f%index(l)  = b%index(k)
               fidxb(k) = l
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            if(j .gt. lena .and. k .gt. lenb) exit
            if(l .ge. LDsize) exit 
         enddo
         f%index(0) = l

         do j=1, lena
            f%deriv(fidxa(j)) = f%value *b%value/a%value*a%deriv(j)
         enddo

         do j=1, lenb
            f%deriv(fidxb(j)) = f%deriv(fidxb(j)) + &
                 & f%value*log(a%value)*b%deriv(j)
         enddo

      endif
    end function power5

    !power6: (real(spr),ADF95_dpr)
    elemental function power6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in)  :: a
      type(ADF95_dpr), intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lenb
      
      f%value = real(a,dpr) ** b%value

      if(a .ne. 0.0_spr) then
         lenb = b%index(0)
         f%deriv(1:lenb) = f%value * real(log(a),dpr) * b%deriv(1:lenb)
         f%index(0:lenb) = b%index(0:lenb)
      endif
    end function power6

    !power7: (ADF95_dpr, real(spr))
    elemental function power7(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in)  :: a
      real(spr)      , intent(in)  :: b
      type(ADF95_dpr)              :: f
      integer(ipr)                 :: lena

      lena = a%index(0)
      
      f%value = a%value ** real(b,dpr)

      f%deriv(1:lena) = real(b,dpr) &
           & * a%value ** (real(b,dpr)-1.0_dpr)*a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function power7

    !OPERATOR(*)

    !multiply1: (int(ipr),ADF95_dpr)
    elemental function multiply1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in)    :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = a * b%value 

      f%deriv(1:lenb) = a * b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function multiply1

    !multiply2: (ADF95_dpr,int(ipr))
    elemental function multiply2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      integer(ipr)   , intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value * b

      f%deriv(1:lena) = a%deriv(1:lena) * b
      f%index(0:lena) = a%index(0:lena)
    end function multiply2

    !multiply3: (real(dpr),ADF95_dpr)
    elemental function multiply3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in)    :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = a * b%value 

      f%deriv(1:lenb) = a * b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function multiply3

    !multiply4: (ADF95_dpr,real(dpr))
    elemental function multiply4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      real(dpr)       , intent(in)   :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value * b

      f%deriv(1:lena) = a%deriv(1:lena) * b
      f%index(0:lena) = a%index(0:lena)
    end function multiply4

    !multiply5: (ADF95_dpr,ADF95_dpr)
    elemental function multiply5(a, b) result(f)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), intent(in)       :: a, b
      type(ADF95_dpr)                   :: f
      integer(ipr), dimension(1:LDsize) :: fidxa, fidxb
      integer(ipr)                      :: j, k, l, lena, lenb

      lena = a%index(0)
      lenb = b%index(0)

      f%value = a%value * b%value

      j=1; k=1; l=0
      do
         if((a%index(j) .le. b%index(k) .and. &
              & a%index(j) .ne. 0) .or. (b%index(k) .eq. 0)) then
            l = l + 1
            f%index(l) = a%index(j)
            fidxa(j) = l
            if(a%index(j) .eq. b%index(k)) then
               fidxb(k) = l
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            j = j + 1
            if(j .gt. LDsize) exit 
         else
            l = l + 1
            f%index(l)  = b%index(k)
            fidxb(k) = l
            k = k + 1
            if(k .gt. LDsize) exit 
         endif
         if(j .gt. lena .and. k .gt. lenb) exit
         if(l .ge. LDsize) exit 
      enddo
      f%index(0) = l

      do j=1, lena
         f%deriv(fidxa(j)) = b%value * a%deriv(j)
      enddo

      do j=1, lenb
         f%deriv(fidxb(j)) = f%deriv(fidxb(j)) + a%value * b%deriv(j)
      enddo
    end function multiply5

    !multiply6: (real(spr),ADF95_dpr)
    elemental function multiply6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in)    :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = real(a,dpr) * b%value 

      f%deriv(1:lenb) = real(a,dpr) * b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function multiply6

    !multiply7: (ADF95_dpr,real(spr))
    elemental function multiply7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      real(spr)       , intent(in)   :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value * real(b,dpr)

      f%deriv(1:lena) = a%deriv(1:lena) * real(b,dpr)
      f%index(0:lena) = a%index(0:lena)
    end function multiply7

    !OPERATOR(/)

    !divide1: (int(ipr),ADF95_dpr)
    elemental function divide1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in)    :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = a / b%value 

      f%deriv(1:lenb) = -f%value / b%value * b%deriv(1:lenb) 
      f%index(0:lenb) = b%index(0:lenb)
    end function divide1

    !divide2: (ADF95_dpr,int(ipr))
    elemental function divide2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      integer(ipr)   , intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value / b

      f%deriv(1:lena) = a%deriv(1:lena) / b
      f%index(0:lena) = a%index(0:lena)
    end function divide2

    !divide3: (real(dpr),ADF95_dpr)
    elemental function divide3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)       , intent(in)   :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = a / b%value 

      f%deriv(1:lenb) = -f%value / b%value * b%deriv(1:lenb) 
      f%index(0:lenb) = b%index(0:lenb)
    end function divide3

    !divide4: (ADF95_dpr,real(dpr))
    elemental function divide4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      real(dpr)       , intent(in)   :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value / b

      f%deriv(1:lena) = a%deriv(1:lena) / b
      f%index(0:lena) = a%index(0:lena)
    end function divide4

    !divide5: (ADF95_dpr,ADF95_dpr)
    elemental function divide5(a, b) result(f)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), intent(in)       :: a, b
      type(ADF95_dpr)                   :: f
      integer(ipr), dimension(1:LDsize) :: fidxa, fidxb
      integer(ipr)                      :: j, k, l, lena, lenb

      lena = a%index(0)
      lenb = b%index(0)

      f%value = a%value / b%value
      
      j=1; k=1; l=0
      do
         if((a%index(j) .le. b%index(k) .and. &
              & a%index(j) .ne. 0) .or. (b%index(k) .eq. 0)) then
            l = l + 1
            f%index(l) = a%index(j)
            fidxa(j) = l
            if(a%index(j) .eq. b%index(k)) then
               fidxb(k) = l
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            j = j + 1
            if(j .gt. LDsize) exit 
         else
            l = l + 1
            f%index(l)  = b%index(k)
            fidxb(k) = l
            k = k + 1
            if(k .gt. LDsize) exit 
         endif
         if(j .gt. lena .and. k .gt. lenb) exit
         if(l .ge. LDsize) exit 
      enddo
      f%index(0) = l

      do j=1, lena
         f%deriv(fidxa(j)) = a%deriv(j) / b%value
      enddo

      do j=1, lenb
         f%deriv(fidxb(j)) = f%deriv(fidxb(j)) - b%deriv(j) * f%value/b%value
      enddo
    end function divide5

    !divide6: (real(spr),ADF95_dpr)
    elemental function divide6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)       , intent(in)   :: a
      type(ADF95_dpr), intent(in)    :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lenb

      lenb = b%index(0)

      f%value = real(a,dpr) / b%value 

      f%deriv(1:lenb) = -f%value / b%value * b%deriv(1:lenb) 
      f%index(0:lenb) = b%index(0:lenb)
    end function divide6

    !divide7: (ADF95_dpr,real(spr))
    elemental function divide7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)    :: a
      real(spr)       , intent(in)   :: b
      type(ADF95_dpr)                :: f
      integer(ipr)                   :: lena

      lena = a%index(0)
      
      f%value = a%value / real(b,dpr)

      f%deriv(1:lena) = a%deriv(1:lena) / real(b,dpr)
      f%index(0:lena) = a%index(0:lena)
    end function divide7

    !OPERATOR(+)

    !add1: (int(ipr),ADF95_dpr)
    elemental function add1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lenb

      lenb = b%index(0)

      f%value = a + b%value

      f%deriv(1:lenb) = b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function add1

    !add2: (ADF95_dpr,int(ipr))
    elemental function add2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value + b

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function add2

    !add3: (real(dpr),ADF95_dpr)
    elemental function add3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lenb

      lenb = b%index(0)

      f%value = a + b%value

      f%deriv(1:lenb) = b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function add3

    !add4: (ADF95_dpr,real(dpr))
    elemental function add4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value + b

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function add4

    !add5: (ADF95_dpr,ADF95_dpr)
    elemental function add5(a, b) result(f)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), intent(in)       :: a, b
      type(ADF95_dpr)                   :: f
      integer(ipr), dimension(1:LDsize) :: fidxa, fidxb
      integer(ipr)                      :: j, k, l, lena, lenb

      lena = a%index(0)
      lenb = b%index(0)

      f%value = a%value + b%value
      
      j=1; k=1; l=0
      do
         if((a%index(j) .le. b%index(k) .and. &
              & a%index(j) .ne. 0) .or. (b%index(k) .eq. 0)) then
            l = l + 1
            f%index(l) = a%index(j)
            fidxa(j) = l
            if(a%index(j) .eq. b%index(k)) then
               fidxb(k) = l
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            j = j + 1
            if(j .gt. LDsize) exit 
         else
            l = l + 1
            f%index(l)  = b%index(k)
            fidxb(k) = l
            k = k + 1
            if(k .gt. LDsize) exit 
         endif
         if(j .gt. lena .and. k .gt. lenb) exit
         if(l .ge. LDsize) exit 
      enddo
      f%index(0) = l

      do j=1, lena
         f%deriv(fidxa(j)) = a%deriv(j)
      enddo

      do j=1, lenb
         f%deriv(fidxb(j)) = f%deriv(fidxb(j)) + b%deriv(j)
      enddo
    end function add5

    !add6: (real(spr),ADF95_dpr)
    elemental function add6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lenb

      lenb = b%index(0)

      f%value = real(a,dpr) + b%value

      f%deriv(1:lenb) = b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function add6

    !add7: (ADF95_dpr,real(spr))
    elemental function add7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value + real(b,dpr)

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function add7

    !unary_plus1: (ADF95_dpr)
    elemental function unary_plus1(a) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function unary_plus1

    !unary_minus1: (ADF95_dpr)
    elemental function unary_minus1(a) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)
      
      f%value = - a%value

      f%deriv(1:lena) = - a%deriv(1:lena)
      f%index(0:lena) =   a%index(0:lena)
    end function unary_minus1

    !OPERATOR(-)

    !minus1: (int(ipr),ADF95_dpr)
    elemental function minus1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                      :: lenb

      lenb = b%index(0)

      f%value = a - b%value

      f%deriv(1:lenb) = - b%deriv(1:lenb)
      f%index(0:lenb) = b%index(0:lenb)
    end function minus1

    !minus2: (ADF95_dpr,int(ipr))
    elemental function minus2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)
      
      f%value = a%value - b

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function minus2

    !minus3: (real(dpr),ADF95_dpr)
    elemental function minus3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lenb

      lenb = b%index(0)

      f%value = a - b%value

      f%deriv(1:lenb) = - b%deriv(1:lenb)
      f%index(0:lenb) =   b%index(0:lenb)
    end function minus3

    !minus4: (ADF95_dpr,real(dpr))
    elemental function minus4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value - b

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function minus4

    !minus5: (ADF95_dpr,ADF95_dpr)
    elemental function minus5(a, b) result(f)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), intent(in)       :: a, b
      type(ADF95_dpr)                   :: f
      integer(ipr), dimension(1:LDsize) :: fidxa, fidxb
      integer(ipr)                      :: j, k, l, lena, lenb

      lena = a%index(0)
      lenb = b%index(0)

      f%value = a%value - b%value
      
      j=1; k=1; l=0
      do
         if((a%index(j) .le. b%index(k) .and. &
              & a%index(j) .ne. 0) .or. (b%index(k) .eq. 0)) then
            l = l + 1
            f%index(l) = a%index(j)
            fidxa(j) = l
            if(a%index(j) .eq. b%index(k)) then
               fidxb(k) = l
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            j = j + 1
            if(j .gt. LDsize) exit 
         else
            l = l + 1
            f%index(l)  = b%index(k)
            fidxb(k) = l
            k = k + 1
            if(k .gt. LDsize) exit 
         endif
         if(j .gt. lena .and. k .gt. lenb) exit
         if(l .ge. LDsize) exit 
      enddo
      f%index(0) = l

      do j=1, lena
         f%deriv(fidxa(j)) = a%deriv(j)
      enddo

      do j=1, lenb
         f%deriv(fidxb(j)) = f%deriv(fidxb(j)) - b%deriv(j)
      enddo
    end function minus5

    !minus6: (real(spr),ADF95_dpr)
    elemental function minus6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lenb

      lenb = b%index(0)

      f%value = real(a,dpr) - b%value

      f%deriv(1:lenb) = - b%deriv(1:lenb)
      f%index(0:lenb) =   b%index(0:lenb)
    end function minus6

    !minus7: (ADF95_dpr,real(spr))
    elemental function minus7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = a%value - real(b,dpr)

      f%deriv(1:lena) = a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function minus7

    !operators: comparison

    !OPERATOR(.gt.)

    !gt1: (int(ipr),ADF95_dpr)
    elemental function gt1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .gt. b%value
    end function gt1

    !gt2: (ADF95_dpr,int(ipr))
    elemental function gt2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .gt. b
    end function gt2

    !gt3: (real(dpr),ADF95_dpr)
    elemental function gt3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .gt. b%value
    end function gt3

    !gt4: (ADF95_dpr,real(dpr))
    elemental function gt4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .gt. b
    end function gt4

    !gt5: (ADF95_dpr,ADF95_dpr)
    elemental function gt5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .gt. b%value
    end function gt5

    !gt6: (real(spr),ADF95_dpr)
    elemental function gt6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .gt. b%value
    end function gt6

    !gt7: (ADF95_dpr,real(spr))
    elemental function gt7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .gt. real(b,dpr)
    end function gt7

    !OPERATOR(.ge.)

    !ge1: (int(ipr),ADF95_dpr)
    elemental function ge1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .ge. b%value
    end function ge1

    !ge2: (ADF95_dpr,int(ipr))
    elemental function ge2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .ge. b
    end function ge2

    !ge3: (real(dpr),ADF95_dpr)
    elemental function ge3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .ge. b%value
    end function ge3

    !ge4: (ADF95_dpr,real(dpr))
    elemental function ge4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .ge. b
    end function ge4

    !ge5: (ADF95_dpr,ADF95_dpr)
    elemental function ge5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .ge. b%value
    end function ge5

    !ge6: (real(spr),ADF95_dpr)
    elemental function ge6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .ge. b%value
    end function ge6

    !ge7: (ADF95_dpr,real(spr))
    elemental function ge7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .ge. real(b,dpr)
    end function ge7

    !OPERATOR(.eq.)

    !eq1: (int(ipr),ADF95_dpr)
    elemental function eq1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .eq. b%value
    end function eq1

    !eq2: (ADF95_dpr,int(ipr))
    elemental function eq2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .eq. b
    end function eq2

    !eq3: (real(dpr),ADF95_dpr)
    elemental function eq3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .eq. b%value
    end function eq3

    !eq4: (ADF95_dpr,real(dpr))
    elemental function eq4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .eq. b
    end function eq4

    !eq5: (ADF95_dpr,ADF95_dpr)
    elemental function eq5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .eq. b%value
    end function eq5

    !eq6: (real(spr),ADF95_dpr)
    elemental function eq6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .eq. b%value
    end function eq6

    !eq7: (ADF95_dpr,real(spr))
    elemental function eq7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .eq. real(b,dpr)
    end function eq7

    !OPERATOR(.ne.)

    !ne1: (int(ipr),ADF95_dpr)
    elemental function ne1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .ne. b%value
    end function ne1

    !ne2: (ADF95_dpr,int(ipr))
    elemental function ne2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .ne. b
    end function ne2

    !ne3: (real(dpr),ADF95_dpr)
    elemental function ne3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .ne. b%value
    end function ne3

    !ne4: (ADF95_dpr,real(dpr))
    elemental function ne4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .ne. b
    end function ne4

    !ne5: (ADF95_dpr,ADF95_dpr)
    elemental function ne5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .ne. b%value
    end function ne5

    !ne6: (real(spr),ADF95_dpr)
    elemental function ne6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .ne. b%value
    end function ne6

    !ne7: (ADF95_dpr,real(spr))
    elemental function ne7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .ne. real(b,dpr)
    end function ne7

    !OPERATOR(.le.)

    !le1: (int(ipr),ADF95_dpr)
    elemental function le1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .le. b%value
    end function le1

    !le2: (ADF95_dpr,int(ipr))
    elemental function le2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .le. b
    end function le2

    !le3: (real(dpr),ADF95_dpr)
    elemental function le3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .le. b%value
    end function le3

    !le4: (ADF95_dpr,real(dpr))
    elemental function le4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .le. b
    end function le4

    !le5: (ADF95_dpr,ADF95_dpr)
    elemental function le5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .le. b%value
    end function le5

    !le6: (real(spr),ADF95_dpr)
    elemental function le6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .le. b%value
    end function le6

    !le7: (ADF95_dpr,real(spr))
    elemental function le7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .le. real(b,dpr)
    end function le7

    !OPERATOR(.lt.)

    !lt1: (int(ipr),ADF95_dpr)
    elemental function lt1(a, b) result(f)
      implicit none
      integer(ipr)   , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .lt. b%value
    end function lt1

    !lt2: (ADF95_dpr,int(ipr))
    elemental function lt2(a, b) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: b
      logical                     :: f

      f = a%value .lt. b
    end function lt2

    !lt3: (real(dpr),ADF95_dpr)
    elemental function lt3(a, b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = a .lt. b%value
    end function lt3

    !lt4: (ADF95_dpr,real(dpr))
    elemental function lt4(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      logical                     :: f

      f = a%value .lt. b
    end function lt4

    !lt5: (ADF95_dpr,ADF95_dpr)
    elemental function lt5(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      logical                     :: f

      f = a%value .lt. b%value
    end function lt5

    !lt6: (real(spr),ADF95_dpr)
    elemental function lt6(a, b) result(f)
      use mod_precision
      implicit none
      real(spr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      logical                     :: f

      f = real(a,dpr) .lt. b%value
    end function lt6

    !lt7: (ADF95_dpr,real(spr))
    elemental function lt7(a, b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(spr)      , intent(in) :: b
      logical                     :: f

      f = a%value .lt. real(b,dpr)
    end function lt7

    !operators: elemental mathematical functions

    !log_e: (ADF95_dpr)
    elemental function log_e(a) result(f)
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)
      
      f%value = log(a%value)

      f%deriv(1:lena) = a%deriv(1:lena) / a%value
      f%index(0:lena) = a%index(0:lena)
    end function log_e

    !log_10: (ADF95_dpr)
    elemental function log_10(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)
      
      f%value = log10(a%value)

      f%deriv(1:lena) = a%deriv(1:lena) / (a%value * log(10.0_dpr))
      f%index(0:lena) = a%index(0:lena)
    end function log_10

    !exp: (ADF95_dpr)
    elemental function exponent1(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = exp(a%value)

      f%deriv(1:lena) = f%value * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function exponent1

    !sqrt: (ADF95_dpr)
    elemental function square_root(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      f%value = sqrt(a%value)

      !Note: derivative not defined at
      !      a%value .eq. 0.0_dpr: produces Inf (Infinity)
      lena = a%index(0)
      f%deriv(1:lena) = a%deriv(1:lena) / (2.0_dpr * f%value)
      f%index(0:lena) = a%index(0:lena)
    end function square_root

    !sin: (ADF95_dpr)
    elemental function sinus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = sin(a%value)

      f%deriv(1:lena) = cos(a%value) * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function sinus

    !cos: (ADF95_dpr)
    elemental function cosinus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = cos(a%value)

      f%deriv(1:lena) = -sin(a%value) * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function cosinus

    !tan: (ADF95_dpr)
    elemental function tangens(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = tan(a%value)

      f%deriv(1:lena) = a%deriv(1:lena)  / cos(a%value)**2
      f%index(0:lena) = a%index(0:lena)
    end function tangens

    !sin_hyperbolicus: (ADF95_dpr)
    elemental function sinus_hyperbolicus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = sinh(a%value)

      f%deriv(1:lena) = cosh(a%value) * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function sinus_hyperbolicus

    !cos_hyperbolicus: (ADF95_dpr)
    elemental function cosinus_hyperbolicus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = cosh(a%value)

      f%deriv(1:lena) = sinh(a%value) * a%deriv(1:lena)
      f%index(0:lena) = a%index(0:lena)
    end function cosinus_hyperbolicus

    !tan_hyperbolicus: (ADF95_dpr)
    elemental function tangens_hyperbolicus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = tanh(a%value)

      if(abs(a%value) .lt. 2.0_dpr*range(a%value)) then
         f%deriv(1:lena) = a%deriv(1:lena) / cosh(a%value)**2
      else
         f%deriv(1:lena) = a%deriv(1:lena) * 4.0_dpr*exp(-2.0_dpr*abs(a%value))
      endif
      f%index(0:lena) = a%index(0:lena)
    end function tangens_hyperbolicus

    !arcus_sinus: (ADF95_dpr)
    elemental function arcus_sinus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = asin(a%value)

      !Note: derivative not defined at
      !      abs(a%value) .eq. 1.0_dpr: produces Inf (Infinity)
      f%deriv(1:lena) = a%deriv(1:lena) / sqrt(1.0_dpr - a%value**2)
      f%index(0:lena) = a%index(0:lena)
    end function arcus_sinus

    !arcus_cosinus: (ADF95_dpr)
    elemental function arcus_cosinus(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = acos(a%value)

      !Note: derivative not defined at
      !      abs(a%value) .eq. 1.0_dpr: produces -Inf (Infinity)
      f%deriv(1:lena) = - a%deriv(1:lena) / sqrt(1.0_dpr - a%value**2)
      f%index(0:lena) =   a%index(0:lena)
    end function arcus_cosinus

    !arcus_tangens: (ADF95_dpr)
    elemental function arcus_tangens(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)

      f%value = atan(a%value)

      f%deriv(1:lena) = a%deriv(1:lena) / (1.0_dpr + a%value**2)
      f%index(0:lena) = a%index(0:lena)
    end function arcus_tangens

    !arcus2_tangens1: (ADF95_dpr, ADF95_dpr)
    elemental function arcus2_tangens1(a,b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      type(ADF95_dpr)             :: f

      !Note: derivative undefined for b=0
      !      produces -NaN (Not a Number)
      f       = atan(a/b)
      f%value = atan2(a%value,b%value)
    end function arcus2_tangens1

    !arcus2_tangens2: (ADF95_dpr, dpr)
    elemental function arcus2_tangens2(a,b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      type(ADF95_dpr)             :: f

      !Note: derivative undefined for b=0
      !      produces -NaN (Not a Number)
      f       = atan(a/b)
      f%value = atan2(a%value,b)
    end function arcus2_tangens2

    !arcus2_tangens3: (dpr, ADF95_dpr)
    elemental function arcus2_tangens3(a,b) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: f

      !Note: derivative undefined for b=0
      !      produces -NaN (Not a Number)
      f       = atan(a/b)
      f%value = atan2(a,b%value)
    end function arcus2_tangens3

    !operators: inquiry functions for any type
    
    !kind1: (ADF95_dpr)
    elemental function kind1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = kind(a%value)
    end function kind1

    !operators: transformal functions that reduce arrays

    !sum_vec: (ADF95_dpr)
    pure function sum_vec(vecin,msk) result(sumvec)
      use mod_precision
      use mod_ADF95defs, only : LDsize
      implicit none
      type(ADF95_dpr), dimension(:)          , intent(in):: vecin
      type(ADF95_dpr), dimension(:), allocatable         :: vec
      type(ADF95_dpr)                                    :: sumvec, f2
      integer(ipr), dimension(1:LDsize,1:size(vecin))    :: fidx
      logical     , dimension(1:size(vecin)), intent(in), optional :: msk
      integer(ipr)                                   :: nsize, i, j, k, l

      nsize = size(vecin)

      !Absoft (Version 8.2a) doesn't like it in argument list
      if(.not. allocated(vec)) allocate(vec(1:nsize))

      if(present(msk)) then
         where(msk(1:nsize))
            vec(1:nsize) = vecin(1:nsize)
         endwhere
      else
         vec(1:nsize) = vecin(1:nsize)
      endif

      sumvec%value = sum(vec(1:nsize)%value)
            
      sumvec%index(0:vec(1)%index(0)) = vec(1)%index(0:vec(1)%index(0))

      do i=2, nsize
         j=1; k=1; l=0
         do
            if((vec(i)%index(j) .le. sumvec%index(k) .and. &
                 & vec(i)%index(j) .ne. 0) .or. (sumvec%index(k) .eq. 0)) then
               l = l + 1
               f2%index(l) = vec(i)%index(j)
               if(vec(i)%index(j) .eq. sumvec%index(k)) then
                  k = k + 1
                  if(k .gt. LDsize) exit 
               endif
               j = j + 1
               if(j .gt. LDsize) exit 
            else
               l = l + 1
               f2%index(l)  = sumvec%index(k)
               k = k + 1
               if(k .gt. LDsize) exit 
            endif
            if(vec(i)%index(j) .eq. 0 .and. sumvec%index(k) .eq. 0) exit
            if(l .ge. LDsize) exit 
         enddo
         f2%index(0)       = l
         sumvec%index(0:l) = f2%index(0:l)
      enddo

      do i=1, nsize
         do k=1, sumvec%index(0)
            do j=1, vec(i)%index(0)
               if(sumvec%index(k) .eq. vec(i)%index(j)) fidx(j,i) = k
            enddo
         enddo
      enddo

      do i=1, nsize
         do j=1, vec(i)%index(0)
            sumvec%deriv(fidx(j,i)) = sumvec%deriv(fidx(j,i)) + vec(i)%deriv(j)
         enddo
      enddo
      
      if(allocated(vec)) deallocate(vec)

    end function sum_vec

    !product_vec: (ADF95_dpr)
    pure function product_vec(vec,msk) result(productvec)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:), intent(in) :: vec
      type(ADF95_dpr)                           :: productvec
      logical     , dimension(1:size(vec)), intent(in), optional :: msk
      integer(ipr)                                   :: nsize, i

      nsize = size(vec)

      productvec = 1.0_dpr
      if(present(msk)) then
         do i=1, nsize
            if(msk(i)) then
               productvec = productvec * vec(i)
            endif
         enddo
      else
         do i=1, nsize
            productvec = productvec * vec(i)
         enddo
      endif
    end function product_vec

    !minimum_value: (ADF95_dpr)
    pure function minimum_value(vec, msk) result(minvec)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:), intent(in)                     :: vec
      type(ADF95_dpr)                                               :: minvec
      logical        , dimension(1:size(vec)), intent(in), optional :: msk
      integer(ipr)   , dimension(1:1)                               :: loc
      integer(ipr) :: i, nsize, lenc

      nsize    = size(vec)
      if(present(msk)) then
         loc(1:1) = minloc(vec(1:nsize)%value,mask=msk(1:nsize))
      else
         loc(1:1) = minloc(vec(1:nsize)%value)
      endif
      
      minvec%value = vec(loc(1))%value

      minvec%deriv(1:vec(loc(1))%index(0)) = vec(loc(1))%deriv(1:vec(loc(1))%index(0))
      minvec%index(0:vec(loc(1))%index(0)) = vec(loc(1))%index(0:vec(loc(1))%index(0))

      !search for undifined situations
      do i=1, nsize
         if((i .ne. loc(1)) .and. (minvec%value .eq. vec(i)%value)) then
            if(vec(i)%index(0) .gt. vec(loc(1))%index(0)) then
               lenc = vec(i)%index(0)
               minvec%deriv(1:lenc) = vec(i)%deriv(1:lenc)
               minvec%index(0:lenc) = vec(i)%index(0:lenc)
               !Note: derivative not defined: return NaN
               minvec%deriv(vec(loc(1))%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
               lenc = vec(loc(1))%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  minvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            else if(vec(i)%index(0) .lt. vec(loc(1))%index(0)) then
               !Note: derivative not defined: return NaN
               minvec%deriv(vec(i)%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
               lenc = vec(i)%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  minvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            else
               lenc = vec(i)%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  minvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            endif
         endif
      enddo
    end function minimum_value

    !minimum_loc1: (ADF95_dpr)
    pure function minimum_loc1(vec, msk) result(minl)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:)          , intent(in)           :: vec
      logical        , dimension(1:size(vec)), intent(in), optional :: msk
      integer(ipr)   , dimension(1:1)                               :: minl
      integer :: nsize

      nsize = size(vec)
      if(present(msk)) then
         minl(1:1) = minloc(vec(1:nsize)%value,mask=msk(1:nsize))
      else
         minl(1:1) = minloc(vec(1:nsize)%value)
      endif
    end function minimum_loc1

    !maximum_value: (ADF95_dpr)
    pure function maximum_value(vec, msk) result(maxvec)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:), intent(in)                     :: vec
      type(ADF95_dpr)                                               :: maxvec
      logical        , dimension(1:size(vec)), intent(in), optional :: msk
      integer(ipr)   , dimension(1:1)                               :: loc
      integer(ipr) :: i, nsize, lenc

      nsize    = size(vec)
      if(present(msk)) then
         loc(1:1) = maxloc(vec(1:nsize)%value,mask=msk(1:nsize))
      else
         loc(1:1) = maxloc(vec(1:nsize)%value)
      endif
      
      maxvec%value = vec(loc(1))%value

      maxvec%deriv(1:vec(loc(1))%index(0)) = vec(loc(1))%deriv(1:vec(loc(1))%index(0))
      maxvec%index(0:vec(loc(1))%index(0)) = vec(loc(1))%index(0:vec(loc(1))%index(0))

      !search for undifined situations
      do i=1, nsize
         if((i .ne. loc(1)) .and. (maxvec%value .eq. vec(i)%value)) then
            if(vec(i)%index(0) .gt. vec(loc(1))%index(0)) then
               lenc = vec(i)%index(0)
               maxvec%deriv(1:lenc) = vec(i)%deriv(1:lenc)
               maxvec%index(0:lenc) = vec(i)%index(0:lenc)
               !Note: derivative not defined: return NaN
               maxvec%deriv(vec(loc(1))%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
               lenc = vec(loc(1))%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  maxvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            else if(vec(i)%index(0) .lt. vec(loc(1))%index(0)) then
               !Note: derivative not defined: return NaN
               maxvec%deriv(vec(i)%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
               lenc = vec(i)%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  maxvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            else
               lenc = vec(i)%index(0)
               where(vec(i)%deriv(1:lenc) .ne. vec(loc(1))%deriv(1:lenc))
                  !Note: derivative not defined: return NaN
                  maxvec%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
               endwhere
            endif
         endif
      enddo
    end function maximum_value

    !maximal_loc1: (ADF95_dpr)
    pure function maximal_loc1(vec, msk) result(maxl)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:)          , intent(in)           :: vec
      logical        , dimension(1:size(vec)), intent(in), optional :: msk
      integer(ipr)   , dimension(1:1)                               :: maxl

      if(present(msk)) then
         maxl(1:1) = maxloc(vec(1:size(vec))%value,mask=msk(1:size(vec)))
      else
         maxl(1:1) = maxloc(vec(1:size(vec))%value)
      endif
    end function maximal_loc1

    !operators: elemental functions that do not convert

    !dim1: (ADF95_dpr, ADF95_dpr)
    elemental function dim1(a, b) result(dima)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: dima

      dima = max(a-b,0.0_dpr)
    end function dim1

    !min_value1: (ADF95_dpr, ADF95_dpr)
    elemental function min_value1(a, b) result(minc)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)     :: a
      type(ADF95_dpr), intent(in)     :: b
      type(ADF95_dpr)                 :: minc
      integer(ipr)                    :: lenc

      if(a%value .lt. b%value) then
         lenc       = a%index(0)
         minc%value = a%value
         minc%deriv(1:lenc) = a%deriv(1:lenc)
         minc%index(0:lenc) = a%index(0:lenc)
      else if(a%value .gt. b%value) then
         lenc       = b%index(0)
         minc%value = b%value
         minc%deriv(1:lenc) = b%deriv(1:lenc)
         minc%index(0:lenc) = b%index(0:lenc)
      else
         minc%value = a%value
         if(a%index(0) .gt. b%index(0)) then
            lenc = a%index(0)
            minc%deriv(1:lenc) = a%deriv(1:lenc)
            minc%index(0:lenc) = a%index(0:lenc)
            !Note: derivative not defined: return NaN
            minc%deriv(b%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
            lenc = b%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               minc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         else if(a%index(0) .lt. b%index(0)) then
            lenc = b%index(0)
            minc%deriv(1:lenc) = b%deriv(1:lenc)
            minc%index(0:lenc) = b%index(0:lenc)
            !Note: derivative not defined: return NaN
            minc%deriv(a%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
            lenc = a%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               minc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         else
            lenc = a%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               minc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         endif
      endif
    end function min_value1

    !min_value2: (ADF95_dpr, dpr)
    elemental function min_value2(a, b) result(minc)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)     :: a
      real(dpr)      , intent(in)     :: b
      type(ADF95_dpr)                 :: minc
      integer(ipr)                    :: lenc

      if(a%value .lt. b) then
         lenc       = a%index(0)
         minc%value = a%value
         minc%deriv(1:lenc) = a%deriv(1:lenc)
         minc%index(0:lenc) = a%index(0:lenc)
      else if(a%value .gt. b) then
         minc%value    = b
         minc%index(0) = 0
      else
         minc%value = a%value
         lenc = a%index(0)
         minc%value = a%value
         minc%deriv(1:lenc) = a%deriv(1:lenc)
         minc%index(0:lenc) = a%index(0:lenc)
         where(a%deriv(1:lenc) .ne. 0.0_dpr)
            !Note: derivative not defined: return NaN
            minc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
         endwhere
      endif
    end function min_value2

    !min_value3: (dpr, ADF95_dpr)
    elemental function min_value3(a, b) result(minc)
      use mod_precision
      implicit none
      real(dpr)      , intent(in)     :: a
      type(ADF95_dpr), intent(in)     :: b
      type(ADF95_dpr)                 :: minc
      integer(ipr)                    :: lenc

      if(a .lt. b%value) then
         minc%value    = a
         minc%index(0) = 0
      else if(a .gt. b%value) then
         lenc       = b%index(0)
         minc%value = b%value
         minc%deriv(1:lenc) = b%deriv(1:lenc)
         minc%index(0:lenc) = b%index(0:lenc)
      else
         minc%value = b%value
         lenc = b%index(0)
         minc%value = b%value
         minc%deriv(1:lenc) = b%deriv(1:lenc)
         minc%index(0:lenc) = b%index(0:lenc)
         where(b%deriv(1:lenc) .ne. 0.0_dpr)
            !Note: derivative not defined: return NaN
            minc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
         endwhere
      endif
    end function min_value3

    !max_value1: (ADF95_dpr, ADF95_dpr)
    elemental function max_value1(a, b) result(maxc)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)     :: a
      type(ADF95_dpr), intent(in)     :: b
      type(ADF95_dpr)                 :: maxc
      integer(ipr)                    :: lenc

      if(a%value .gt. b%value) then
         lenc       = a%index(0)
         maxc%value = a%value
         maxc%deriv(1:lenc) = a%deriv(1:lenc)
         maxc%index(0:lenc) = a%index(0:lenc)
      else if(a%value .lt. b%value) then
         lenc       = b%index(0)
         maxc%value = b%value
         maxc%deriv(1:lenc) = b%deriv(1:lenc)
         maxc%index(0:lenc) = b%index(0:lenc)
      else
         maxc%value = a%value
         if(a%index(0) .gt. b%index(0)) then
            lenc = a%index(0)
            maxc%deriv(1:lenc) = a%deriv(1:lenc)
            maxc%index(0:lenc) = a%index(0:lenc)
            !Note: derivative not defined: return NaN
            maxc%deriv(b%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
            lenc = b%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               maxc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         else if(a%index(0) .lt. b%index(0)) then
            lenc = b%index(0)
            maxc%deriv(1:lenc) = b%deriv(1:lenc)
            maxc%index(0:lenc) = b%index(0:lenc)
            !Note: derivative not defined: return NaN
            maxc%deriv(a%index(0)+1:lenc) = -sqrt(asin(-1.0_dpr))
            lenc = a%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               maxc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         else
            lenc = a%index(0)
            where(a%deriv(1:lenc) .ne. b%deriv(1:lenc))
               !Note: derivative not defined: return NaN
               maxc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
            endwhere
         endif
      endif
    end function max_value1

    !max_value2: (ADF95_dpr, dpr)
    elemental function max_value2(a, b) result(maxc)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)     :: a
      real(dpr)      , intent(in)     :: b
      type(ADF95_dpr)                 :: maxc
      integer(ipr)                    :: lenc

      if(a%value .gt. b) then
         lenc       = a%index(0)
         maxc%value = a%value
         maxc%deriv(1:lenc) = a%deriv(1:lenc)
         maxc%index(0:lenc) = a%index(0:lenc)
      else if(a%value .lt. b) then
         maxc%value    = b
         maxc%index(0) = 0
      else
         maxc%value = a%value
         lenc = a%index(0)
         maxc%value = a%value
         maxc%deriv(1:lenc) = a%deriv(1:lenc)
         maxc%index(0:lenc) = a%index(0:lenc)
         where(a%deriv(1:lenc) .ne. 0.0_dpr)
            !Note: derivative not defined: return NaN
            maxc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
         endwhere
      endif
    end function max_value2

    !max_value3: (dpr, ADF95_dpr)
    elemental function max_value3(a, b) result(maxc)
      use mod_precision
      implicit none
      real(dpr)      , intent(in)     :: a
      type(ADF95_dpr), intent(in)     :: b
      type(ADF95_dpr)                 :: maxc
      integer(ipr)                    :: lenc

      if(a .gt. b%value) then
         maxc%value    = a
         maxc%index(0) = 0
      else if(a .lt. b%value) then
         lenc       = b%index(0)
         maxc%value = b%value
         maxc%deriv(1:lenc) = b%deriv(1:lenc)
         maxc%index(0:lenc) = b%index(0:lenc)
      else
         maxc%value = b%value
         lenc = b%index(0)
         maxc%value = b%value
         maxc%deriv(1:lenc) = b%deriv(1:lenc)
         maxc%index(0:lenc) = b%index(0:lenc)
         where(b%deriv(1:lenc) .ne. 0.0_dpr)
            !Note: derivative not defined: return NaN
            maxc%deriv(1:lenc) = -sqrt(asin(-1.0_dpr))
         endwhere
      endif
    end function max_value3

    !mod1: (ADF95_dpr, ADF95_dpr)
    elemental function mod1(a, p) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr), intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p%value .ne. 0.0_dpr) then
         div = int(a%value/p%value)
      else
         div = 0
      endif
      
      f       = a - div * p
      f%value = mod(a%value,p%value)
    end function mod1

    !mod2: (real(dpr), ADF95_dpr)
    elemental function mod2(a, p) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p%value .ne. 0.0_dpr) then
         div = int(a/p%value)
      else
         div = 0
      endif
      
      f       = a - div * p
      f%value = mod(a,p%value)
    end function mod2

    !mod3: (ADF95_dpr, real(dpr))
    elemental function mod3(a, p) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p .ne. 0.0_dpr) then
         div = int(a%value/p)
      else
         div = 0
      endif
      
      f       = a - div * p
      f%value = mod(a%value,p)
    end function mod3

    !modulo1: (ADF95_dpr, ADF95_dpr)
    elemental function modulo1(a, p) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr), intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p%value .ne. 0.0_dpr) then
         div = floor(a%value/p%value)
      else
         div = 0
      endif
      
      f = a - div * p
      f%value = modulo(a%value,p%value)
    end function modulo1

    !modulo2: (real(dpr), ADF95_dpr)
    elemental function modulo2(a, p) result(f)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p%value .ne. 0.0_dpr) then
         div = floor(a/p%value)
      else
         div = 0
      endif
      
      f = a - div * p
      f%value = modulo(a,p%value)
    end function modulo2

    !modulo3: (ADF95_dpr, real(dpr))
    elemental function modulo3(a, p) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: p
      type(ADF95_dpr)             :: f
      integer(ipr)                :: div

      if(p .ne. 0.0_dpr) then
         div = floor(a%value/p)
      else
         div = 0
      endif
      
      f = a - div * p
      f%value = modulo(a%value,p)
    end function modulo3

    !signum1: (ADF95_dpr, dpr)
    elemental function signum1(a, b) result(signum)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: b
      type(ADF95_dpr)             :: signum
      integer(ipr)                :: lena

      lena = a%index(0)

      signum%value = sign(a%value, b)

      signum%deriv(1:lena) = sign(1.0_dpr,a%value*b) * a%deriv(1:lena)
      signum%index(0:lena) = a%index(0:lena)
    end function signum1

    !signum2: (real(dpr), ADF95_dpr)
    elemental function signum2(a, b) result(signum)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: b
      type(ADF95_dpr)             :: signum

      signum%value    = sign(a, b%value)
      signum%index(0) = 0
    end function signum2

    !signum3: (ADF95_dpr, ADF95_dpr)
    elemental function signum3(a, b) result(signum)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a, b
      type(ADF95_dpr)             :: signum
      integer(ipr)                :: lena

      lena = a%index(0)

      signum%value = sign(a%value, b%value)

      signum%deriv(1:lena) = sign(1.0_dpr,a%value*b%value) * a%deriv(1:lena)
      signum%index(0:lena) = a%index(0:lena)
    end function signum3

    !operators: elemental functions that may convert

    !abs1: (ADF95_dpr)
    elemental function absolute_value1(a) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr)             :: f
      integer(ipr)                :: lena

      lena = a%index(0)
      
      f%value = abs(a%value)

      if(a%value .eq. 0.0_dpr) then
         where(a%deriv(1:lena) .eq. 0.0_dpr)
            f%deriv(1:lena) = sign(1.0_dpr,a%value) * a%deriv(1:lena)
         elsewhere
            !Note: Undefined derivative at a%value = 0; discontinous function
            !assign NaN (Not a Number) Sign information returned!
            f%deriv(1:lena) = -sqrt(asin(-1.0_dpr))*sign(1.0_dpr,a%value)
         endwhere
      else
         f%deriv(1:lena) = sign(1.0_dpr,a%value) * a%deriv(1:lena)
      endif
      f%index(0:lena) = a%index(0:lena)
    end function absolute_value1

    !aint1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function aint1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in)           :: a
      real(dpr)                             :: b

      b = aint(a%value)
    end function aint1

    !anint1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function anint1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = anint(a%value)
    end function anint1

    !ceiling1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function ceiling1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = ceiling(a%value)
    end function ceiling1

    !floor1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function floor1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = floor(a%value)
    end function floor1

    !int1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function int1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = int(a%value)
    end function int1

    !nint1: (ADF95_dpr); optional kind cannot be inplemented
    elemental function nint1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = nint(a%value)
    end function nint1

    !operator: numerical inquiry functions

    !digits1: (ADF95_dpr)
    elemental function digits1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = digits(a%value)
    end function digits1

    !epsilon1: (ADF95_dpr)
    elemental function epsilon1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = epsilon(a%value)
    end function epsilon1

    !huge1: (ADF95_dpr)
    elemental function huge1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = huge(a%value)
    end function huge1

    !maxexponent1: (ADF95_dpr)
    elemental function maxexponent1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = maxexponent(a%value)
    end function maxexponent1

    !minexponent1: (ADF95_dpr)
    elemental function minexponent1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = minexponent(a%value)
    end function minexponent1

    !precision1: (ADF95_dpr)
    elemental function precision1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = precision(a%value)
    end function precision1

    !radix1: (ADF95_dpr)
    elemental function radix1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = radix(a%value)
    end function radix1

    !range1: (ADF95_dpr)
    elemental function range1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = range(a%value)
    end function range1

    !tiny1: (ADF95_dpr)
    elemental function tiny1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = tiny(a%value)
    end function tiny1

    !operator: elemental functions to manipulate reals

    !exponent_iqr1: (ADF95_dpr)
    elemental function exponent_iqr1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)                :: b

      b = exponent(a%value)
    end function exponent_iqr1

    !fraction1: (ADF95_dpr)
    elemental function fraction1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = fraction(a%value)
    end function fraction1

    !nearest1: (ADF95_dpr, ADF95_dpr)
    elemental function nearest1(a, s) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      type(ADF95_dpr), intent(in) :: s
      real(dpr)                   :: b

      b = nearest(a%value,s%value)
    end function nearest1

    !nearest2: (ADF95_dpr, dpr)
    elemental function nearest2(a, s) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)      , intent(in) :: s
      real(dpr)                   :: b

      b = nearest(a%value,s)
    end function nearest2

    !nearest3: (dpr, ADF95_dpr)
    elemental function nearest3(a, s) result(b)
      use mod_precision
      implicit none
      real(dpr)      , intent(in) :: a
      type(ADF95_dpr), intent(in) :: s
      real(dpr)                   :: b

      b = nearest(a,s%value)
    end function nearest3

    !rrspacing1: (ADF95_dpr)
    elemental function rrspacing1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = rrspacing(a%value)
    end function rrspacing1

    !scale1: (ADF95_dpr)
    elemental function scale1(a, i) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: i
      type(ADF95_dpr)             :: b
      integer                     :: lena

      lena = a%index(0)

      b%value = scale(a%value,i)

      b%deriv(1:lena) = scale(a%deriv(1:lena),i)
      b%index(0:lena) = a%index(0:lena)
    end function 

    !set_exponent1: (ADF95_dpr)
    elemental function set_exponent1(a, i) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      integer(ipr)   , intent(in) :: i
      type(ADF95_dpr)             :: b
      integer                     :: lena

      lena = a%index(0)

      b%value = set_exponent(a%value,i)

      b%deriv(1:lena) = set_exponent(a%deriv(1:lena),i)
      b%index(0:lena) = a%index(0:lena)
    end function set_exponent1

    !spacing1: (ADF95_dpr)
    elemental function spacing1(a) result(b)
      use mod_precision
      implicit none
      type(ADF95_dpr), intent(in) :: a
      real(dpr)                   :: b

      b = spacing(a%value)
    end function spacing1

    !operator: vector and matrix multipliacation functions

    !does not work with Absoft (Version 8.2a): BEGIN comment out
    !dot_product1: (ADF95_dpr, ADF95_dpr)
    pure function dot_product1(a,b) result(dotp)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:)        , intent(in)  :: a
      type(ADF95_dpr), dimension(:)        , intent(in)  :: b
      type(ADF95_dpr)                                    :: dotp

      dotp = sum(a*b)

    end function dot_product1

    !matmul1: (ADF95_dpr, ADF95_dpr)
    pure function matmul1(a,b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:,:), intent(in)         :: a
      type(ADF95_dpr), dimension(:,:), intent(in)         :: b
      type(ADF95_dpr), dimension(1:size(a,1),1:size(b,2)) :: f
      integer :: i, j

      forall(i=1:size(a,1))
         forall(j=1:size(b,2))
            f(i,j) = sum(a(i,:)*b(:,j))
         end forall
      end forall
    end function matmul1

    !matmul2: (ADF95_dpr, ADF95_dpr)
    pure function matmul2(a,b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:)  , intent(in) :: a
      type(ADF95_dpr), dimension(:,:), intent(in) :: b
      type(ADF95_dpr), dimension(1:size(b,2))     :: f
      integer :: j

      forall(j=1:size(b,2))
         f(j) = sum(a(:)*b(:,j))
      end forall
    end function matmul2

    !matmul3: (ADF95_dpr, ADF95_dpr)
    pure function matmul3(a,b) result(f)
      use mod_precision
      implicit none
      type(ADF95_dpr), dimension(:,:), intent(in) :: a
      type(ADF95_dpr), dimension(:)  , intent(in) :: b
      type(ADF95_dpr), dimension(1:size(a,1))     :: f
      integer :: i

      forall(i=1:size(a,1))
         f(i) = sum(a(i,:)*b(:))
      end forall
    end function matmul3
    !matmul does not work with Absoft (Version 8.2a): END comment out

end module mod_ADF95types

!interfaces for user functions
module mod_ADF95interfaces

  use mod_ADF95types
  
  implicit none

  !independent
  interface ADF95_independent
     module procedure ADF95_independent1, ADF95_independent2, &
          & ADF95_independent3
  end interface

  !value
  interface ADF95_value
     module procedure ADF95_value_ADF, ADF95_value_DPR
  end interface

  !first derivative
  interface ADF95_deriv
     module procedure ADF95_deriv1
  end interface

  !non-zero element fill-in
  interface ADF95_fillin
     module procedure ADF95_fillin1
  end interface

  !jacobian
  interface ADF95_jacobian
     module procedure ADF95_jacobian1
  end interface

  contains

    !independent
    elemental subroutine ADF95_independent1(icc,x,val)
      use mod_precision
      use mod_ADF95types
      type(ADF95_dpr), intent(inout) :: x
      integer(ipr)   , intent(in)    :: icc
      real(dpr)      , intent(in)    :: val

      x%value    = val
      x%deriv(1) = 1.0_dpr
      x%index(1) = icc
      x%index(0) = 1
    end subroutine ADF95_independent1

    elemental subroutine ADF95_independent2(icc,x,val)
      use mod_precision
      use mod_ADF95types
      type(ADF95_dpr), intent(inout) :: x
      integer(ipr)   , intent(in)    :: icc
      real(spr)      , intent(in)    :: val

      x%value    = real(val,dpr)
      x%deriv(1) = 1.0_dpr
      x%index(1) = icc
      x%index(0) = 1
    end subroutine ADF95_independent2

    elemental subroutine ADF95_independent3(icc,x,val)
      use mod_precision
      use mod_ADF95types
      type(ADF95_dpr), intent(inout) :: x
      integer(ipr)   , intent(in)    :: icc
      integer(ipr)   , intent(in)    :: val

      x%value    = real(val,dpr)
      x%deriv(1) = 1.0_dpr
      x%index(1) = icc
      x%index(0) = 1
    end subroutine ADF95_independent3

    !value
    elemental function ADF95_value_ADF(x) result(f)
      use mod_precision
      use mod_ADF95types
      implicit none
      type(ADF95_dpr), intent(in) :: x
      real(dpr)                   :: f

      f = x%value
    end function ADF95_value_ADF

    elemental function ADF95_value_DPR(x) result(f)
      use mod_precision
      use mod_ADF95types
      implicit none
      real(dpr), intent(in) :: x
      real(dpr)             :: f

      f = x
    end function ADF95_value_DPR

    !first derivative
    elemental function ADF95_deriv1(x, icc) result(df)
      use mod_precision
      use mod_ADF95types, only : ADF95_dpr
      implicit none
      type(ADF95_dpr), intent(in) :: x
      integer(ipr)   , intent(in) :: icc
      real(dpr)                   :: df
  
      integer(ipr) :: i

      df = 0.0_dpr
      do i=1, x%index(0)
         if(x%index(i) .eq. icc) then
            df = x%deriv(i)
            exit
         endif
      enddo
    end function ADF95_deriv1

    !find minimal number of upper and lower bands with non-zero entries
    pure subroutine ADF95_fillin1(xvec, LDsize_opt, ml, mu)
      use mod_precision
      use mod_ADF95types, only : ADF95_dpr
      implicit none
      type(ADF95_dpr), dimension(:), intent(in)  :: xvec
      integer(ipr)                 , intent(out) :: LDsize_opt
      integer(ipr)   , optional    , intent(out) :: ml, mu
      integer(ipr) :: i, j, nsize, num
  
      nsize = size(xvec)

      if(present(ml)) then
         ml = 0
         do i=1, nsize
            num = 0
            do j=1, xvec(i)%index(0)
               if(xvec(i)%index(j) .lt. i) num = num + 1
            enddo
            ml = max(ml,num)
         enddo
      endif
  
      if(present(mu)) then
         mu = 0
         do i=1, nsize
            num = 0
            do j=1, xvec(i)%index(0)
               if(xvec(i)%index(j) .gt. i) num = num + 1
            enddo
            mu = max(mu,num)
         enddo
      endif

      LDsize_opt = maxval(xvec(1:nsize)%index(0))
    end subroutine ADF95_fillin1

    !ADF95_jacobian
    pure subroutine ADF95_jacobian1(f, df, ir, ic, nz)
      use mod_precision
      use mod_ADF95types
      implicit none
      type(ADF95_dpr), dimension(:), intent(in)  :: f
      real(dpr)      , dimension(:), intent(out) :: df
      integer(ipr)   , dimension(:), intent(out) :: ir, ic
      integer(ipr)                 , intent(out) :: nz
      integer(ipr) :: i, j, k, sf, sdf

      sf  = size(f)
      sdf = min(size(df),size(ir),size(ic))

      nz = -1
      k  =  0
      do i=1, sf
         do j=1, f(i)%index(0)
            k = k + 1
            if(k .gt. sdf) return
            df(k) = f(i)%deriv(j)
            ir(k) = i
            ic(k) = f(i)%index(j)
         enddo
      enddo

      nz = k

    end subroutine ADF95_jacobian1

end module mod_ADF95interfaces

!gather modules
module mod_ADF95

  use mod_ADF95types
  use mod_ADF95interfaces

end module mod_ADF95
