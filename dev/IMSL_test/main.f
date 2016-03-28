      INCLUDE 'link_fnl_static.h'
      USE IMSL_LIBRARIES
      IMPLICIT NONE
      INTEGER LDSUR, NDATA, NXOUT, NYOUT
      PARAMETER (NDATA=20, NXOUT=3, NYOUT=3, LDSUR=NXOUT)
!
      INTEGER I, J, NOUT
      REAL ABS, COS, F, FDATA(NDATA), FLOAT, PI,
     &     SIN, SUR(LDSUR,NYOUT), X, XOUT(NXOUT),
     &     XYDATA(2,NDATA), Y, YOUT(NYOUT)
      INTRINSIC ABS, COS, FLOAT, SIN
!     Define function
      F(X,Y) = 3.0 + 7.0*X + 2.0*Y
!     Get value for PI
      PI = CONST('PI')
!     Set up X, Y, and F data on a circle
      DO 10 I=1, NDATA
      XYDATA(1,I) = 3.0*SIN(2.0*PI*FLOAT(I-1)/FLOAT(NDATA))
      XYDATA(2,I) = 3.0*COS(2.0*PI*FLOAT(I-1)/FLOAT(NDATA))
      FDATA(I) = F(XYDATA(1,I),XYDATA(2,I))
   10 CONTINUE
!     Set up XOUT and YOUT data on [0,1] by
!     [0,1] grid.
      DO 20 I=1, NXOUT
      XOUT(I) = FLOAT(I-1)/FLOAT(NXOUT-1)
   20 CONTINUE
      DO 30 I=1, NXOUT
      YOUT(I) = FLOAT(I-1)/FLOAT(NYOUT-1)
   30 CONTINUE
!     Interpolate scattered data
      CALL SURF (XYDATA, FDATA, XOUT, YOUT, SUR)
!     Get output unit number
      CALL UMACH (2, NOUT)
!     Write heading
      WRITE (NOUT,99998)
!     Print results
      DO 40 I=1, NYOUT
      DO 40 J=1, NXOUT
      WRITE (NOUT,99999) XOUT(J), YOUT(I), SUR(J,I),
     &        F(XOUT(J),YOUT(I)),
     &        ABS(SUR(J,I)-F(XOUT(J),YOUT(I)))
   40 CONTINUE
99998 FORMAT (' ', 10X, 'X', 11X, 'Y', 9X, 'SURF', 6X, 'F(X,Y)', 7X,
     &         'ERROR', /)
99999 FORMAT (1X, 5F12.4)
      END