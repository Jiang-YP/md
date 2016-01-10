      INTEGER FUNCTION MDVAR(NVARS,ICOL,IROW)
      PARAMETER (IPARAM_NVARS=32)
      INTEGER IVARIABLES(8,IPARAM_NVARS), IVRSN
      DATA IVARIABLES/
     1  4HOPT1,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT2,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT3,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT4,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT5,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT6,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT7,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT8,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT9,  4H    ,  4H    ,  4H    ,  0,  -1,  0,  1,
     1  4HOPT1,  4H0   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H1   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H2   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H3   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H4   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H5   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HOPT1,  4H6   ,  4H    ,  4H    ,  0,  -1,  0,  2,
     1  4HHFNU,  4HMBER,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHFLE,  4HNGTH,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHFIN,  4HNERD,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHFOU,  4HTERD,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHOUS,  4HEID ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HHOUS,  4HEOD ,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HTHER,  4HMCND,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HPORO,  4HSITY,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HPORE,  4HSIZE,  4H    ,  4H    ,  1,  -1,  0,  1,
     1  4HRE1 ,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HRE2 ,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HNU1 ,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HNU2 ,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HJM  ,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HVAR1,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2,
     1  4HVAR2,  4H    ,  4H    ,  4H    ,  1,  -1,  0,  2/

      DATA IVRSN/1451940121/
      NVARS = IPARAM_NVARS
      IF(IROW .EQ. 0) THEN
          MDVAR=IVRSN
      ELSE IF(IROW .LE. NVARS) THEN
          MDVAR=IVARIABLES(ICOL,IROW)
      ENDIF
      RETURN
      END
