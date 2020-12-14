*DECK LHS
      SUBROUTINE LHS
       !DEC$ ATTRIBUTES DLLEXPORT::LHS
c this program is to generate the random number for the desired parameter. There are different
c method for that. e.g. Uniform distribution, normal distribution,etc.
c  we don't need change this main program, but we need change the input file LHSinput.dat.
c parameters for changing includes: 1--NOBS (number of random ), usually set to 200 250
c 2--distribution function, for example, uniform or normal, lognormal, etc, and then use 
c two parameters to control the shape as shown in the LHSinput.dat
c 3--Define the correlation coefficient between different parameters.
c 4--the last step is to get the random array from the output file "LHSout.dat", and 
c then name the file as "paralist.txt" for the next program.
C***********************************************************************
C           
C MODULE:     LHS2_LHS.FOR
C VERSION:    2.51
C DATE:       10/20/03
C TITLE:      Latin Hypercube or random sampling program
C
C  ABSTRACT
C  --------
C       The purpose of the LHS program is to sample distributions of input
C  parameters using either normal Monte Carlo sampling or efficient Latin 
C  Hypercube Sampling.  LHS permits correlations (restricted pairings) between
C  parameters.  Latin Hypercube sampling reduces the minimum number of sample 
C  sets [nV] required to about 4/3 * nV where nV is number of uncertain 
C  (varying) parameters (Iman and Shortencarier, 1984).  
CLHSinput
C  AUTHORS
C  -------
C  R. L. Iman                         M.J. Shortencarier
C  Department 6415                    Department 7223
C  Sandia National Laboratories       Sandia National Laboratories
C  Albuquerque, NM  87185             Albuquerque, NM  87185
C  (505) 844-8334                     (505) 846-1662
C
C  SPONSOR                            CONSULTANT
C  -------                            ----------
C  Lanny Smith                        Harrold Iuzzolino
C  (505) 766-9629                     (505) 766-9629
C
C  SOURCE
C  ------
C  Sandia National Laboratories
C  WIPP Performance Assessment
C  Department 6749
C  Albuquerque, NM  87185-5800
C
C  LANGUAGE
C  --------
C  FORTRAN 77
C
C  HARDWARE / SOFTWARE PLATFORM
C  ------------------- --------
C  DEC Alpha OpenVMS AXP version 6.1
C
C  INPUT / OUTPUT FILES
C  -------------- -----
C  LHS2_DBG$OUTPUT    1   Output debug output file
C  LHS2_UIF$INPUT     5   Input ASCII user input file
C  LHS2_NO2$SCRATCH   2   Temporary scratch file
C  LHS2_NO3$SCRATCH   3   Temporary scratch file
C  LHS2_NO4$SCRATCH   4   Temporary scratch file
C  LHS2_OUT$OUTPUT    6   ASCII text output file
C  LHS2_NO7$SCRATCH   7   Temporary scratch file
C  LHS2_NO8$SCRATCH   8   Temporary scratch file
C  LHS2_NO9$SCRATCH   9   Temporary scratch file
C  LHS2_PLOT$OUTPUT  20   SPLAT plot output file
C
C  PRIMARY REFERENCE
C  ------- ---------
C       Ronald L. Iman and Michael J. Shortencarier, 1984.  A FORTRAN 77 
C  Program and User's Guide for the Generation of Latin Hypercube and Random 
C  Samples for Use With Computer Models.  NUREG/CR-3624, SAND83-2365.
C  Albuquerque, NM: Sandia National Laboratories.
C
C  DISCLAIMER
C  ----------
C      This computer program was prepared as an account of work
C  sponsored by an agency of the United States Government. Neither
C  the United States Government nor any agency thereof, nor any of
C  their employees, nor any of their contractors, subcontractors,
C  or their employees, makes any warranty, express or implied, or
C  assumes any legal liability or responsibility for the accuracy,
C  completeness, or usefulness of any information, apparatus,
C  product, or process disclosed, or represents that its use would
C  not infringe privately owned rights. Reference herein to any
C  specific commercial product, process, or service by trade name,
C  trademark, manufacturer, or otherwise, does not necessarily
C  constitute or imply its endorsement, recommendation, or
C  favoring by the United States Government, any agency thereof
C  or any of their contractors or subcontractors. The views and
C  opinions expressed herein do not necessarily state or reflect
C  those of the United States Government, any agency thereof or
C  any of their contractors or subcontractors.
C
C  UPDATE HISTORY
C  ------ -------
C  04/09/97  LNS  Added capablilty to create an output file by specification
C                 of a keyword, PLOT, that provides sample values and 
C                 cumulative probabilities for each distribution evaluated.
C                 This output file is in a format suitable for use by the SPLAT
C                 plotting program.
C            HJI  Fixed beta distribution bug by replacing routines
C                 BETALN and GAMALN with new versions. 
C                 Replaced subroutine FINVER (Beta function inverse) by 
C                 new subroutine INVBETA to allow very small P and Q to
C                 be used in Beta functions.
C  03/06/96  LNS  Modified coding associated with Student and LogStudent
C                 distributions to force all sampled points to fall between
C                 the upper and lower bounds specified for the distribution. 
C  01/31/96  EJD  Changed the comment line containing the version number of
C                 the code, because the code was linked with new versions
C                 of the CAMDAT, CAMCON, and CAMSUPES libraries.  No executable
C                 lines in the code were modified.
C  09/11/95  LNS  Modified code to bring it into compliance with PRETEST
C                 criteria prior to installation into the WIPP PA CMS.
C  08/13/93  LNS  Repaired bug induced code execution errors associated with 
C                 the introduction of the Student distribution into the code.
C  05/21/93  HJI  Added log Student distribution.  The keyword is 'LOGSTUDENT'
C  05/18/93  HJI  Added Student distribution.  The keyword is 'STUDENT '.
C                 Replaced call to RDCMLN by call to FILCMDLIN.
C  02/25/93  HJI  Modified RAYLEIGH distribution so that the constant in
C                 the exponent is the maximum intrusion rate (lambda)
C                 instead of the rate of change of the intrusion rate.
C  05/20/92  HJI  Changed SMPLHS output format to 10 columns of width
C                 11 columns to prevent negative numbers from being
C                 adjacent to the previous column.
C  08/22/91  HJI  Changed maximum number of vectors (NMAX) from 1000 to 10000.
C  03/24/91  HJI  Changed minimum correlation matrix keyword from "C" to "CORR".
C  08/15/90  HJI  Changed endpoints of the NORMAL and LOGNORMAL distributions 
C                 to be the 0.01 and 0.99 percentile points.
C  07/30/90  HJI  Implemented both increasing and decreasing exponential
C                 distributions: f(X)= EXP(C*X)*C/(EXP(C*B)-EXP(C*A)),
C                 where C can be positive or negative.
C  01/03/90  HJI  Typed TITLE as a CHARACTER variable, moved it to
C                 COMMON block HEADNG, and shortened it to 100
C                 characters (for FORTRAN 77 standard compatibility).
C  12/15/89  HJI  Added exponential, rayleigh, and rayleigh-exponential 
C                 distributions.  The keywords are 'EXPONENTIAL', 'RAYLEIGH',
C                 and 'RAYLEXP'.
C  07/20/88  HJI  Added READONLY to OPEN of input file.
C  05/09/88  HJI  Added OPEN statements for files 1,2,3,4,6,7,8, and 9.
C  04/21/88  HJI  Added call to RDCMLN to allow the input file name to be 
C                 placed on the command line.
C <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      COMMON/WORKX/XX
      REAL*4       XX(NMAX*NVAR)
      COMMON/WORKC/Q, S
      REAL*4       Q((NVAR*(NVAR+1))/2),S((NVAR*(NVAR+1))/2)

      INTEGER*4 LOC, I, J, LOC1, ITER, NVX, LB, K, IREP, IDT
      INTEGER*4 ISTUD
c************************************************
      character*40 str, strr,chr
c*************************************************      

      LOC(I,J)=(J-1)*N+I
      LOC1(I,J)=J+(I*I-I)/2

C
C   Subroutine FILCMDLIN reads the input file name from the command line.
C   If there is no file name on the command line, input will be from
C   SYS$INPUT.
C
C      CALL FILCMDLIN(3,NFILES,FILESP)
C...Now open the file
C      IF (FILESP(1) .NE. ' ') THEN
C         OPEN (UNIT=5,FILE=FILESP(1),STATUS='OLD',READONLY)
C      ELSE
c**********change for a new input file************************
c	 OPEN (UNIT=5,FILE='LHS2_UIF$INPUT',STATUS='OLD',READONLY)
c      write(*,*) 'Which input file do you want to use?'
c      read(*,'(a40)') str
      open (unit=5, file='LHSinput.dat', status='old')
c************************************************************	 
C      ENDIF
C      IF (FILESP(3) .NE. ' ') THEN
C         OPEN (UNIT=1,FILE=FILESP(3),STATUS='NEW')
C      ELSE
c**********change for a new DBG file*********************
c	 OPEN (UNIT=1,FILE='LHS2_DBG$OUTPUT',STATUS='NEW')
c	 write(*,*) 'Which DBG file do you want to use?'
c	 read(*,'(a40)') strr
	 open (unit=1, file='LHSdbg.out', status='unknown')
c***********************************************************	 
C      ENDIF
      OPEN (UNIT=2,FORM='UNFORMATTED',
     &   STATUS='SCRATCH')
      OPEN (UNIT=3,FORM='UNFORMATTED',
     &   STATUS='SCRATCH')
      OPEN (UNIT=4,FORM='UNFORMATTED',
     &   STATUS='SCRATCH')
C      IF (FILESP(2) .NE. ' ') THEN
C         OPEN (UNIT=6,FILE=FILESP(2),STATUS='NEW')
C      ELSE
c***********change for a new output file******************
c	 OPEN (UNIT=6,FILE='LHS2_OUT$OUTPUT',STATUS='NEW')
c	 write(*,*) 'Which output file do you want to use?'
c	 read(*,'(a40)') chr
	 open (unit=6,file='LHSout.out', status='unknown')
c**********************************************************	 
C      ENDIF
      OPEN (UNIT=7,STATUS='SCRATCH')
      OPEN (UNIT=8,FORM='UNFORMATTED',
     &   STATUS='SCRATCH')
      OPEN (UNIT=9,FORM='UNFORMATTED',
     &   STATUS='SCRATCH')
C
C SUBROUTINE RDPAR READS IN THE PARAMETER STATEMENTS AND DEFINES
C THE VARIABLES IN COMMON /PARAM/
C SUBROUTINES BANNER AND WRTPAR ECHO A DESCRIPTION OF THE SAMPLE
C AS INPUT BY THE USER AND A DESCRIPTION OF THE DISTRIBUTIONS
C USED WITH THE INPUT VARIABLES
C
      ISTUD=1
      CALL RDPAR
      CALL BANNER(1)
      CALL WRTPAR
C
C IF THE USER HAS SPECIFIED A CORRELATION STRUCTURE (ICM=1) THEN
C THE CORRELATION MATRIX IS ECHOED AND CHECKED TO MAKE SURE THAT
C IT IS POSITIVE DEFINITE. THE CHOLESKY FACTORIZATION IS COMPUTED
C TO BE USED LATER AS PART OF THE PROCESS FOR INDUCING THE DESIRED
C CORRELTION STRUCTURE
C
      IF(ICM.NE.0)THEN
	CALL PMTRX(NCM,3)
	IF(IRP.EQ.1)THEN
	  IRP=0
	  WRITE(6,9003)
	ENDIF
	IF(N.LE.NV)THEN
	  ICM=0
	  WRITE(6,9001)
	  GO TO 10
	ENDIF
	CALL POSDEF(ITER)
	IF(ITER.GT.1)THEN
	  WRITE(6,9002)
	  CALL PMTRX(NCM,4)
	ENDIF
	NVX=(NV*(NV+1))/2
	DO 130 I=1,NVX
  130   CORR(I)=0.0
	DO 140 I=1,NV
  140   CORR(LOC1(I,I))=1.0
	LB=1
	REWIND 3
	READ(3)S
	DO 170 I=1,NCM
	  DO 160 K=LB,NCM
  160     CORR(LOC1(LCM(K),LCM(I)))=S(LOC1(K,I))
	  LB=LB+1
  170   CONTINUE
	CALL CHLSKY
	REWIND 3
	WRITE(3)Q       
C
C THE FOLLOWING LINES ARE INCLUDED TO HANDLE THE
C PATHOLOGICAL CASE N=3, NV=2 AND ICM=0
C
      ELSE IF(N.EQ.3.AND.NV.EQ.2)THEN
	CORR(1)=1.0
	CORR(2)=0.0001
	CORR(3)=1.0
	CALL CHLSKY
	REWIND 3
	WRITE(3)Q
      ENDIF
   10 CONTINUE
C
C THIS LOOP IS EXECUTED ONCE FOR EACH REPETITION REQUESTED.
C UNITS 7, 8, AND 9 HAVE BEEN DEFINED IN SUBROUTINE RDPAR
C
      DO 1000 IREP=1,NREP
	REWIND 7
	REWIND 8
	REWIND 9
	IF(IREP.GT.1)CALL BANNER(IREP)
C
C GENERATE THE DISTRIBUTION REQUESTED FOR VARIABLE J
C
	DO 100 J=1,NV
	  IDT=IDIST(J)
	  IF(IDT.EQ.1)THEN
	    CALL BETA(J)
	    GO TO 100
	  ELSE IF(IDT.EQ.2.OR.IDT.EQ.3)THEN
	    CALL NORMAL(J,IDT)
	    GO TO 100
	  ELSE IF(IDT.GE.4.AND.IDT.LE.7)THEN
	    CALL UNIFRM(J,IDT)
	    GO TO 100
	  ELSE IF(IDT.EQ.8)THEN
	    CALL TRIANG(J)
	    GO TO 100
	  ELSE IF(IDT.EQ.9)THEN
	    CALL USRDST(J)
	    GO TO 100
	  ELSE IF (IDT.EQ.10)THEN
	    CALL EXPONT(J)
	    GO TO 100
	  ELSE IF (IDT.EQ.11)THEN
	    CALL RAYLEIG(J)
	    GO TO 100
	  ELSE IF (IDT.EQ.12)THEN
	    CALL RAYLEXP(J)
	    GO TO 100
	  ELSE IF (IDT .EQ. 13 .OR. IDT .EQ. 14) THEN
	    CALL STUDENT(J,IDT,ISTUD)
	    GO TO 100
	  ENDIF
  100   CONTINUE
	IF(N.EQ.1)GO TO 270
C
C IF A RANDOM SAMPLE HAS BEEN GENERATED IT MUST BE SORTED FROM
C FROM SMALLEST TO LARGEST ON EACH VARIABLE AS PART OF THE
C STRUCTURING TO REDUCE THE POSSIBILITY OF SPURIOUS CORRELATIONS
C
	IF (IRS .NE. 1) GO TO 260
	DO 250 I=1,NV
	  DO 240 J=1,N
  240     XV(J)=X(LOC(J,I))
	  CALL SIFT (XV,N)
	  DO 250 J=1,N
	  X(LOC(J,I))=XV(J)
  250   CONTINUE
  260   CONTINUE
	REWIND 4
	WRITE(4)X
C
C SUBROUTINE MIX PAIRS THE SAMPLE OBSERVATIONS TO MATCH THE
C DESIRED CORRELATION STRUCTURE
C
	CALL MIX
C
C THE SAMPLE IS WRITTEN OUT TO UNIT 1
C
  270   CONTINUE
	DO 290 I=1,N
  290   WRITE(1,911)I,NV,(X(LOC(I,J)),J=1,NV)
c*******************************************
  911 format(I5,2x,I5, 300e12.4)
c******************************************
	REWIND 4
	WRITE(4)X
C
C SUBROUTINE DATOUT OUTPUTS THE LATEST SAMPLE AND CORRESPONDING RANKS
C SUBROUTINE HSTOUT OUTPUTS HISTOGRAMS FOR THE CURRENT SAMPLE
C SUBROUTINE COROUT OUTPUTS RAW AND RANK CORRELATIONS FOR THE
C CURRENT SAMPLE
C
	IF(IDATA.NE.0)CALL DATOUT
	IF(IHIST.NE.0)CALL HSTOUT
	IF(ICORR.NE.0)CALL COROUT
 1000 CONTINUE
C
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      IF(ISPLAT.EQ.1) CLOSE(20)
C
      RETURN
 9001 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'SINCE THE SAMPLE SIZE IS LESS THAN OR EQUAL TO THE ',
     2       'NUMBER OF VARIABLES',/,' THIS IS NOT A FULL RANK CASE ',
     3       'SO THE REQUESTED',/,' CORRELATION STRUCTURE, SHOWN ',
     4       'ABOVE, CANNOT BE GENERATED',/,' THEREFORE THE INPUT ',
     5       'MATRIX WILL BE GENERATED AS IF THE',/,' INPUT ',
     6       'VARIABLES WERE INDEPENDENT',//,1X,115('*'))
 9002 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'THE INPUT RANK CORRELATION MATRIX IS NOT POSITIVE ',
     2       'DEFINITE',/,' AN ITERATIVE PROCEDURE HAS BEEN USED TO ',
     3       'PRODUCE A SUBSTITUTE RANK CORRELATION MATRIX',/,' THIS ',
     3       'ADJUSTED RANK CORRELATION MATRIX APPEARS ON THE NEXT ',
     4       'PAGE',/,' THE USER SHOULD EXAMINE THIS MATRIX TO MAKE ',
     5       'SURE THAT THE CORRELATION REQUIREMENTS ARE STILL ',
     6        'SATISFIED',//,1X,115('*'))
 9003 FORMAT('0',3(9('*'),' CAUTION USER PLEASE NOTE '),10('*'),//,1X,
     1       'RANDOM PAIRING CANNOT BE USED WHEN THE USER SPECIFIES ',
     2       'A CORRELATION STRUCTURE',/,' THUS, THE RANDOM PAIRING ',
     3       'PARAMETER HAS BEEN IGNORED',//,1X,115('*'))
      END
      SUBROUTINE BANNER(IREP)
C***********************************************************************
C SUBROUTINE BANNER ECHOES USER INPUT CONCERNING THE SAMPLE
C SETTING. IT IS CALLED AT THE START OF EACH REPETITION OF THE
C SAMPLE AND IS USED TO DELIMIT REPETITIONS IN THE PRINTOUT
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      INTEGER*4 IREP
C
C IF IREP GT 1 THEN RETRIEVE THE CURRENT VALUE OF THE RANDOM
C SEED (ISEED) HERE (IF NECESSARY)
C
      WRITE(6,9001) TITLE,ISEED,NV,N
      IF(NREP.GT.1)WRITE(6,9006)IREP,NREP
      IF(IRP.EQ.1)WRITE(6,9007)
      IF(ICM.EQ.1)WRITE(6,9002)
      IF(IDATA.EQ.1)WRITE(6,9003)
      IF(IHIST.EQ.1)WRITE(6,9004)
      IF(ICORR.EQ.1)WRITE(6,9005)
      IF(ISPLAT.EQ.1)WRITE(6,9008)
      RETURN
 9001 FORMAT('1'//4X,A100//,4X,'RANDOM SEED = ',I11,//,4X,
     1       'NUMBER OF VARIABLES = ',I3,//,4X,'NUMBER ',
     2       'OF OBSERVATIONS = ',I4)
 9002 FORMAT('0',3X,'AN INPUT CORRELATION MATRIX HAS BEEN SPECIFIED')
 9003 FORMAT('0',3X,'THE SAMPLE INPUT VECTORS WILL BE PRINTED ',
     1       'ALONG WITH THEIR CORRESPONDING RANKS')
 9004 FORMAT('0',3X,'HISTOGRAMS OF THE ACTUAL SAMPLE WILL BE PLOTTED ',
     1       'FOR EACH INPUT VARIABLE')
 9005 FORMAT('0',3X,'THE CORRELATION MATRICES (RAW DATA AND RANK ',
     1       'CORRELATIONS) WILL BE PRINTED')
 9006 FORMAT('0',3X,'REPLICATION NUMBER ',I3,' OF ',I3,
     1       ' REPLICATIONS')
 9007 FORMAT('0',3X,'RANDOM PAIRING WILL BE USED')
 9008 FORMAT('0',3X,'A PLOT FILE OF SAMPLED VALUES AND CUMULATIVE ',
     1       'PROBABILITIES WILL BE CREATED')
      END
      SUBROUTINE BETA(J)
C***********************************************************************
C SUBROUTINE BETA IS USED TO GENERATE A BETA DISTRIBUTION ON THE
C INTERVAL (A,B) AND WITH PARAMETERS P AND Q
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      LOGICAL  FAILURE 
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON /PQ/P, Q, NZ
      REAL*4     P, Q
      REAL*4     RNUM
      INTEGER*4  NZ

      REAL*4 CON, PROBINC, A, B, STRTPT, BX, DEL, PROB
      INTEGER*4 LOC, I, J, I10, I0
      EXTERNAL BETAFN
      LOC(I,J) = (J-1)*N+I
      I10 = 10
      I0  = 0
      CALL ERXSET(I10,I0)
      CON = 1.E-06
      PROBINC = 1./FLOAT(N)
      IF(IRS.NE.0)PROBINC=1.0
      READ (8)A,B,P,Q
      STRTPT = 0.
      BX = P/(P+Q)
      DEL = SQRT((P*Q)/((P+Q)**2*(P+Q+1)))
      DO 1 I = 1,N
        RNUM=RAN(ISEED)
	PROB = PROBINC*RNUM+ STRTPT
        CALL INVBETA(PROB,BX,P,Q,CON,CON,FAILURE)
        IF (FAILURE) PRINT '(A,F8.4,A,F10.4,A,F10.4)', 
     1          ' INVBETA FAILURE FOR PROB=',PROB,'   P=',P,'   Q=',Q
	X(LOC(I,J)) = A + (B-A)*BX
	IF(IRS.EQ.0)STRTPT = STRTPT + PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), PROB, RNUM
    1 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE BETAFN(X,FOFX)
C***********************************************************************
C SUBROUTINE BETAFN IS USED IN GENERATING A BETA DISTRIBUTION
C***********************************************************************
      COMMON /PQ/P, Q, NZ
      REAL*4     P, Q
      INTEGER*4  NZ

      REAL*4 X, FOFX, Y(1), OMX
      INTEGER*4 N
      IF(X.LT.0.)X=0.
      IF(X.GT.1.)X=1.
      OMX = 1.-X
      N=1
      CALL BETAIC(X,OMX,P,Q,N,Y,NZ)
      FOFX = Y(1)
      RETURN
      END
C-----------------------------------------------------------------------
      REAL*4 FUNCTION BETAHKH(K)
C   EVALUATE BETA(K/2,1/2)
C   REVISION LOG:
C    4/30/93  HJI  ORIGINAL VERSION.

      IMPLICIT NONE
      INTEGER*4  I, II, K
      REAL*8   R

      IF (MOD(K,2) .EQ. 0) THEN
	 R=2.0D0
	 II=4
      ELSE
	 R=3.141592653589793
	 II=3
      ENDIF
      DO 10 I=II,K,2
   10 R=R*(I-2.)/(I-1.)
      BETAHKH=R
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BETAIC(X,OMX,A,B,N,Y,NZ)
C***********************************************************************
C SUBROUTINE BETAIC IS USED IN GENERATING A BETA DISTRIBUTION
C
C     WRITTEN BY D.E. AMOS AND S.L. DANIEL, JANUARY, 1975.
C
C     REFERENCE SC-DR-69 591
C
C     ABSTRACT
C         BETAIC COMPUTES AN N MEMBER SEQUENCE OF BETA DISTRIBUTIONS
C         Y(K)=I(A+K-1,B,X), K=1,...,N , A.GT.0 , B.GT.0 , AND
C         0.LE.X.LE.1, WHERE I(A,B,X) IS THE INCOMPLETE BETA FUNCTION
C         NORMALIZED TO 1. AT X=1. THE RELATION OF THE INCOMPLETE BETA
C         FUNCTION TO THE GAUSS HYPERGEOMETRIC FUNCTION IS USED OVER
C         VARIOUS PARAMETER RANGES WITH SERIES OR ASYMPTOTIC EXPRESSIONS
C         USED FOR EVALUATION STARTING AT A+N-1 AND BO.GT.0. WITH
C         BO=B-INTEGER PART OF B OR 1. THEN A COMBINATION OF FORWARD
C         RECURSION ON THE PARAMETER TO RAISE BO TO B FOLLOWED BY
C         BACKWARD RECURSION ON THE FIRST PARAMETER TO DECREASE A+N-1 TO
C         A GETS THE REQUIRED SEQUENCE. I(A,B,X) SATISFIES A TWO-TERM
C         RELATION IN BOTH PARAMETERS WHERE ADDITIONS CAN BE USED
C         EXCLUSIVELY TO RETAIN SIGNIFICANT DIGITS. BOTH X AND OMX=1.-X
C         ARE ENTERED IN THE CALL LIST TO AVOID LOSSES OF SIGNIFICANCE
C         IN OMX WHEN AN ANALYTICAL EXPRESSION IS AVAILABLE( SEE THE F
C         AND T DISTRIBUTIONS). BETAIC USES HYPGEO AND BETALN.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C           X      - ARGUMENT, 0.LE.X.LE.1.
C           OMX    - 1.-X
C           A      - START VALUE OF FIRST PARAMETER, A.GT.0.
C           B      - VALUE OF SECOND PARAMETER,      B.GT.0.
C           N      - NUMBER OF BETA FUNCTIONS IN THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           Y      - A VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    Y(K)=I(A+K-1,B,X), K=1,...,N.
C           NZ     - UNDERFLOW FLAG
C                    NZ.EQ.0, A NORMAL RETURN.
C                    NZ.NE.0, UNDERFLOW, Y(K)=0.0, K=N-NZ+1,N RETURNED
C
C     ERROR CONDITIONS
C         IMPROPER INPUT - A FATAL ERROR
C         UNDERFLOW - A NON-FATAL ERROR.
C
C
C
C     BETAIC USES SUBROUTINES BETALN, GAMALN, HYPGEO, ERRCHK, ERRGET,
C                             ERRPRT, ERXSET, ERSTGT
C***********************************************************************
      REAL*4 Y(1), GBZ(23), PP(35), GCOE(66), WCOE(49), GCOE1(59)
      REAL*4 GCOE2(7), X, OMX, A, B, E1, E2, ALB, R1MACH, ELIM, UR
      REAL*4 FLIM, P9LIM, REL2, ALIM, TOL, AA, BB, BO, DA, DB, XLN
      REAL*4 DX, W, OMXLN, SW, A1, P, BP1, T2, T3, T4, COE, BETALN
      REAL*4 HYPGEO, SBP1, OMSA, E, PHY, C, F, RBA, SAPB, RSAPB
      REAL*4 RSAPB2, SUMK, AK, ASUM, BSUM, GAMALN, A2, ATEST, SAM1
      REAL*4 ZZ, RNP1, SN, RN, RSAVE, RTEST, R1, AJ, RA1, DAK, DKBO
      REAL*4 RAT, S, FF, SUM, S1, TERM, AM1, BKK, DAM1, DBM2
      INTEGER*4 N, NZ, I1MACH, MAXR, NSAVE, I, IB, IDX, K, J, IN1
      INTEGER*4 IA, NN, KK, JJ, IN, INP1, KFLAG, NBAR, K1, K2, IS, KL
      EQUIVALENCE(GCOE(1),GCOE1(1)),(GCOE(60),GCOE2(1))
      DATA  WCOE           /-6.41025641025641E-03,-8.33333333333333E-02,
     1-5.00000000000000E-01,-1.46538461538462E+01,-4.58333333333333E+00,
     2-8.25000000000000E+00,-1.10000000000000E+01,-1.10000000000000E+01,
     3-8.25000000000000E+00,-4.58333333333333E+00,-1.83333333333333E+00,
     4-5.00000000000000E-01,-8.33333333333333E-02, 1.91752691752692E-03,
     5 2.10927960927961E-02, 8.72474747474747E-01, 3.16391941391941E-01,
     6 6.32783882783883E-01, 8.85897435897436E-01, 6.63847818847819E+00,
     7 6.32783882783883E-01, 3.16391941391941E-01, 1.05463980463980E-01,
     8 2.10927960927961E-02,-8.41750841750842E-04,-7.57575757575758E-03,
     9-3.03030303030303E-02,-7.07070707070707E-02,-4.31481481481481E+00,
     A-1.06060606060606E-01,-7.07070707070707E-02,-3.03030303030303E-02,
     B-7.57575757575758E-03, 5.95238095238095E-04, 4.16666666666667E-03,
     C 4.29166666666667E-01, 2.08333333333333E-02, 2.08333333333333E-02,
     D 1.25000000000000E-02, 4.16666666666667E-03,-7.93650793650794E-04,
     E-3.96825396825397E-03,-7.93650793650794E-03,-7.93650793650794E-03,
     F-3.96825396825397E-03, 2.77777777777778E-03, 8.33333333333333E-03,
     G 8.33333333333333E-03,-8.33333333333333E-02/
      DATA  GCOE1          /-8.33333333333333E-02, 1.55375180375180E-01,
     1-1.08042328042328E-01, 3.81997721928277E-02,-7.78250845091123E-03,
     2 9.70630021555947E-04,-7.60017468278464E-05, 3.72597920340976E-06,
     3-1.10336922674162E-07, 1.79409630364491E-09,-1.22324747975790E-11,
     4-9.09090909090909E-02, 1.60747354497354E-01,-1.03988095238095E-01,
     5 3.35006062610229E-02,-6.06796859690378E-03, 6.52046360596708E-04,
     6-4.21167695473251E-05, 1.59315751763668E-06,-3.22937334656085E-08,
     7 2.69114445546737E-10,-1.00000000000000E-01, 1.66269841269841E-01,
     8-9.87979497354497E-02, 2.84777336860670E-02,-4.46727109053498E-03,
     9 3.97698045267490E-04,-1.98929398148148E-05, 5.16699735449735E-07,
     A-5.38228891093474E-09,-1.11111111111111E-01, 1.71785714285714E-01,
     B-9.21626984126984E-02, 2.31706532921811E-02,-3.02854938271605E-03,
     C 2.09780092592593E-04,-7.23379629629630E-06, 9.68812003968254E-08,
     D-1.25000000000000E-01, 1.76984126984127E-01,-8.36805555555556E-02,
     E 1.76697530864198E-02,-1.80844907407407E-03, 8.68055555555556E-05,
     F-1.55009920634921E-06,-1.42857142857143E-01, 1.81245000000000E-01,
     G-7.28395061728395E-02, 1.21527777777778E-02,-8.68055555555556E-04,
     H 2.17013888888889E-05,-1.66666666666667E-01, 1.83333333333333E-01,
     I-5.90277777777778E-02, 6.94444444444444E-03,-2.60416666666667E-04,
     J-2.00000000000000E-01, 1.80555555555556E-01,-4.16666666666667E-02/
      DATA  GCOE2          / 2.60416666666667E-03,-2.50000000000000E-01,
     1 1.66666666666667E-01,-2.08333333333333E-02,-3.33333333333333E-01,
     2 1.25000000000000E-01,-5.00000000000000E-01/
C
      E1=ABS(I1MACH(12))
      E2=ABS(I1MACH(13))
      ALB=R1MACH(5)
      ELIM=(MIN(E1,E2)*ALB-3.0)*2.303E0
      UR=MAX(1.0E-14,R1MACH(4))
      FLIM=1.0/(UR*1.0E3)
      P9LIM=1.0-1.0E3*UR
      REL2=UR*1.0E4
      ALIM=UR
      TOL=UR*5.0E3
      MAXR=1000
C
C
C     TESTING OF VARIABLES,  ABS(X+OMX-1).LT.TOL, 0.LE.X.LE.1,
C     A.GT.0,  B.GT.0,  N.GE.1
C
      IF(N.LE.0) GO TO 90
      IF(A.LE.0.0) GO TO 91
      IF(B.LE.0.0) GO TO 92
      NZ=0
      AA=A+FLOAT(N)-1.
      BB=B
      NSAVE=N
      IF (X) 93,25,75
25    DO 30 I=1,N
   30 Y(I)=0.
      RETURN
35    DO 40 I=1,N
   40 Y(I)=1.
      RETURN
75    IF (X.EQ.1.) GO TO 35
      IF(OMX.LT.0.0) GO TO 93
      IF(ABS(X+OMX-1.).GT.TOL) GO TO 93
C
C     COMPUTATION OF BO
C
      IB=BB
      BO=BB-FLOAT(IB)
      IF (BO) 80,80,85
80    BO=1.0
      IB=IB-1
85    DA=AA
      DB=BB
C
C     COMPUTATION OF XLN AND OMXLN BY FORMULA (8)
C
      IF(X.GT.0.9) GO TO 95
      XLN=LOG(X)
      GO TO 105
95    ASSIGN 100 TO IDX
      DX=OMX
      GO TO 120
100   XLN=W
  105 IF(OMX.GT.0.9) GO TO 115
      OMXLN=LOG(OMX)
      GO TO 140
115   ASSIGN 135 TO IDX
      DX=X
C
C     LOGARITHM ROUTINE FOR ARGUMENTS GT 0.9 BY FORMULA (8)
C
120   SW=0.
      A1=DX
      DO 125 K=1,16
      GBZ(K)=A1
125   A1=A1*DX
      K=17
      DO 130 J=1,16
      K=K-1
      SW=SW+GBZ(K)/FLOAT(K)
130   CONTINUE
      W=-SW
      GO TO IDX,(100,135)
135   OMXLN=W
  140 IF(BB.NE.1.0) GO TO 150
C
C     BY FORMULA(1A)
C
      P=EXP(AA*XLN)
      GO TO 395
  150 BP1=B+1.
      IF(BP1*X.LT.(B+.7)/BP1) GO TO 206
      IF(BO.EQ.1.) GO TO 175
      IF(X.GT.0.7) GO TO 215
      IF(X.GT.0.4) GO TO 180
      IF(AA.LE.BB) GO TO 180
C
C     BY FORMULA(4)
C
      T2=1.-DB
      T3=DA+1.
      T4=-X/OMX
      COE=EXP(DA*XLN-T2*OMXLN-BETALN(DA,DB))/DA
      P=COE*HYPGEO(1.,T2,T3,T4,1.)
      GO TO 395
C
C     BY FORMULA(1A)
C
175   P=EXP(AA*XLN)
      GO TO 360
  180 IF(DB.GT.1.) GO TO 200
C
C     BY FORMULA(1)
C
      ASSIGN 395 TO IN1
      GO TO 210
C
  200 ASSIGN 359 TO IN1
      DB=BO
      GO TO 210
  206 ASSIGN 395 TO IN1
C
C     BY FORMULA(1)
C
210   T2=1.-DB
      T3=DA+1.
      COE= EXP(DA*XLN-BETALN(DA,DB))/DA
      P=COE*HYPGEO(DA,T2,T3,X,1.)
      GO TO IN1,(395,359)
  215 IF(AA.GT.BB) GO TO 225
C
C     BY FORMULA(3)
C
      T2=1.-AA
      T3=BB+1.
      T4=-OMX/X
      COE= EXP((AA-1.)*XLN+BB*OMXLN-BETALN(AA,BB))/BB
      P=1.-COE*HYPGEO(1.,T2,T3,T4,1.)
      IF(P.GE.0.1) GO TO 395
C
225   CONTINUE
      IF(AA*OMX.GT.3.) GO TO 280
      IF(AA.LT.20.) GO TO 240
      DA=AA
      GO TO 245
240   IA=AA
      NN=20-IA
      DA=AA+FLOAT(NN+1)
      N=N+NN+1
      IF (DA*OMX.GT.3.) GO TO 280
245   CONTINUE
C
C     P FOR DA*OMX.LE.3. AND DA.GE.20. BY FORMULA (2),(2A),(2B),(2C)
C
      SBP1=BO+1.
      OMSA=1.-DA
      E=HYPGEO(BO,OMSA,SBP1,OMX,0.0)
      PHY=1.+E
      C= EXP(BO*OMXLN-BETALN(DA,BO))/BO
      P=1.-C*PHY
      IF(ABS(P).GE.0.1) GO TO 275
      F=C-1.
      IF(ABS(F).GE.0.1) GO TO 270
C
C     COMPUTE W
C
      RBA=BO/DA
      SAPB=DA+BO
      RSAPB=1./SAPB
      RSAPB2=RSAPB*RSAPB
      SUMK=0.
      AK=1.
      A1=-1.
      DO 260 K=1,13
      A1=-A1*RBA
      SUMK=SUMK+A1/AK
      AK=AK+1.
  260 CONTINUE
      ASUM=0.
      KK=13
      J=1
      DO 261 JJ=1,13,2
      BSUM=0.
      DO 262 K=1,KK
      BSUM=RBA*(BSUM+WCOE(J))
      J=J+1
  262 CONTINUE
      ASUM=RSAPB2*ASUM+BSUM
      KK=KK-2
  261 CONTINUE
      ASUM=RSAPB*ASUM
      ASUM=ASUM+(DA-.5)*SUMK+BO*LOG(SAPB*OMX)-BO-GAMALN(SBP1,0.)
      F=0.
      A2=1.
      AK=1.
      ATEST=ALIM*ABS(ASUM)
      DO 265 K=1,25
      A2=A2*ASUM/AK
      F=F+A2
      IF (ABS(A2).LT.ATEST) GO TO 270
      AK=AK+1.
265   CONTINUE
270   P=-F-E-F*E
  275 IF(BB.LT.1.0) GO TO 395
      GO TO 360
  280 IF(AA*OMX.GT.35.) GO TO 340
      IF(AA.LT.1001.) GO TO 295
      DA=AA
      N=NSAVE
      GO TO 300
  295 IN=1001.-AA
      INP1=IN+1
      N=NSAVE+INP1
      DA=AA+FLOAT(INP1)
      IF (DA*OMX.LT.35.) GO TO 300
      IN=35./OMX-AA+1.
      INP1=IN+1
      DA=AA+FLOAT(INP1)
      N=NSAVE+INP1
      GO TO 345
300   CONTINUE
      SAM1=DA-1.
      ZZ=SAM1*OMX
C
C     BIG GAMMA(BO,ZZ), 3.LT.ZZ.LT.35., ZZ=(DA-1)*OMX BY FORMULA (5)
C
      KFLAG=0
      NBAR=20
  305 IF(NBAR.GT.MAXR) GO TO 94
      RNP1=0.0
      SN=NBAR+1
      DO 310 J=1,NBAR
      SN=SN-1.
      RN=ZZ+(SN-BO)*RNP1/(SN+RNP1)
310   RNP1=RN
      IF(KFLAG.NE.0) GO TO 320
      KFLAG=1
      NBAR=NBAR+NBAR
      RSAVE=RN
      GO TO 305
320   RTEST=ABS((RSAVE-RN)/RN)
      IF(RTEST.LE.REL2) GO TO 330
      RSAVE=RN
      NBAR=NBAR+NBAR
      GO TO 305
330   R1= EXP(BO*LOG(ZZ)-ZZ)
      GBZ(1)=R1/RN
C
C     RECURSION FOR BIG GAMMA(BO+K,ZZ), K=1,2,...23 BY FORMULA (6)
C
      AJ=0.
      DO 335 K=2,23
      GBZ(K)=(AJ+BO)*GBZ(K-1)+R1
      AJ=AJ+1.
      R1=R1*ZZ
335   CONTINUE
C
C     P FOR 0.965.LE.X.LT.0.997 AND DA.GE.1001 BY FORMULA (7)
C
      RA1=1./SAM1
      ASUM=0.
      K1=13
      K2=23
      J=1
      DO 336 KK=1,11
      BSUM=0.
      DO 337 K=K1,K2
      BSUM=BSUM+GCOE(J)*GBZ(K)
      J=J+1
  337 CONTINUE
      K1=K1-1
      K2=K2-2
      ASUM=RA1*(ASUM+BSUM)
  336 CONTINUE
      ASUM=ASUM+GBZ(1)
C
      COE =  EXP(-BETALN(DA,BO) - BO*LOG(SAM1))
      P=COE*ASUM
      GO TO 275
C
C     P FOR 0.7.LT.X.LT.0.965 AND DA*OMX.GT.35. BY FORMULA (4A)
C
340   DA=AA
345   CONTINUE
      A1=1.
      DAK=DA+1.
      DKBO=1.-BO
      RAT=-X/OMX
      DO 350 K=1,35
      A1=(RAT*DKBO/DAK)*A1
      PP(K)=A1
      DAK=DAK+1.
      DKBO=DKBO+1.
350   CONTINUE
      S=0.
      J=36
      DO 355 K=1,35
      J=J-1
355   S=S+PP(J)
      S=S+1.
      COE= EXP(DA*XLN+(BO-1.)*OMXLN-BETALN(DA,BO))/DA
      P=S*COE
      GO TO 275
C
C     FORWARD RECURRENCE BY FORMULA (12A)
C
  359 DB=BB
360   IF (IB.LE.0) GO TO 395
      S=DA*XLN+BO*OMXLN-BETALN(DA,BO)-LOG(BO)
C
C     SCALING FOR UNDERFLOW
C
      IF(S.GE.-ELIM .AND. P.GT.0.0) GO TO 385
  366 FF=1.
      P=0.
      SUM=0.
375   IF (IB.EQ.0) GO TO 395
      IF(FF.GE.FLIM) GO TO 380
      SUM=SUM+FF
      FF=FF*((BO+DA)/(BO+1.))*OMX
      IB=IB-1
      BO=BO+1.
      GO TO 375
380   S1=S
      S=S+LOG(FF)
      IF(S.LT.-ELIM) GO TO 366
      ZZ=S1+LOG(SUM)
      IF(ZZ.LT.-ELIM) GO TO 366
      P=EXP(ZZ)
  385 TERM=EXP(S)
      P=P+TERM
      IB=IB-1
      IF (IB.EQ.0) GO TO 395
      AM1=DA-1.
      BKK=BO
      DO 390 KK=1,IB
      BKK=BKK+1.
      TERM=TERM*((BKK+AM1)/BKK)*OMX
      P=P+TERM
  390 CONTINUE
395   CONTINUE
      IF(N.NE.1) GO TO 405
      Y(1)=P
      IF(P.EQ.0.0) NZ=1
      GO TO 440
C
C     BACKWARD RECURRENCE BY FORMULA (9A)
C
405   DAM1=DA-1.
      DBM2=DB-2.
      S=DAM1*XLN+DB*OMXLN-BETALN(DA,DB)-LOG(DB+DAM1)
C
C     SCALING FOR UNDERFLOW
C
      IF(S.GE.-ELIM .AND. P.GT.0.0) GO TO 430
      IF(P.LT.P9LIM) GO TO 411
      DO 406 I=1,NSAVE
  406 Y(I)=1.
  407 N=NSAVE
      RETURN
  411 FF=1.
      P=0.
      SUM=0.
  420 CONTINUE
      IF(N.EQ.1) GO TO 412
      IF(FF.GE.FLIM) GO TO 425
      SUM=SUM+FF
      FF=FF*(DA-1.)/((DA+DBM2)*X)
      DA=DA-1.
      N=N-1
      GO TO 420
  412 DO 413 I=1,NSAVE
  413 Y(I)=0.
      N =NSAVE
      NZ=N
      RETURN
425   S1=S
      S=S+LOG(FF)
      IF(S.LT.-ELIM) GO TO 411
      ZZ=S1+LOG(SUM)
      IF(ZZ.LT.-ELIM) GO TO 411
      P=EXP(ZZ)
  430 TERM=EXP(S)
      IF(N-NSAVE) 431,434,432
  431 IS=N+1
      DO 429 K=IS,NSAVE
  429 Y(K)=0.
      NZ=NSAVE-N
      Y(N)=P
      IF(N-2) 440,408,409
  408 Y(1)=P+TERM
      GO TO 440
  409 KL=N-1
      GO TO 435
  434 Y(NSAVE)=P
      KL=NSAVE-1
      IF(KL.GT.1) GO TO 435
  441 Y(1)=P+TERM
      GO TO 440
  432 KL=N-NSAVE
      DO 433 K=1,KL
      P=P+TERM
      TERM=TERM*(DA-1.)/((DA+DBM2)*X)
      DA=DA-1.
  433 CONTINUE
      Y(NSAVE)=P
      KL=NSAVE-1
      IF(KL-1) 440,441,435
  435 P=P+TERM
      Y(KL)=P
      KL=KL-1
      KK=KL
      DO 436 K=1,KL
      TERM=TERM*(DA-1.)/((DA+DBM2)*X)
      DA=DA-1.
      Y(KK)=TERM+Y(KK+1)
  436 KK=KK-1
440   N=NSAVE
      DO 442 I=1,N
      IF(Y(I).GT.1.) Y(I)=1.
  442 CONTINUE
      RETURN
C
   90 CALL ERRCHK('IN BETAIC, IMPROPER INPUT FOR N.')
      RETURN
   91 CALL ERRCHK('IN BETAIC, IMPROPER INPUT FOR A.')
      RETURN
   92 CALL ERRCHK('IN BETAIC, IMPROPER INPUT FOR B.')
      RETURN
   93 CALL ERRCHK('IN BETAIC, IMPROPER INPUT FOR X OR OMX.')
      RETURN
   94 CALL ERRCHK(
     1  'IN BETAIC, NO CONVERGENCE IN CONTINUED FRACTION FOR GAMMA'//
     2  ' FUNCTION.')
      RETURN
      END
C=======================================================================
      FUNCTION BETALN(A,B)
C***********************************************************************
C
C FUNCTION BETALN IS USED IN GENERATING A BETA DISTRIBUTION BY
C COMPUTING THE NATURAL LOG OF BETA FUNCTION 
C
C   REVISION LOG:
C    2/07/97  HJI  INITIAL VERSION.
C
C***********************************************************************
C
      IMPLICIT NONE
C
C   DECLARATION FOR EXTERNAL FUNCTION
C
      REAL*4 BETALN, GAMALN 
C
      REAL*4 A,B
C
C   THE DUMMY SECOND ARGUMENT IS FOR COMPATIBILITY WITH DON AMOS' GAMALN
C   ROUTINE AS USED IN LHS.
C
      BETALN = GAMALN(A,0.) + GAMALN(B,0.) - GAMALN(A+B,0.)
C
      RETURN
      END 
C
C=======================================================================
      FUNCTION GAMALN(XX,YY)
C***********************************************************************
C
C   FUNCTION GAMALN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION 
C   SOURCE:  NUMERICAL RECIPES, PRESS, FLANNERY, ET AL, PAGE 157. 
C
C   REVISION LOG:  
C   12/18/89  HJI  ORIGINAL VERSION
C   THE DUMMY SECOND ARGUMENT YY IS FOR COMPATIBILITY WITH DON AMOS' 
C   GAMALN ROUTINE AS USED IN LHS.
C
C***********************************************************************
C
      INTEGER*4   I
      REAL*4      XX, YY
      REAL*8      COEF(6), STP, SUM, TMP, X
C  
      DATA   COEF/  76.18009173,  -86.50532033,     24.01409822, 
     1              -1.231739516,   0.120858003D-2, -0.536382D-5 /
C
      DATA   STP/   2.50662827465/ 
C
      X=XX-1.0D0
      TMP=X+5.5D0
      TMP=(X+0.5)*DLOG(TMP)-TMP      
      SUM=1.0D0
C
      DO 10 I=1,6
      X=X+1.
      SUM=SUM+COEF(I)/X
   10 CONTINUE
C
      GAMALN=TMP+DLOG(SUM*STP)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE CHKDAT(PAR,A,B,P,Q)
C***********************************************************************
C SUBROUTINE CHKDAT CHECKS DISTRIBUTION PARAMETERS FOR CONSISTENCY
C***********************************************************************
      REAL*4 A, B, P, Q
      CHARACTER PAR*(*), PLOG*3, PBETA*4
      PARAMETER (PLOG='LOG',PBETA='BETA')
      IF(A.GE.B)THEN
	WRITE(6,9001)PAR,A,B
	RETURN
      ELSE IF(PAR(1:3).EQ.PLOG.AND.(A.LE.0.0 .OR. B.LE.0.0))THEN
	WRITE(6,9002)PAR,A,B
	RETURN
      ELSE IF(PAR.EQ.PBETA.AND.(P.LT.0.0 .OR. Q.LT.0.0))THEN
	WRITE(6,9003)P,Q
	RETURN
      ELSE IF(PAR.EQ.PBETA)THEN
	WRITE(8)A,B,P,Q
      ELSE
	WRITE(8)A,B
      ENDIF
      RETURN
 9001 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE LOWER LIMIT A ',
     1       G20.10,/,6X,'IS GREATER THAN OR EQUAL TO THE UPPER ',
     2       'LIMIT B ',G20.10)
 9002 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE LOWER LIMIT A ',
     1       G20.10,/,6X,'OR THE UPPER LIMIT B ',G20.10,/,6X,
     2       'SHOULD NOT BE LESS THAN OR EQUAL TO ZERO')
 9003 FORMAT('1',5X,'FOR THE BETA DISTRIBUTION THE PARAMETER P ',
     1       G20.10,/,6X,'OR THE PARAMETER Q ',G20.10,/,6X,'SHOULD ',
     2       'NOT BE LESS THAN 0.0 OR NUMERICAL INSTABILITIES WILL ',
     3       'OCCUR AS',/,6X,'A RESULT OF TRYING TO INVERT THE BETA ',
     4       'DISTRIBUTION ON PORTIONS OF THE CURVE',/,6X,'WITH ',
     5       'EXTREME SLOPES')
      END
C=======================================================================
      SUBROUTINE CHKDIM(IOPT,NREQ,NMX,PCHR,LCHR)
C***********************************************************************
C SUBROUTINE CHKDIM CHECKS TO SEE IF PROGRAM DIMENSIONS HAVE BEEN
C EXCEEDED
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (LENC=80)
      INTEGER*4 IOPT, NREQ, NMX, IERR, I
      CHARACTER*(*) PCHR, LCHR
      COMMON/CHRCRD/CRDSTR
      CHARACTER     CRDSTR(NVAR)*(LENC)
      COMMON/OBSTR/NSTR, NOBSTR
      INTEGER*4    NSTR, NOBSTR(NVAR)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      IF(IOPT.EQ.2)GO TO 100
      IF(NREQ.GT.NMX)THEN
	WRITE(6,9001)PCHR,NREQ,LCHR,NMX,LCHR,LCHR
	RETURN
      ENDIF
      RETURN
  100 CONTINUE
      IF(NREQ.GT.NMX)THEN
	WRITE(6,9002)NREQ,NMX
	RETURN
      ENDIF
      IERR=0
      DO 200 I=1,NSTR
	IF(NOBSTR(I).NE.N)THEN
	  WRITE(6,9003)NOBSTR(I),N,CRDSTR(I)
	  IERR=1
	ENDIF
  200 CONTINUE
      IF(IERR.EQ.1)RETURN
      RETURN
 9001 FORMAT('1',5X,'THE PARAMETER CARD ',A,'REQUESTED ',I4,1X,A,/,
     1       6X,'ONLY ',I4,1X,A,' ARE CURRENTLY PERMITTED',/,6X,
     2       'PLEASE CONSULT THE USER MANUAL FOR INSTRUCTIONS ON ',
     3       'HOW TO ALLOW MORE ',A)
 9002 FORMAT('1',5X,'THE NUMBER OF VARIABLES REQUESTED ',I3,/,6X,
     1       ' EXCEEDS THE MAXIMUM NUMBER OF VARIABLES CURRENTLY ',
     2       'PERMITTED ',I3,/,6X,'PLEASE CONSULT THE USER MANUAL FOR',
     3       'INSTRUCTIONS ON HOW TO ALLOW MORE VARIABLES')
 9003 FORMAT('1',5X,'THE FOLLOWING DISTRIBUTION CARD REQUESTED ',I4,
     1       ' OBSERVATIONS',/,6X,'HOWEVER THE NOBS PARAMETER CARD ',
     2       'REQUESTED ',I4,' OBSERVATIONS',/,6X,'THIS DISCREPANCY ',
     3       'MUST BE RESOLVED BEFORE PROCESSING CAN CONTINUE',//,3X,
     4       '***',A,'***')
      END
      SUBROUTINE CHKEXP(A,B,C)
C***********************************************************************
C SUBROUTINE CHKEXP CHECKS PARAMETERS OF THE EXPONENTIAL
C DISTRIBUTION FOR CONSISTENCY
C***********************************************************************
      REAL*4 A, B, C
      IF(A.GE.B)THEN
	WRITE(6,9010)A,B,C
	RETURN
      ELSE
	WRITE(8)A,B,C
      ENDIF
      RETURN
 9010 FORMAT('1',5X,'FOR THE EXPONENTIAL DISTRIBUTION THE PARAMETERS ',
     1       'A,B ',2G20.10,/,6X,'HAVE BEEN INPUT IN THE INCORRECT ',
     2       'ORDER',/,6X,'THE INTERVAL START (A) MUST BE LESS THAN ',
     3       'THE INTERVAL END (B)')
      END
      SUBROUTINE CHKRAY(A,B,C)
C***********************************************************************
C SUBROUTINE CHKRAY CHECKS PARAMETERS OF THE RAYLEIGH
C DISTRIBUTION FOR CONSISTENCY
C***********************************************************************
      REAL*4 A, B, C
      IF(A.GE.B)THEN
	WRITE(6,9010)A,B
	RETURN
      ELSE
	WRITE(8)A,B,C
      ENDIF
      RETURN
 9010 FORMAT('1',5X,'FOR THE RAYLEIGH DISTRIBUTION THE PARAMETERS ',
     1       'A,B ',2G20.10/6X,'HAVE BEEN INPUT IN THE INCORRECT ',
     2       'ORDER',/,6X,'THE INTERVAL START (A) MUST BE LESS THAN ',
     3       'THE INTERVAL END (B)')
      END
      SUBROUTINE CHKREX(A,B,C,D)
C***********************************************************************
C SUBROUTINE CHKREX CHECKS PARAMETERS OF THE RAYLEIGH-EXPONENTIAL
C DISTRIBUTION FOR CONSISTENCY
C***********************************************************************
      REAL*4 A, B, C, D
      IF(A .GE. B .OR. B .GE. C) THEN
	WRITE(6,9010)A,B,C
	RETURN
      ELSE
	WRITE(8)A,B,C,D
      ENDIF
      RETURN
 9010 FORMAT('1',5X,
     1  'FOR THE RAYLEIGH-EXPONENTIAL DISTRIBUTION, THE PARAMETERS '/
     2  5X,'A,B,C',3G20.10/6X,'HAVE BEEN INPUT IN AN INCORRECT ORDER.'/
     3  5X,'THE INTERVAL START (A) MUST BE LESS THAN THE MIDPOINT (B)'/
     4  5X,'AND THE MIDPOINT (B) MUST BE LESS THAN THE ENDPOINT (C).')
      END
      SUBROUTINE CHKSTR(PAR,CARD)
C***********************************************************************
C SUBROUTINE CHKSTR CHECKS PARAMETERS OF THE UNIFORM* AND THE
C LOGUNIFORM* DISTRIBUTIONS FOR CONSISTENCY
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (NINTMX=50)
      PARAMETER (LENC=80)
      CHARACTER PAR*(*), CARD*(LENC), PLOG*3
      COMMON/STAR/NSUBOB, NINT, SUBINT
      INTEGER*4   NSUBOB(NINTMX), NINT
      REAL*4      SUBINT(NINTMX+1)
      COMMON/CHRCRD/CRDSTR
      CHARACTER     CRDSTR(NVAR)*(LENC)
      COMMON/OBSTR/NSTR, NOBSTR
      INTEGER*4    NSTR, NOBSTR(NVAR)
C
      INTEGER*4 I
      PARAMETER (PLOG='LOG')
C
      IF(NINT.EQ.0)THEN
	WRITE(6,9001)PAR
	RETURN
      ELSE IF(NINT.GT.NINTMX)THEN
	WRITE(6,9002)PAR,NINT,NINTMX
	RETURN
      ELSE
	WRITE(8)NINT
      ENDIF
      NSTR=NSTR+1
      DO 100 I=1,NINT
	IF(PAR(1:3).EQ.PLOG.AND.SUBINT(I).LE.0.0)THEN
	  WRITE(6,9003)PAR,I,SUBINT(I)
	  RETURN
	ELSE IF(SUBINT(I).GE.SUBINT(I+1))THEN
	  WRITE(6,9004)PAR,I,SUBINT(I),SUBINT(I+1)
	  RETURN
	ELSE IF(NSUBOB(I).LT.0)THEN
	  WRITE(6,9005)PAR,I
	  RETURN
	ELSE
	  NOBSTR(NSTR)=NOBSTR(NSTR)+NSUBOB(I)
	  WRITE(8)NSUBOB(I),SUBINT(I),SUBINT(I+1)
	ENDIF
  100 CONTINUE
      CRDSTR(NSTR)=CARD
      RETURN
 9001 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE NUMBER OF ',
     1       'SUBINTERVALS IS ZERO')
 9002 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE NUMBER OF ',
     1       'SUBINTERVALS REQUESTED ',I3,/,6X,'IS GREATER THAN THE ',
     2       'MAXIMUM NUMBER OF SUBINTERVALS CURRENTLY PERMITTED ',I3,
     3       /,6X,'PLEASE CONSULT THE USER MANUAL FOR INSTRUCTIONS ',
     4       'ON HOW TO ALLOW MORE SUBINTERVALS')
 9003 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION THE SUBINTERVAL ',
     1       'LIMIT FOR SUBINTERVAL ',I3,/,6X,'IS LESS THAN OR ',
     2       'EQUAL TO ZERO ',G20.10)
 9004 FORMAT('1',5X,'ON THE ',A,'DISTRIBUTION FOR SUBINTERVAL ',
     1       I3,' THE LOWER LIMIT ',G20.10,/,6X,'IS GREATER ',
     2       'THAN OR EQUAL TO THE UPPER LIMIT ',G20.10)
 9005 FORMAT('1',5X,'FOR THE ',A,'DISTRIBUTION SUBINTERVAL ',I3,
     1       ' REQUESTED A NEGATIVE NUMBER OF OBSERVATIONS')
      END
      SUBROUTINE CHKTRI(A,B,C)
C***********************************************************************
C SUBROUTINE CHKTRI CHECKS PARAMETERS OF THE TRIANGULAR
C DISTRIBUTION FOR CONSISTENCY
C***********************************************************************
      REAL*4 A, B, C
      IF(A.GT.B.OR.B.GT.C.OR.A.EQ.C)THEN
	WRITE(6,9000)A,B,C
	RETURN
      ELSE
	WRITE(8)A,B,C
      ENDIF
      RETURN
 9000 FORMAT('1',5X,'FOR THE TRIANGULAR DISTRIBUTION THE PARAMETERS ',
     1       'A,B,C ',3G20.10,/,6X,'HAVE BEEN INPUT IN THE INCORRECT ',
     2       'ORDER',/,6X,'PLEASE CONSULT THE USER MANUAL FOR ',
     3       'INSTRUCTIONS ON THE ORDER OF THE PARAMETERS')
      END
      SUBROUTINE CHKZRO(N,NV,IRSET)
C***********************************************************************
C SUBROUTINE CHKZRO CHECKS TO MAKE SURE THAT THE MINIMUM
C REQUIREMENTS FOR A SAMPLE HAVE BEEN MET
C***********************************************************************
      INTEGER*4 N, NV, IRSET
      IF(N.EQ.0)THEN
	WRITE(6,9001)
	RETURN
      ELSE IF(NV.EQ.0)THEN
	WRITE(6,9002)
	RETURN
      ELSE IF(IRSET.EQ.0)THEN
	WRITE(6,9003)
	RETURN
      ENDIF
      RETURN
 9001 FORMAT('1',5X,'THE NUMBER OF OBSERVATIONS HAS NOT BEEN ',
     1       'SPECIFIED')
 9002 FORMAT('1',5X,'NO VARIABLES HAVE BEEN SPECIFIED')
 9003 FORMAT('1',5X,'A RANDOM SEED HAS NOT BEEN SPECIFIED')
      END
      SUBROUTINE CHLSKY
C***********************************************************************
C SUBROUTINE CHLSKY COMPUTES THE CHOLESKY FACTORIZATION OF A
C CORRELATION MATRIX STORED IN LOWER TRIANGULAR SYMMETRIC FORM
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/WORKC/Q, S
      REAL*4       Q((NVAR*(NVAR+1))/2),S((NVAR*(NVAR+1))/2)

      INTEGER*4 LOC1, I, J, NVX, I1, I2, I1MIN, IPLUS, K
      LOC1(I,J)=J+(I*I-I)/2
      NVX=(NV*(NV+1))/2
      DO 70 I1=1,NVX
   70 Q(I1)=0.0
      I1=1
      DO 80 I2=1,NV
   80 Q(LOC1(I2,I1))=CORR(LOC1(I2,I1))
   90 I1=I1+1
      I1MIN=I1-1
      DO 100 I2=1,I1MIN
  100 Q(LOC1(I1,I1))=Q(LOC1(I1,I1))+Q(LOC1(I1,I2))**2
      Q(LOC1(I1,I1))=SQRT(1.0-Q(LOC1(I1,I1)))
      IF (I1.GE.NV) GO TO 130
      IPLUS=I1+1
      DO 120 I2=IPLUS,NV
      DO 110 K=1,I1MIN
      Q(LOC1(I2,I1))=Q(LOC1(I2,I1))+Q(LOC1(I2,K))*Q(LOC1(I1,K))
  110 CONTINUE
  120 Q(LOC1(I2,I1))=(CORR(LOC1(I2,I1))-Q(LOC1(I2,I1)))/Q(LOC1(I1,I1))
      GO TO 90
  130 CONTINUE
      RETURN
      END
      SUBROUTINE CMCRD
C***********************************************************************
C SUBROUTINE CMCRD PROCESSES THE CORRELATION MATRIX PARAMETER CARD
C***********************************************************************
      PARAMETER (NVAR=100)
c**************************
      PARAMETER (NCVAR=1000)
c***************************
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/UICORR/ICVAR, JCVAR, CVAR, NCV
      INTEGER*4     ICVAR(NCVAR), JCVAR(NCVAR), NCV
      REAL*4        CVAR(NCVAR)

      REAL*4 RIJ(NCVAR*2)
      INTEGER*4 IJCVAR(2*NCVAR), LOC1, I, J, NCV2, ICV, NSIZE
      INTEGER*4 KCM, IM, JM
      LOC1(I,J)=J+(I*I-I)/2
      NCV2=NCV*2
      DO 100 ICV=1,NCV
	I=ICVAR(ICV)
	J=JCVAR(ICV)
	IF(ABS(CVAR(ICV)).GE.1.0)THEN
	  WRITE(6,9001)I,J,CVAR(ICV)
	  RETURN
	ELSE IF(I.EQ.J.AND.CVAR(ICV).NE.1.0)THEN
	  WRITE(6,9002)I,J,CVAR(ICV)
	  RETURN
	ELSE IF(I.GT.NV.OR.J.GT.NV)THEN
	  WRITE(6,9003)I,J,CVAR(ICV),NV
	  RETURN
	ELSE
	  IJCVAR(ICV)=I
	  IJCVAR(NCV+ICV)=J
	ENDIF
  100 CONTINUE
      DO 110 I=1,NCV2
  110 RIJ(I)=IJCVAR(I)
      CALL SIFT(RIJ,NCV2)
      DO 120 I=1,NCV2
  120 IJCVAR(I)=RIJ(I)
      NCM=1
      LCM(NCM)=IJCVAR(1)
      DO 200 I=2,NCV2
	IF(IJCVAR(I).NE.LCM(NCM))THEN
	  NCM=NCM+1
	  LCM(NCM)=IJCVAR(I)
	ENDIF
  200 CONTINUE
      NSIZE=(NCM*(NCM+1))/2
      DO 300 I=1,NSIZE
  300 CORR(I)=0.0
      DO 400 I=1,NCM
  400 CORR(LOC1(I,I))=1.0
      DO 600 ICV=1,NCV
	I=ICVAR(ICV)
	J=JCVAR(ICV)
	DO 500 KCM=1,NCM
	  IF(I.EQ.LCM(KCM))IM=KCM
	  IF(J.EQ.LCM(KCM))JM=KCM
  500   CONTINUE
	IF(IM.GT.JM)THEN
	  CORR(LOC1(IM,JM))=CVAR(ICV)
	ELSE
	  CORR(LOC1(JM,IM))=CVAR(ICV)
	ENDIF
  600 CONTINUE
      RETURN
 9001 FORMAT('1',3X,'THE CORRELATION BETWEEN VARIABLE ',I3,' AND ',
     1       'VARIABLE ',I3,' IS GREATER THAN ONE IN ABSOLUTE ',
     2       'VALUE ',F5.2)
 9002 FORMAT('1',3X,'THE CORRELATION BETWEEN VARIABLE ',I3,' AND ',
     1       'VARIABLE ',I3,' IS NOT EQUAL TO ONE ',F5.2)
 9003 FORMAT('1',3X,'THE CORRELATION BETWEEN VARIABLE ',I3,' AND ',
     2       'VARIABLE ',I3,' IS ',F5.2,/,4X,'HOWEVER ONLY ',I3,
     3       ' VARIABLES HAVE BEEN DEFINED')
      END
      SUBROUTINE CORCAL
C***********************************************************************
C SUBROUTINE CORCAL COMPUTES A CORRELATION MATRIX FROM THE SAMPLE
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM

      REAL*4 XM(NVAR), SSQ(NVAR), FN
      LOGICAL ERRS
      CHARACTER ERRMSG*128, ERMSG(2)*64
      INTEGER*4 LOC, I, J, LOC1, K, IL1
      EQUIVALENCE (ERRMSG, ERMSG(1))
      DATA ERMSG /
     1'0ERROR IN CORCAL - ALL VALUES FOR VARIABLE     ARE EQUAL (',
     2'    )    TO AVOID DIVISION BY ZERO EXECUTION IS BEING STOPPED.'
     3 /
      LOC(I,J)=(J-1)*N+I
      LOC1(I,J)=J+(I*I-I)/2
      ERRS = .False.
      IF(NV.EQ.1)THEN
	CORR(1)=1.0
	RETURN
      ENDIF
      DO 10 I=1,NV
	XM(I)=0.0
	SSQ(I)=0.0
   10 CONTINUE
      K=(NV*(NV+1))/2
      DO 20 I=1,K
   20 CORR(I)=0.0
      FN=N
C
C     COMPUTE THE MEAN
C
      DO 30 J=1,NV
      DO 30 K=1,N
   30 XM(J)=XM(J)+X(LOC(K,J))
      DO 40 J=1,NV
   40 XM(J)=XM(J)/FN
C
C     SUBTRACT THE MEAN FROM ALL OF THE ELEMENTS OF THE MATRIX
C     AND COMPUTE THE SUM OF SQUARES
C
      DO 50 J=1,NV
      DO 50 K=1,N
	X(LOC(K,J))=X(LOC(K,J))-XM(J)
	SSQ(J)=SSQ(J)+X(LOC(K,J))*X(LOC(K,J))
   50 CONTINUE
      DO 55 J = 1, NV
	IF (SSQ(J) .EQ. 0.0) THEN
          ERRS = .TRUE.
          WRITE (ERRMSG(44:46), '(I3)') J
          WRITE (ERRMSG(59:68), '(E10.4)') XM(J)
          CALL XERPRT(ERRMSG)
        ENDIF
   55 CONTINUE
      IF (ERRS) RETURN
C
C     COMPUTE THE CORRELATION
C
      DO 70 I=2,NV
      IL1=I-1
      DO 70 J=1,IL1
	DO 60 K=1,N
   60   CORR(LOC1(I,J))=CORR(LOC1(I,J))+X(LOC(K,I))*X(LOC(K,J))
   70 CONTINUE
      DO 80 I=2,NV
      IL1=I-1
      DO 80 J=1,IL1
   80 CORR(LOC1(I,J))=CORR(LOC1(I,J))/(SQRT(SSQ(I))*SQRT(SSQ(J)))
      DO 90 I=1,NV
   90 CORR(LOC1(I,I))=1.0
      RETURN
      END
      SUBROUTINE COROUT
C***********************************************************************
C SUBROUTINE COROUT WILL OUTPUT THE RAW AND RANK CORRELATION
C MATRICES OF THE SAMPLE IF REQUESTED
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      INTEGER*4 LOC, I, J
      LOC(I,J)=(J-1)*N+I
      REWIND 4
      READ(4)X
      DO 600 I=1,NV
  600 LCM(I)=I
      CALL CORCAL
      CALL PMTRX(NV,1)
      IF(N.GT.NV.AND.ICM.EQ.0)CALL VIF
      DO 620 J=1,NV
	DO 610 I=1,N
  610   XV(I)=X(LOC(I,J))
	CALL RANKER
	DO 620 I=1,N
	X(LOC(I,J))=RXV(I)
  620 CONTINUE
      CALL CORCAL
      CALL PMTRX(NV,2)
      IF(N.GT.NV.AND.ICM.EQ.0)CALL VIF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION CUMSTT(K,Z)
C   CUMULATIVE STUDENT T DISTRIBUTION (INTEGRAL OF INCOMPLETE BETA
C   FUNCTION).
C   K = # DEGREES OF FREEDOM
C   Z = ARGUMENT (-1.E34 < Z**K < 1.E34)
C   REVISION LOG:
C    4/30/93  HJI  ORIGINAL VERSION.
C    5/11/93  HJI  CHANGED TO DOUBLE PRECISION TO AVOID LOSS OF
C                  SIGNIFICANCE ON SMALL ARGUMENTS (X<-30).

      IMPLICIT NONE
      REAL*4     CUMSTT, BETAHKH
      INTEGER*4  I, K, N
      REAL*4     BETAI, Z
      REAL*8   C, DK, F, HALFPI, SUM, TERM, X

      DATA HALFPI/ 1.570796326794896/

      DK=K
      X=Z/DSQRT(DK)
      F=1./(1.+X**2)
      SUM=0.

      IF (X .LE. -40.) THEN
	 BETAI=DABS(1./X)**K/DK
	 GO TO 100
      ENDIF

      IF (MOD(K,2) .EQ. 0) THEN
	 C=1./DK
	 DO 10 I=K-2,0,-2
	    C=(I+2.)*C/(I+1.)
	    SUM=SUM*F+C
   10    CONTINUE
	 BETAI=X*SQRT(F)*SUM + C
      ELSE
	 N=(K+1)/2
	 C=1.
	 DO 50 I=1,2*N-2,2
   50    C=I*C/(I+1.)
	 IF (N .GT. 1) THEN
	    TERM=F
	    SUM=TERM
	    DO 60 I=2,2*N-3,2
	    TERM=F*TERM*I/(I+1)
   60       SUM=SUM+TERM
	 ENDIF
	 BETAI=C*(DATAN(X)+HALFPI+X*SUM)
CX         WRITE (2,'(I4,F10.4,1P2E16.6)')  K,X,DATAN(X)+HALFPI,X*SUM
      ENDIF
  100 CUMSTT=BETAI/BETAHKH(K)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DATOUT
C**************************************************************************
C SUBROUTINE DATOUT WILL OUTPUT THE SAMPLE AND ITS RANKS IF REQUESTED
C**************************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      INTEGER*4 LOC, I, J
      LOC(I,J)=(J-1)*N+I
      REWIND 4
      READ(4)X
      CALL OUTDAT(0)
      DO 620 J=1,NV
	DO 610 I=1,N
  610   XV(I)=X(LOC(I,J))
	CALL RANKER
	DO 620 I=1,N
	X(LOC(I,J))=RXV(I)
  620 CONTINUE
      CALL OUTDAT(1)
      RETURN
      END
      SUBROUTINE DATSQZ(CARD,CRDTYP,TMPCRD)
C***********************************************************************
C SUBROUTINE DATSQZ PROCESSES PARAMETER CARDS WHICH REQUIRE
C CONVERTING CHARACTER DATA TO INTEGER DATA
C***********************************************************************
      PARAMETER (LENT=11)
      INTEGER*4 IZERO, ININE, LENC, IC, IBEG, ITEST, IEND, ILEN
      INTEGER*4 IT, I
      CHARACTER CARD*(*), CRDTYP*(*), TMPCRD*(LENT), BLANK*1, MINUS*1
      PARAMETER (BLANK=' ',MINUS='-')
      IZERO=ICHAR('0')
      ININE=ICHAR('9')
      LENC=LEN(CARD)
      IC=0
C SEARCH FOR BEGINNING OF NON-BLANK CHARACTER STRING
   10 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)THEN
	WRITE(6,9001)CRDTYP
	RETURN
      ENDIF
      IF(CARD(IC:IC).EQ.BLANK)GO TO 10
      IBEG=IC
      IF(CARD(IC:IC).EQ.MINUS)GO TO 20
      ITEST=ICHAR(CARD(IC:IC))
      IF(ITEST.LT.IZERO.OR.ITEST.GT.ININE)THEN
	WRITE(6,9002)CRDTYP,CARD(IC:IC)
	RETURN
      ENDIF
C SEARCH FOR ENDING OF NON-BLANK CHARACTER STRING
   20 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)GO TO 30
      IF(CARD(IC:IC).EQ.BLANK)GO TO 30
      ITEST=ICHAR(CARD(IC:IC))
      IF(ITEST.LT.IZERO.OR.ITEST.GT.ININE)THEN
	WRITE(6,9002)CRDTYP,CARD(IC:IC)
	RETURN
      ENDIF
      GO TO 20
C MOVE NON-BLANK CHARACTER STRING INTO TMPCRD RIGHT-JUSTIFIED
   30 CONTINUE
      IEND=IC-1
      ILEN=IEND-IBEG+1
      IF(ILEN.GT.LENT)THEN
	WRITE(6,9003)CRDTYP,ILEN,LENT
	RETURN
      ENDIF
      TMPCRD=BLANK
      IT=LENT-ILEN
      DO 40 I=IBEG,IEND
	IT=IT+1
	TMPCRD(IT:IT)=CARD(I:I)
   40 CONTINUE
      RETURN
 9001 FORMAT('1',5X,'THE PARAMETER CARD ',A,'CONTAINS NO DATA')
 9002 FORMAT('1',5X,'THE PARAMETER CARD ',A,'CONTAINS THE ',
     1       'NON-NUMERIC CHARACTER ',A)
 9003 FORMAT('1',5X,'THE DATA ON PARAMETER CARD ',A,'CONTAINS ',I2,
     1       ' DIGITS',/,6X,'THE MAXIMUM NUMBER OF DIGITS ALLOWED ',
     2       'IS ',I2)
      END
      SUBROUTINE DMFSD (N,IPARM)
C***********************************************************************
C SUBROUTINE DMFSD IS USED IN INVERTING A CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM

      REAL*4 TOL, DSUM, DPIV
      INTEGER*4 N, IPARM, KPIV, K, IND, LEND, I, L, LANF
      INTEGER*4 LIND, KT
      KPIV=0
      DO 60 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
      TOL=ABS(.01*(CORR(KPIV)))
      DO 60 I=K,N
      DSUM=0.0
      IF (LEND.EQ.0) GO TO 20
      DO 10 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
   10 DSUM=DSUM+CORR(LANF)*CORR(LIND)
   20 DSUM=CORR(IND)-DSUM
      IF (I.NE.K) GO TO 50
      IF (DSUM-TOL) 30,30,40
   30 IF (DSUM.LE.0.0) GO TO 70
      KT=K-1
      WRITE(6,80)KT
   40 DPIV=SQRT(DSUM)
      CORR(KPIV)=DPIV
      DPIV=1.0/DPIV
      GO TO 60
   50 CORR(IND)=DSUM*DPIV
   60 IND=IND+I
      RETURN
   70 WRITE(6,90)K
      IPARM=-K
      RETURN
C
   80 FORMAT(20X,'ROUNDING ERROR IN ROW ',I2)
   90 FORMAT(20X,'MATRIX IS SINGULAR AT ROW ',I2)
      END
      SUBROUTINE DSINV (N)
C***********************************************************************
C SUBROUTINE DSINV IS USED IN INVERTING A CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM

      REAL*4 DIN, WORK
      INTEGER*4 N, IPARM, IPIV, IND, I, MIN, KEND, LANF, J, K
      INTEGER*4 LHOR, LVER, L
      IPARM=0
      CALL DMFSD (N,IPARM)
      IF (IPARM.LT.0) RETURN
      IPIV=N*(N+1)/2
      IND=IPIV
      DO 40 I=1,N
      DIN=1.0/CORR(IPIV)
      CORR(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF (KEND.LE.0) GO TO 30
      J=IND
      DO 20 K=1,KEND
      WORK=0.0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
      DO 10 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
   10 WORK=WORK+CORR(LVER)*CORR(LHOR)
      CORR(J)=-WORK*DIN
   20 J=J-MIN
   30 IPIV=IPIV-MIN
   40 IND=IND-1
      DO 60 I=1,N
      IPIV=IPIV+I
      J=IPIV
      DO 60 K=I,N
      WORK=0.0
      LHOR=J
      DO 50 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+CORR(LHOR)*CORR(LVER)
   50 LHOR=LHOR+L
      CORR(J)=WORK
   60 J=J+K
      RETURN
      END
      SUBROUTINE ERRCHK(MESSG)
C***********************************************************************
C SUBROUTINE ERRCHK IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2642
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO 87115
C
C     SIMPLIFIED VERSION FOR STAND-ALONE USE.     APRIL 1977
C
C     ABSTRACT
C         THE ROUTINES ERRCHK, ERXSET, AND ERRGET TOGETHER PROVIDE
C         A UNIFORM METHOD WITH SEVERAL OPTIONS FOR THE PROCESSING
C         OF DIAGNOSTICS AND WARNING MESSAGES WHICH ORIGINATE
C         IN THE MATHEMATICAL PROGRAM LIBRARY ROUTINES.
C         ERRCHK IS THE CENTRAL ROUTINE, WHICH ACTUALLY PROCESSES
C         MESSAGES.
C
C     DESCRIPTION OF ARGUMENTS
C         MESSG  - NAME OF CHARACTER VARIABLE CONTAINING THE MESSAGE,
C                  OR ELSE A CHARACTER CONSTANT CONTAINING
C                  THE MESSAGE.  BY CONVENTION, ALL MESSAGES SHOULD
C                  BEGIN WITH *IN SUBNAM, ...*, WHERE SUBNAM IS THE
C                  NAME OF THE ROUTINE CALLING ERRCHK.
C
C     EXAMPLE
C            CALL ERRCHK('IN QUAD, INVALID VALUE OF ERR.')
C
C
C
C     ERRCHK USES SUBROUTINES ERRGET, ERRPRT, ERXSET, ERSTGT
C
C***********************************************************************
      INTEGER*4 NF, NT
      CHARACTER MESSG*(*)
C
      CALL ERRGET(NF,NT)
C     IF MESSAGES ARE TO BE SUPPRESSED, RETURN
      IF (NF.EQ.0) RETURN
C     PRINT MESSAGE
      CALL ERRPRT(MESSG)
C     IF LAST MESSAGE, SAY SO
      IF (NF.EQ.1) WRITE(6,10)
   10 FORMAT (' ERRCHK MESSAGE LIMIT REACHED.')
C     DECREMENT MESSAGE COUNT
      IF (NF.GT.0) NF = NF-1
      CALL ERXSET(NF,NT)
C     IF ALL IS WELL, RETURN
      IF (NF.GE.0) RETURN
C     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL,
C     THEN EXPLAIN ERXSET USAGE.
      WRITE(6,15)
   15 FORMAT (/' *** NOTE ***'
     1/' TO MAKE THE ERROR MESSAGE PRINTED ABOVE BE NONFATAL,',
     2/' OR TO SUPPRESS THE MESSAGE COMPLETELY,',
     3/' INSERT AN APPROPRIATE CALL TO ERXSET',
     4 ' AT THE START OF YOUR PROGRAM.',
     5/' FOR EXAMPLE, TO PRINT UP TO 10 NONFATAL WARNING MESSAGES, USE',
     6/'          CALL ERXSET(10,0)'    )
      WRITE(6,20)
   20 FORMAT (/' PROGRAM ABORT DUE TO ERROR.')
      RETURN
      END
      SUBROUTINE ERRGET(NFATAL,NTRACE)
C***********************************************************************
C SUBROUTINE ERRGET IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
C     ABSTRACT
C         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.
C         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES
C         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE
C         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.
C
C     DESCRIPTION OF ARGUMENTS
C         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER.
C         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.)
C         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.)
C
C***********************************************************************
      INTEGER*4 NFATAL, NTRACE
      CALL ERSTGT(1,NFATAL,NTRACE)
      RETURN
      END
      SUBROUTINE ERRPRT(MESSG)
C***********************************************************************
C SUBROUTINE ERRPRT IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
C
C     UTILITY ROUTINE TO SIMPLY PRINT THE CHARACTER MESSAGE IN MESSG
C
C***********************************************************************
      INTEGER*4 FIRST, LAST
      CHARACTER MESSG*(*)
    1 FORMAT (1X,A)
      FIRST = 1
      DO WHILE (FIRST .LE. LEN(MESSG))
        LAST = FIRST + 127
        IF (LAST .GT. LEN(MESSG)) LAST = LEN(MESSG)
        WRITE(6,1) MESSG(FIRST:LAST)
        FIRST = LAST + 1
      END DO
      RETURN
      END
      SUBROUTINE ERSTGT(K,NFATAL,NTRACE)
C***********************************************************************
C SUBROUTINE ERSTGT IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
C
C     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS
C     THE FLAGS AS LOCAL VARIABLES.
C
C     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN
C     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE
C     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING
C     VALUES IN THE MAIN PROGRAM.
C
C***********************************************************************
      INTEGER*4 K, NFATAL, NTRACE, LNF, LNT
      DATA LNF/-1/,LNT/0/
      IF (K.LE.0) LNF = NFATAL
      IF (K.LE.0) LNT = NTRACE
      IF (K.GT.0) NFATAL = LNF
      IF (K.GT.0) NTRACE = LNT
      RETURN
      END
      SUBROUTINE ERXSET(NFATAL,NTRACE)
C***********************************************************************
C SUBROUTINE ERXSET IS AN ERROR CHECKING ROUTINE USED IN
C GENERATING A BETA DISTRIBUTION
C
C
C     ABSTRACT
C         ERXSET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK.
C         ERXSET ASSIGNS THE VALUES OF NFATAL AND NTRACE RESPECTIVELY
C         TO NF AND NT IN COMMON BLOCK MLBLK0 THEREBY SPECIFYING THE
C         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK.
C
C     DESCRIPTION OF ARGUMENTS
C         BOTH ARGUMENTS ARE INPUT ARGUMENTS OF DATA TYPE INTEGER.
C         NFATAL - IS A FATAL-ERROR / MESSAGE-LIMIT FLAG. A NEGATIVE
C                  VALUE DENOTES THAT DETECTED DIFFICULTIES ARE TO BE
C                  TREATED AS FATAL ERRORS.  NONNEGATIVE MEANS NONFATAL.
C                  A NONNEGATIVE VALUE IS THE MAXIMUM NUMBER OF NONFATAL
C                  WARNING MESSAGES WHICH WILL BE PRINTED BY ERRCHK,
C                  AFTER WHICH NONFATAL MESSAGES WILL NOT BE PRINTED.
C                  (DEFAULT VALUE IS -1.)
C         NTRACE - .GE.1 WILL CAUSE A TRACE-BACK TO BE GIVEN,
C                        IF THIS FEATURE IS IMPLEMENTED ON THIS SYSTEM.
C                  .LE.0 WILL SUPPRESS ANY TRACE-BACK, EXCEPT FOR
C                        CASES WHEN EXECUTION IS TERMINATED.
C                  (DEFAULT VALUE IS 0.)
C
C         *NOTE* -- SOME CALLS TO ERRCHK WILL CAUSE UNCONDITIONAL
C         TERMINATION OF EXECUTION.  ERXSET HAS NO EFFECT ON SUCH CALLS.
C
C     EXAMPLES
C         1. TO PRINT UP TO 100 MESSAGES AS NONFATAL WARNINGS USE
C            CALL ERXSET(100,0)
C         2. TO SUPPRESS ALL MATHLIB WARNING MESSAGES USE
C            CALL ERXSET(0,0)
C
C     ERXSET USES SUBROUTINES ERSTGT
C
C***********************************************************************
      INTEGER*4 NFATAL, NTRACE, I0
      I0 = 0
      CALL ERSTGT(I0,NFATAL,NTRACE)
      RETURN
      END
      SUBROUTINE EXPONT(J)
C***********************************************************************
C SUBROUTINE EXPONT GENERATES AN EXPONENTIAL DISTRIBUTION
C
C         f(X) = C EXP(C*X)/(EXP(C*B)-EXP(C*A)),             A <= X <= B
C              = C EXP(C*(X-A)/(EXP(C*(B-A) - 1)
C
C         F(X) = (EXP(C*X) - EXP(C*A))/(EXP(C*B)-EXP(C*A)),  A <= X <= B
C              = (EXP(C*(X-A)) - 1)/(EXP(C*(B-A)) - 1)
C
C         C CAN BE EITHER POSITIVE OR NEGATIVE.
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 LAMBDA, PROBINC, A, B, EXBA, STRTPT, R
      REAL*4 RNUM
      INTEGER*4 J, LOC, I
      LOC(I,J)=(J-1)*N+I
      PROBINC=1./FLOAT(N)
      IF(IRS.EQ.1)PROBINC=1.0
      READ(8)A,B,LAMBDA
C  THE FOLLOWING CAUSES NORMALIZATION (F(T)=1) AT T = B.
      EXBA=EXP(LAMBDA*(B-A))-1.
      STRTPT=0.
      DO 10 I=1,N
        RNUM=RAN(ISEED)
	R=PROBINC*RNUM+STRTPT
	X(LOC(I,J))=A+ALOG(1.+R*EXBA)/LAMBDA
	IF(IRS.EQ.0)STRTPT=STRTPT+PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   10 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE FINDIT(NP,M,EIG,ICONV)
C***********************************************************************
C SUBROUTINE FINDIT IS USED IN THE POSITIVE DEFINITE CHECK
C OF THE CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/PDMAT/Z, D
      REAL*4       Z(NVAR,NVAR), D(NVAR)
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 EIG
      INTEGER*4 NP, M, ICONV, LOC, I, J, NEV, L1, L2, K, KI
      LOC(I,J)=(J-1)*N+I
      NEV = 0
      DO 9 I = 1,NP
    9 IF(D(I).LT.0.0)NEV=NEV+1
      IF(NEV.EQ.0)THEN
	ICONV=1
	RETURN
      ELSE
	DO 11 I = 1,NEV
   11   D(I)=EIG
	L1=NEV+1
	L2=NEV+NEV
	DO 12 I = L1,L2
	  IF(D(I).LT.EIG)D(I)=EIG
   12   CONTINUE
	DO 4 I = 1,M
	DO 4 J = 1,M
    4   X(LOC(I,J)) = 0.0
	DO 5 I = 1,NP
	DO 5 J = 1,NP
	DO 5 K = 1,NP
    5   X(LOC(I,J)) = X(LOC(I,J)) + Z(I,K)*D(K)*Z(J,K)
	DO 30 I = 1,NP
   30   X(LOC(I,I)) = 1.0
	KI = 0
	DO 10 I = 1,NP
	DO 10 J = 1,I
	  KI = KI + 1
	  CORR(KI) = X(LOC(I,J))
   10   CONTINUE
      ENDIF
      RETURN
      END
      FUNCTION FINVNOR (X)
C***********************************************************************
C SUBROUTINE FINVNOR IS USED IN GENERATING THE NORMAL AND
C LOGNORMAL DISTRIBUTIONS
C***********************************************************************
      REAL*4 FINVNOR, X, F, Y, RIERFC1
      IF (X-0.5) 20,10,30
   10 FINVNOR=0.0
      RETURN
   20 F=-1.0
      Y=X
      GO TO 40
   30 F=1.0
      Y=1.0-X
   40 FINVNOR=SQRT(2.0)*F*RIERFC1(2.*Y)
      RETURN
      END
      SUBROUTINE HISTO
C***********************************************************************
C SUBROUTINE HISTO GENERATES HISTOGRAMS OF THE VARIABLES
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      PARAMETER (IFLG=20)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)

      REAL*4 SUM, SUMSQ, XMEAN, XVAR, FMED, R, SIZE, POWER, CON, TEMP
      REAL*4 CELL, POINT
      INTEGER*4 I, N2, LOC, M, KSUM, K, KK, J
      CHARACTER KX*1
C     THIS IS A MODIFICATION OF DISTAT -- IT MAKES THE HISTOGRAM
C     AND PRINTS THE MEAN AND THE MEDIAN
C      N = NO. OF OBSERVATIONS
C      XV = ARRAY NAME
      DATA KX/ 'X' /
      IF (N-1) 10,20,30
   10 WRITE(6,260)
      RETURN
   20 WRITE(6,270)XV(1)
      RETURN
   30 CONTINUE
      CALL SIFT (XV,N)
      SUM=0.0
      SUMSQ=0.0
      DO 40 I=1,N
      SUM=SUM+XV(I)
      SUMSQ=SUMSQ+XV(I)*XV(I)
   40 CONTINUE
      XMEAN=SUM/FLOAT(N)
      XVAR=(SUMSQ-FLOAT(N)*XMEAN*XMEAN)/FLOAT(N)
      IF(IRS.NE.0)XVAR=XVAR*FLOAT(N)/FLOAT(N-1)
      N2=N/2
      LOC=(N+1)/2-N2
      FMED=(XV(N2+LOC)+XV((N+1-LOC)/2+1))/2.0
      R=XV(N)-XV(1)
      IF (R.NE.0.0) GO TO 50
      WRITE(6,280)
      GO TO 250
   50 CONTINUE
      M=1
      SIZE=R/IFLG
      POWER=LOG10(SIZE)
      IF (POWER.GE.0.) GO TO 60
      POWER=AINT(POWER)
      GO TO 70
   60 POWER=AINT(POWER)+1.
   70 SIZE=SIZE/10.**POWER
      CON=.01
   80 IF (SIZE.LE.(CON+.005)) GO TO 90
      CON=CON+.01
      GO TO 80
   90 SIZE=CON*10.**POWER
      TEMP=XV(1)/SIZE
      IF (TEMP) 100,100,110
  100 TEMP=AINT(TEMP-1.0)
      GO TO 120
  110 TEMP=AINT(TEMP)
  120 IF (XV(1)-TEMP*SIZE.GT.0.0) GO TO 130
      TEMP=TEMP-0.5
  130 CELL=TEMP*SIZE+SIZE
      POINT=TEMP*SIZE+SIZE*0.5
      I=0
      KSUM=0
      WRITE(6,290)
      K=-1
      M=1
  140 I=I+1
      K=K+1
  150 IF (I-N) 160,160,190
  160 IF (XV(I)-CELL) 140,140,170
  170 IF (K) 180,180,200
  180 WRITE(6,310)POINT,K
      GO TO 240
  190 M=2
  200 IF (K-90) 210,210,220
  210 KK=K
      GO TO 230
  220 KK=90
  230 WRITE(6,310)POINT,K,(KX,J=1,KK)
      KSUM=KSUM+K
      IF (M.GT.1) GO TO 250
      K=0
  240 CELL=CELL+SIZE
      POINT=POINT+SIZE
      GO TO 150
  250 WRITE(6,320)KSUM
      WRITE(6,300)XV(1),XV(N),R,XMEAN,FMED,XVAR
      RETURN
C
  260 FORMAT(' N IS ZERO',//)
  270 FORMAT(' ONE OBS. ',E17.8,//)
  280 FORMAT(' NO HISTOGRAM - RANGE =0',/)
  290 FORMAT(/,5X,'MIDPOINT',10X,'FREQ.',/)
  300 FORMAT(//,6X,'MIN',12X,'MAX',11X,'RANGE',11X,'MEAN',10X,
     1       'MEDIAN',8X,'VARIANCE',//,1X,6G15.7,/)
  310 FORMAT(1X,G15.7,5X,I5,4X,90A1)
  320 FORMAT('0',20X,I5)
      END
      SUBROUTINE HPSRT
C***********************************************************************
C SUBROUTINE HPSRT IS USED IN THE RANKING OF THE SAMPLED DATA
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      REAL*4 XHOLD, YHOLD
      INTEGER*4 R, L, J, I
      L=N/2+1
      R=N
   10 IF (L.LE.1) GO TO 70
      L=L-1
      XHOLD=XV(L)
      YHOLD=RXV(L)
   20 J=L
   30 I=J
      J=2*J
      IF (J-R) 40,50,60
   40 IF (XV(J).LT.XV(J+1)) J=J+1
   50 IF (XHOLD.GE.XV(J)) GO TO 60
      XV(I)=XV(J)
      RXV(I)=RXV(J)
      GO TO 30
   60 XV(I)=XHOLD
      RXV(I)=YHOLD
      GO TO 10
   70 XHOLD=XV(R)
      YHOLD=RXV(R)
      XV(R)=XV(1)
      RXV(R)=RXV(1)
      R=R-1
      IF (R.GT.1) GO TO 20
      XV(1)=XHOLD
      RXV(1)=YHOLD
      RETURN
      END
      SUBROUTINE HSTOUT
C***********************************************************************
C SUBROUTINE HSTOUT IS USED TO GENERATE HISTOGRAMS OF THE
C VARIABLES
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      PARAMETER (LEND=14)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      INTEGER*4 LOC, I, J, IDT
      CHARACTER*15 DIST(LEND)
      DATA DIST(1)/'BETA'/,DIST(2)/'NORMAL'/,DIST(3)/'LOGNORMAL'/,
     1     DIST(4)/'UNIFORM'/,DIST(5)/'LOGUNIFORM'/,
     2     DIST(6)/'UNIFORM*'/,DIST(7)/'LOGUNIFORM*'/,
     3     DIST(8)/'TRIANGULAR'/,DIST(9)/'USER SUPPLIED'/,
     4     DIST(10)/'EXPONENTIAL'/, DIST(11)/'RAYLEIGH '/,
     5     DIST(12)/'RAYL-EXPON '/, DIST(13)/'STUDENT '/,
     6     DIST(14)/'LOGSTUDENT '/
      LOC(I,J)=(J-1)*N+I
      REWIND 4
      READ(4)X
      DO 590 I=1,NV
	IDT=IDIST(I)
	WRITE(6,9001) TITLE
	WRITE(6,9002)I,DIST(IDT)
	DO 530 J=1,N
  530   XV(J)=X(LOC(J,I))
	CALL HISTO
  590 CONTINUE
 9001 FORMAT('1   ',A100)
 9002 FORMAT('0  HISTOGRAM FOR VARIABLE NO.',I3,5X,A,' DISTRIBUTION')
      RETURN
      END
      FUNCTION HYPGEO(A,B,C,X,SUM)
C***********************************************************************
C FUNCTION HYPGEO IS USED IN GENERATING A BETA DISTRIBUTION
C
C     FUNCTION HYPGEO SUMS THE SERIES FOR F(A,B,C,X) WHEN SUM=1. OR
C     F(A,B,C,X)-1. WHEN SUM=0. WHERE F(A,B,C,X) IS THE GAUSS
C     HYPERGEOMETRIC FUNCTION.
C     UPDATED, JANUARY, 1975
C
C***********************************************************************
      REAL*4 HYPGEO, A, B, C, X, SUM, TOL, R1MACH, FSUM, TK, TLIM, AK
      INTEGER*4 KMAX, M, K     
      KMAX=1000
      TOL=MAX(1.0E-14,R1MACH(4))
      M=2.-SUM
      FSUM=X*A*(B/C)
      TK=FSUM
      TLIM=TOL
      IF(M.EQ.2) TLIM=TLIM*ABS(FSUM)
      AK=1.
      DO 10 K=2,KMAX
      TK=TK*((A+AK)/(AK+1.))*((B+AK)/(C+AK))*X
      FSUM=FSUM+TK
      IF(ABS(TK).LE.TLIM) GO TO 20
      AK=AK+1.
   10 CONTINUE
      CALL ERRCHK('IN HYPGEO, NO CONVERGENCE OF SERIES')
C      HYPGEO=0.
      RETURN
C      RETURN
   20 HYPGEO=FSUM+SUM
      RETURN
      END
      FUNCTION I1MACH(I)
C***********************************************************************
C FUNCTION I1MACH SETS INTEGER MACHINE DEPENDENT CONSTANTS
C
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  801001   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   These machine constant routines must be activated for
C   a particular environment.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
C***********************************************************************
      INTEGER*4 I1MACH, IMACH(16), OUTPUT, I
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C      DATA IMACH/       7,    2,    2,    2,   36,    4,    2,   33,
C     +         Z1FFFFFFFF,    2,   24, -256,  255,   60, -256,  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C      DATA IMACH/       5,    6,    7,    6,   48,    6,    2,   39,
C     +  O0007777777777777,    8,   13,  -50,   76,   26,  -50,   76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C      DATA IMACH/       5,    6,    7,    6,   48,    6,    2,   39,
C     +  O0007777777777777,    8,   13,  -50,   76,   26,-32754,32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C      DATA IMACH/       5,    6,    7,6LOUTPUT,60,   10,    2,   48,
C     +  00007777777777777777B,2,   47, -929, 1070,   94, -929, 1069 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C      DATA IMACH/       5,    6,  102,    6,   64,    8,    2,   63,
C     + 777777777777777777777B,2,   47,-8189, 8190,   94,-8099, 8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C      DATA IMACH/      11,   12,    8,   10,   16,    2,    2,   15,
C     +              32767,   16,    6,  -64,   63,   14,  -64,   63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C      DATA IMACH/       5,    6,    0,    6,   24,    3,    2,   23,
C     +            8388607,    2,   23, -127,  127,   38, -127,  127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C      DATA IMACH/       5,    6,   43,    6,   36,    6,    2,   35,
C     +      O377777777777,    2,   27, -127,  127,   63, -127,  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C      DATA IMACH/       5,    6,    4,    1,   16,    2,    2,   15,
C     +              32767,    2,   23, -128,  127,   39, -128,  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C      DATA IMACH/       5,    6,    4,    1,   16,    2,    2,   15,
C     +              32767,    2,   23, -128,  127,   55, -128,  127 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C      DATA IMACH/       5,    6,    7,    6,   32,    4,   16,   31,
C     +          Z7FFFFFFF,   16,    6,  -64,   63,   14,  -64,   63 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C      DATA IMACH/       5,    6,    5,    6,   36,    5,    2,   35,
C     +      "377777777777,    2,   27, -128,  127,   54, -101,  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C      DATA IMACH/       5,    6,    5,    6,   36,    5,    2,   35,
C     +      "377777777777,    2,   27, -128,  127,   62, -128,  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C      DATA IMACH/       5,    6,    5,    6,   32,    4,    2,   31,
C     +         2147483647,    2,   24, -127,  127,   56, -127,  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C      DATA IMACH/       5,    6,    5,    6,   16,    2,    2,   15,
C     +              32767,    2,   24, -127,  127,   56, -127,  127 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C      DATA IMACH/       5,    6,    1,    6,   36,    4,    2,   35,
C     +      O377777777777,    2,   27, -128,  127,   60,-1024, 1023 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FOR COMPILER
C      DATA IMACH/       5,    6,    7,    6,   36,    6,    2,   35,
C     +      O377777777777,    2,   27, -128,  127,   60,-1024, 1023 /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
      DATA IMACH/        5,    6,    5,    6,   32,    4,    2,   31,
     +          2147483647,    2,   24, -127,  127,   56, -127,  127 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780. G-FLOATING OPTION
C      DATA IMACH/       5,    6,    5,    6,   32,    4,    2,   31,
C     +         2147483647,    2,   24, -127,  127,   54, -511,  511 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
C
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH=IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE(OUTPUT,9000)
9000  FORMAT('1ERROR    1 IN I1MACH - I OUT OF BOUNDS' )
C
      RETURN
      END
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
C***********************************************************************
C SUBROUTINE IMTQL2 IS USED IN THE POSTIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  IMTQL2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EISPACK,EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH ET AL
C***PURPOSE  COMPUTES EIGENVALUES AND EIGENVECTORS OF SYMMETRIC
C            TRIDIAGONAL MATRIX USING IMPLICIT QL METHOD.
C***DESCRIPTION
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG(A,B) FOR SQRT(A**2 + B**2).
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  *MATRIX EIGENSYSTEM ROUTINES-EISPACKGUIDE*,
C                 B.T.SMITH,J.M.BOYLE,J.J.DONGARRA,B.S.GARBOW,
C                 Y.I.KEBE,V.C.KLEMA,C.B.MOLER,SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  IMTQL2
C
C***********************************************************************
      INTEGER*4 I,J,K,L,M,N,II,NM,MML,IERR
      REAL*4 D(N),E(N),Z(NM,N)
      REAL*4 B,C,F,G,P,R,S,S1,S2
      REAL*4 PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  IMTQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 240 L = 1, N
	 J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
	    IF (M .EQ. N) GO TO 120
	    S1 = ABS(D(M)) + ABS(D(M+1))
	    S2 = S1 + ABS(E(M))
	    IF (S2 .EQ. S1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
	 IF (M .EQ. L) GO TO 240
	 IF (J .EQ. 30) GO TO 1000
	 J = J + 1
C     .......... FORM SHIFT ..........
	 G = (D(L+1) - P) / (2.0E0 * E(L))
	 R = PYTHAG(G,1.0E0)
	 G = D(M) - P + E(L) / (G + SIGN(R,G))
	 S = 1.0E0
	 C = 1.0E0
	 P = 0.0E0
	 MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
	 DO 200 II = 1, MML
	    I = M - II
	    F = S * E(I)
	    B = C * E(I)
	    IF (ABS(F) .LT. ABS(G)) GO TO 150
	    C = G / F
	    R = SQRT(C*C+1.0E0)
	    E(I+1) = F * R
	    S = 1.0E0 / R
	    C = C * S
	    GO TO 160
  150       S = F / G
	    R = SQRT(S*S+1.0E0)
	    E(I+1) = G * R
	    C = 1.0E0 / R
	    S = S * C
  160       G = D(I+1) - P
	    R = (D(I) - G) * S + 2.0E0 * C * B
	    P = S * R
	    D(I+1) = G + P
	    G = C * R - B
C     .......... FORM VECTOR ..........
	    DO 180 K = 1, N
	       F = Z(K,I+1)
	       Z(K,I+1) = S * Z(K,I) + C * F
	       Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
C
  200    CONTINUE
C
	 D(L) = D(L) - P
	 E(L) = G
	 E(M) = 0.0E0
	 GO TO 105
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
	 I = II - 1
	 K = I
	 P = D(I)
C
	 DO 260 J = II, N
	    IF (D(J) .GE. P) GO TO 260
	    K = J
	    P = D(J)
  260    CONTINUE
C
	 IF (K .EQ. I) GO TO 300
	 D(K) = D(I)
	 D(I) = P
C
	 DO 280 J = 1, N
	    P = Z(J,I)
	    Z(J,I) = Z(J,K)
	    Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C-----------------------------------------------------------------------  
      SUBROUTINE INVBETA(BVALUE,X,P,Q,ABSERR,RELERR,FAILURE)
C
C   INVERT THE INCOMPLETE BETA FUNCTION: FIND X SUCH THAT 
C     BETAIC(X,1-X,P,Q,1,ANSWER,UNDERFLOW) RETURNS ANSWER=BVALUE.
C     NOTE: BETAIC REQUIRES [0 <= X <= 1].
C
C   REVISION LOG:
C    2/14/97  HJI  INITIAL VERSION.
C
C   ALGORITHM:  SINCE THE INCOMPLETE BETA FUNCTION INCREASES 
C   MONOTONICALLY FROM 0 T0 1, X CAN BE FOUND BY BISECTION. 
C
      IMPLICIT NONE
      INTEGER*4  I, MAXITER, NITER, NRET, UNDERFLOW 
      REAL*4     ABSERR, BVALUE, FI, FF, FMID(1), OMXMID, P, Q, RELERR,
     &           X,XI, XF, XMID 
      LOGICAL  FAILURE

C   150 ITERATIONS SHOULD ALWAYS GIVE AN ANSWER UNLESS THE BETAIC 
C   FUNCTION IS CORRUPT.
      DATA MAXITER/150/ 

      FAILURE=.FALSE.  

C   SPECIAL CASES:
      IF (BVALUE .EQ. 0.) THEN
         X=0.
         RETURN
      ELSE IF (BVALUE .EQ. 1) THEN
         X=1.
         RETURN
      ENDIF 

C   CHECK FOR ILLEGAL VALUE OF BVALUE, P, OR Q:
      IF (BVALUE .LT. 0. .OR. BVALUE .GT. 1. .OR.
     &    P .LE. 0. .OR. Q .LE. 0.) THEN
         FAILURE=.TRUE.
         RETURN 
      ENDIF 

      XI=0.
      FI=0.
      XF=1.
      FF=1.
      NRET=1
      NITER=0

   20 XMID=0.5*(XI+XF)
      OMXMID=1.-XMID 
      CALL BETAIC(XMID,OMXMID,P,Q,NRET,FMID(1),UNDERFLOW) 
      X=XMID 
      
      IF (ABS(FMID(1)-BVALUE) .LE. ABSERR  .OR.
     &    ABS(FMID(1)-BVALUE) .LE. RELERR*BVALUE) RETURN

C   IF THE FUNCTION IS TOO FLAT TO ALLOW THE DESIRED ACCURACY,
C   USE THE BEST POSSIBLE VALUE.
      IF (XMID .EQ. XI) THEN 
         IF (ABS(FMID(1)-BVALUE) .GT. ABS(FF-BVALUE)) X=XF
         RETURN
      ELSE IF (XMID .EQ. XF) THEN 
         IF (ABS(FMID(1)-BVALUE) .GT. ABS(FI-BVALUE)) X=XI
         RETURN
      ENDIF

      NITER=NITER+1
      IF (FMID(1) .GT. BVALUE) THEN
         FF=FMID(1)
         XF=XMID
      ELSE 
         FI=FMID(1)
         XI=XMID
      ENDIF
      IF (NITER .LT. MAXITER) GO TO 20 
      PRINT '(A,I4,A)', ' INVBETA FAILURE: > ',MAXITER,' ITERATIONS'
      PRINT '(A,1PE10.2,A,E10.2/)', 
     &  ' ABSERR=',ABS(FMID(1)-BVALUE),
     &  '  RELERR=',ABS(FMID(1)-BVALUE)/BVALUE
      FAILURE=.TRUE.
      RETURN
      END                                                                   
C======================================================================
      REAL*4 FUNCTION INVCST(K,CUMP)
C   INVERST CUMULATIVE STUDENT T DISTRIBUTION
C   REVISION LOG:
C    5/11/93  HJI  ORIGINAL VERSION.  (ITERATIVE SEARCH ALGORITHM).
      IMPLICIT NONE
      REAL*4   CUMSTT
      INTEGER*4  J, K
      REAL*4     CUMP, F, FX, F1, F2, U, X, X1, X2
      REAL*4   PI

      DATA PI    / 3.141592653589793/

      IF (K .EQ. 1) THEN
	 X=TAN(PI*(CUMP-0.5))
	 GO TO 100
      ELSE IF (K .EQ. 2) THEN
	 X=(CUMP-0.5)*SQRT(2./(CUMP*(1.-CUMP)))
	 GO TO 100
      ENDIF

      X=0.
      IF (CUMP .EQ. 0.5) GO TO 100

C   ESTIMATE X (INVERSE CUMULATIVE PROBABILITY) AND UPPER AND LOWER
C   BOUNDS X2 AND X1.
      U=(CUMP-0.5)*SQRT(2./(CUMP*(1.-CUMP)))
      X=SIGN(ABS(U)**((K+1.80)/(3*K-3.)),U)
      IF (U .GT. 0.) THEN
	 X2=1.25*X+1.
	 X1=AMAX1(0.,SQRT(X)-1.)
      ELSE
	 X1=1.25*X-1.
	 X2=AMIN1(0.,-SQRT(ABS(X))+1.)
      ENDIF

      F=CUMP
      F1=CUMSTT(K,X1)
      F2=CUMSTT(K,X2)

C   Iteratively improve the bounds on X by linear interpolation for
C   one bound and then by halving the interval (usually) for the other
C   bound.
      J=0
   10 CONTINUE
      J=J+1
      IF (F2 .LT. F1) THEN
	 PRINT '(A,I4,A,I4,A,F10.2,A,F10.2)',
     &         ' F2 < F1 FOR K=',K,', J=',J,', X1=',X1,', X2=',X2
      ENDIF
      IF (F1 .GT. F .OR. F .GT. F2) THEN
	 PRINT '(A,I4,A,I4,A,F10.2,A,F10.2)',
     &   ' VALUE NOT BOUNDED FOR K=',K,', J=',J,', X1=',X1,', X2=',X2
      ENDIF
      IF (ABS(F1-F2) .LT. 5.E-10) GO TO 100
      X=X1+(X2-X1)*(F-F1)/(F2-F1)
      FX=CUMSTT(K,X)
      IF (ABS(F-FX) .LT. 1.E-6*ABS(F+FX)) GO TO 100
      IF (FX .LT. F) THEN
	 X1=X
	 F1=FX
	 X=0.5*(X1+X2)
	 FX=CUMSTT(K,X)
	 IF (FX .GT. F .AND. FX .LT. F2) THEN
	    X2=X
	    F2=FX
	 ELSE IF (FX .LT. F .AND. FX .GT. F1) THEN
	    X1=X
	    F1=FX
	 ENDIF
      ELSE
	 X2=X
	 F2=FX
	 X=0.5*(X1+X2)
	 FX=CUMSTT(K,X)
	 IF (FX .LT. F .AND. FX .GT. F1) THEN
	    X1=X
	    F1=FX
	 ELSE IF (FX .GT. F .AND. FX .LT. F2) THEN
	    X2=X
	    F2=FX
	 ENDIF
      ENDIF
      IF (X2-X1 .GT. 1.E-5*ABS(X)) GO TO 10

  100 INVCST=X
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***********************************************************************
C FUNCTION J4SAVE IS USED IN THE POSTIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  J4SAVE
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  SAVES AND RECALLS SOME GLOBAL VARIABLES NEEDED BY LIBRARY
C            ERROR HANDLING ROUTINES
C***DESCRIPTION
C     ABSTRACT
C        J4SAVE SAVES AND RECALLS SEVERAL GLOBAL VARIABLES NEEDED
C        BY THE LIBRARY ERROR HANDLING ROUTINES.
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        IWHICH - INDEX OF ITEM DESIRED.
C                 = 1 REFERS TO CURRENT ERROR NUMBER.
C                 = 2 REFERS TO CURRENT ERROR CONTROL FLAG.
C                 = 3 REFERS TO CURRENT UNIT NUMBER TO WHICH ERROR
C                     MESSAGES ARE TO BE SENT.  (0 MEANS USE STANDARD.)
C                 = 4 REFERS TO THE MAXIMUM NUMBER OF TIMES ANY
C                     MESSAGE IS TO BE PRINTED (AS SET BY XERMAX).
C                 = 5 REFERS TO THE TOTAL NUMBER OF UNITS TO WHICH
C                     EACH ERROR MESSAGE IS TO BE WRITTEN.
C                 = 6 REFERS TO THE 2ND UNIT FOR ERROR MESSAGES
C                 = 7 REFERS TO THE 3RD UNIT FOR ERROR MESSAGES
C                 = 8 REFERS TO THE 4TH UNIT FOR ERROR MESSAGES
C                 = 9 REFERS TO THE 5TH UNIT FOR ERROR MESSAGES
C        IVALUE - THE VALUE TO BE SET FOR THE IWHICH-TH PARAMETER,
C                 IF ISET IS .TRUE. .
C        ISET   - IF ISET=.TRUE., THE IWHICH-TH PARAMETER WILL BE
C                 GIVEN THE VALUE, IVALUE.  IF ISET=.FALSE., THE
C                 IWHICH-TH PARAMETER WILL BE UNCHANGED, AND IVALUE
C                 IS A DUMMY PARAMETER.
C      --OUTPUT--
C        THE (OLD) VALUE OF THE IWHICH-TH PARAMETER WILL BE RETURNED
C        IN THE FUNCTION VALUE, J4SAVE.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     ADAPTED FROM BELL LABORATORIES PORT LIBRARY ERROR HANDLER
C     LATEST REVISION ---  23 MAY 1979
C***REFERENCES  JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR
C                 HANDLING PACKAGE*, SAND78-1189, SANDIA LABORATORIES,
C                 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
C
C***********************************************************************
      LOGICAL ISET
      INTEGER*4 J4SAVE, IWHICH, IVALUE, IPARAM(9)
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,1,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION LENTEXT(TEXT)
C   RETURN THE LENGTH OF A STRING TO THE FIRST NULL OR ELSE TO THE
C   LAST NON-BLANK CHARACTER.
C
C   REVISION LOG:
C   10/24/90  HJI  ORIGINAL VERSION.

      IMPLICIT  NONE
      CHARACTER*(*) CHAR,TEXT
      INTEGER*4  LENTEXT, I, K, L, LEN

C   LOOK FOR A NULL IN THE STRING
      L=INDEX(TEXT,CHAR(0))
      IF (L .GT. 0) THEN
	 LENTEXT=L-1
      ELSE
	 L=LEN(TEXT)
C   FIND THE LAST NON-BLANK.  RETURN 1 IF THE STRING IS BLANK.
	 DO 10 I=1,L
	 K=L+1-I
	 IF (TEXT(K:K) .NE. ' ') GO TO 20
   10    CONTINUE
   20    LENTEXT=K
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATINV
C***********************************************************************
C SUBROUTINE MATINV IS USED TO INVERT A LOWER TRIANGULAR MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/WORKC/Q, S
      REAL*4       Q((NVAR*(NVAR+1))/2),S((NVAR*(NVAR+1))/2)

      REAL*4 TEMP
      INTEGER*4 LOC1, I, J, K, JPLUS
      LOC1(I,J)=J+(I*I-I)/2
      DO 10 I=1,NV
   10 Q(LOC1(I,I))=1.0/Q(LOC1(I,I))
      K=NV
   20 J=K-1
   30 TEMP=0.0
      JPLUS=J+1
      DO 40 I=JPLUS,K
   40 TEMP=TEMP+Q(LOC1(K,I))*Q(LOC1(I,J))
      Q(LOC1(K,J))=-TEMP*Q(LOC1(J,J))
      J=J-1
      IF (J.GE.1) GO TO 30
      K=K-1
      IF (K.GT.1) GO TO 20
      RETURN
      END
      SUBROUTINE MIX
C***********************************************************************
C SUBROUTINE MIX IMPLEMENTS THE SCHEME DESCRIBED IN IMAN AND CONOVER
C (1982) -SEE USERS GUIDE- FOR INDUCING A DESIRED CORRELATION STRUCTURE
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      COMMON/WORKX/XX
      REAL*4       XX(NMAX*NVAR)
      COMMON/WORKC/Q, S
      REAL*4       Q((NVAR*(NVAR+1))/2),S((NVAR*(NVAR+1))/2)

      REAL*4 CMXOLD, DN, FINVNOR, CMX, ELEM
      INTEGER*4 LOC1, I, J, LOC, NVN, IJK, IR, JJ, KK, NDX, IMIN
      INTEGER*4 NVX, K, KX, IT, ITX, IW
      LOC1(I,J)=J+(I*I-I)/2
      LOC(I,J)=(J-1)*N+I
C
C     IF ICM EQUALS 0 THEN A MIX OF THE GENERATED DATA IS
C     OBTAINED.  OTHERWISE, A MIX OF RANDOM INTEGERS IS
C     GENERATED AND MULTIPLIED BY THE CORRELATION MATRIX.
C     THE RESULT IS RANKED AND THE RANKS ARE USED TO GENERATE
C     THE RANDOM SAMPLE.
C
C     FOR THE FULL RANK CASE(N.GT.NV)THE FOLLOWING PROCEDURE IS USED,
C     GENERATE A SAMPLE SUCH THAT NONE OF THE OFF-DIAGONAL ELEMENTS OF
C     THE CORRESPONDING CORRELATION MARTIX ARE EQUAL TO ONE IN ABSOLUTE
C     VALUE.CHOLESKY FACTORIZATION IS DONE NEXT.
C     FOR THE LESS THAN FULL RANK CASE(N.LE.NV) 25 SAMPLES ARE GENERATED
C     AND THAT SAMPLE WHOSE CORRESPONDING CORRELATION MATRIX HAS THE
C     SMALLEST MAXIMUM(IN ABSOLUTE VALUE)OFF-DIAGONAL ELEMENT IS THE
C     SAMPLE THAT IS CHOSEN. CHOLESKY FACTORIZATION IS SKIPPED IN
C     THIS CASE.
C
      NVN=NV*N
      CMXOLD=9999.0
      DN=N+1
      DO 1 I = 1,N
      XV(I) = FINVNOR(FLOAT(I)/DN)
    1 CONTINUE
      DO 50 IJK=1,25
   10 CONTINUE
      DO 2 I=1,N
    2 IWK(I)=I
      DO 30 I=1,NV
      DO 20 J=2,N
      IR=N+2-J
      JJ=IR*RAN(ISEED)+1
      KK=IWK(IR)
      IWK(IR)=IWK(JJ)
   20 IWK(JJ)=KK
      DO 30 J=1,N
      NDX=IWK(J)
      X(LOC(J,I))=XV(NDX)
   30 CONTINUE
      IF(NV.EQ.1.OR.IRP.EQ.1)GO TO 235
C
C     COMPUTE THE CORRELATION MATRIX ON THE VDW SCORES
C
      CALL CORCAL
      IF(N.LE.NV)GO TO 37
      DO 33 I=2,NV
      IMIN=I-1
      DO 33 J=1,IMIN
   33 IF(ABS(CORR(LOC1(I,J))).GE.0.9999)GO TO 10
      GO TO 60
   37 CONTINUE
      CMX=0.0
      DO 40 I=2,NV
      IMIN=I-1
      DO 40 J=1,IMIN
      ELEM=ABS(CORR(LOC1(I,J)))
      IF(ELEM.GT.CMX)CMX=ELEM
   40 CONTINUE
      IF(CMX.GE.CMXOLD)GO TO 50
      CMXOLD=CMX
      REWIND 2
      WRITE(2)(X(I),I=1,NVN)
   50 CONTINUE
      REWIND 2
      READ(2)(X(I),I=1,NVN)
C
C     X IS AN N X NV MATRIX OF VDW SCORES. IF THIS IS A FULL RANK CASE
C     (N.GT.NV) THE RANK CORRELATION MATRIX OF X MUST BE FOUND AND
C     THE CHOLESKY FACTORIZATION PERFORMED. IF THIS IS NOT A FULL RANK
C     CASE (N.LE.NV) THIS MATRIX OF VDW SCORES WILL BE USED TO MIX THE
C     INTERVALS OF THE RAW DATA MATRIX XX
C
      GO TO 235
   60 CONTINUE
      NVX=(NV*(NV+1))/2
      CALL CHLSKY
C
C     INVERT THE LOWER TRIANGULAR MATRIX Q RESULTING FROM THE
C     CHOLESKY FACTORIZATION
C
      CALL MATINV
      IF(ICM.EQ.0.AND.(N.NE.3.OR.NV.NE.2))GO TO 160
C
C     IF THE USER HAS SPECIFIED HIS OWN RANK CORRELATION STRUCTURE
C     THIS MATRIX MULTIPLICATION LOOP IS NEEDED
C
      REWIND 3
      READ(3)CORR
      DO 150 K=1,NV
      KX=K
      DO 140 IT=1,KX
      XV(IT)=0.0
      ITX=IT
      DO 140 IW=ITX,KX
  140 XV(IT)=CORR(LOC1(K,IW))*Q(LOC1(IW,IT))+XV(IT)
      DO 150 IW=1,NV
  150 S(LOC1(K,IW))=XV(IW)
      GO TO 180
C
C     THIS ASSUMES THE DESIRED CORRELATION STRUCTURE IS THE IDENTITY
C     MATRIX, I.E. CLOSE TO ORTHOGONAL
C
  160 DO 170 I=1,NVX
  170 S(I)=Q(I)
  180 CONTINUE
      DO 210 K=1,N
      DO 190 IT=1,NV
      XV(IT)=0.0
      ITX=IT
      DO 190 IW=1,ITX
  190 XV(IT)=X(LOC(K,IW))*S(LOC1(IT,IW))+XV(IT)
      DO 200 IW=1,NV
  200 X(LOC(K,IW))=XV(IW)
  210 CONTINUE
C
C     SINCE X NO LONGER CONTAINS INTEGERS WE MUST RANK EACH COLUMN
C     IN ORDER TO GET A NEW MATRIX TO USE FOR MIXING THE INTERVALS
C     OF THE RAW DATA MATRIX XX
C
  235 CONTINUE
      REWIND 4
      READ(4)XX
      DO 230 J=1,NV
      DO 220 I=1,N
  220 XV(I)=X(LOC(I,J))
      CALL RANKER
      DO 230 I=1,N
  230 X(LOC(I,J))=RXV(I)
      DO 240 J=1,NV
      DO 240 I=1,N
      K=INT(X(LOC(I,J))+0.01)
  240 X(LOC(I,J))=XX(LOC(K,J))
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION NEXTBL(LINE,KI,KF)
C   RETURNS K WHERE LINE(K:K) IS THE FIRST BLANK IN LINE(KI:KF).
C   RETURNS 0 IF THERE ARE NO BLANKS IN LINE(KI:KF).

C   REVISION LOG:
C    09/19/90  HJI  ORIGINAL VERSION.

      INTEGER*4 NEXTBL, INDEX, KI, KF
      CHARACTER*(*) LINE

      NEXTBL=INDEX(LINE(KI:KF),' ')
      IF (NEXTBL .GT. 0) NEXTBL=NEXTBL+KI-1
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION NEXTNB(LINE,KI,KF)
C   RETURNS K WHERE LINE(K:K) IS THE FIRST NON-BLANK IN LINE(KI:KF).
C   RETURNS 0 IF LINE(KI:KF) IS BLANK.

C   REVISION LOG:
C    09/19/90  HJI  ORIGINAL VERSION.

      INTEGER*4 NEXTNB, K, KI, KF
      CHARACTER*(*) LINE

      NEXTNB=0
      DO 10 K=KI,KF
      IF (LINE(K:K) .NE. ' ') THEN
	 NEXTNB=K
	 GO TO 20
      ENDIF
   10 CONTINUE
   20 RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NORMAL(J,IDT)
C***********************************************************************
C SUBROUTINE NORMAL GENERATES NORMAL AND LOGNORMAL DISTRIBUTIONS
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 PROBINC, A, B, PMU, SIG, STRTPT, R, FINVNOR
      REAL*4 RNUM
      INTEGER*4 J, IDT, LOC, I
      LOC(I,J)=(J-1)*N+I
      PROBINC=1./FLOAT(N)
      IF(IRS.EQ.1)PROBINC=1.0
C
C     A IS ASSUMED TO BE THE LOWER .01 QUANTILE OF THE NORMAL
C     DISTRIBUTION.
C     B IS ASSUMED TO BE THE UPPER .01 QUANTILE OF THE NORMAL
C     DISTRIBUTION.
C
      READ(8)A,B
      IF(IDT.EQ.3)THEN
	A=LOG(A)
	B=LOG(B)
      ENDIF
      PMU=(A+B)/2.
c*********change from 99% to 99.9%********************
c      SIG=(B-PMU)/3.09023
      SIG=(B-PMU)/2.3263478741
c***************************************************
      STRTPT=0.
      DO 10 I=1,N
        RNUM=RAN(ISEED)
	R=PROBINC*RNUM+STRTPT
c************change from 99% to 99.9%*******************
c      IF(R.GE.0.999)R=0.999
c	IF(R.LE.0.001)R=0.001
	IF(R.GE.0.99)R=0.99
	IF(R.LE.0.01)R=0.01
c*******************************************************
	X(LOC(I,J))=FINVNOR(R)*SIG+PMU
	IF(IDT.EQ.3)X(LOC(I,J))=EXP(X(LOC(I,J)))
	IF(IRS.EQ.0)STRTPT=STRTPT+PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   10 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE OUTCRD(CARD)
C***********************************************************************
C SUBROUTINE OUTCRD PROCESSES THE OUTPUT PARAMETER OPTIONS
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (LENC=80)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP

      INTEGER*4 IC, IE
      CHARACTER CARD*(LENC),PDATA*4,PHIST*4,PCORR*4,PSPLAT*4,BLANK
      PARAMETER (PDATA='DATA',PHIST='HIST',PCORR='CORR',PSPLAT='PLOT',
     1           BLANK=' ')
      IC=7
   10 CONTINUE
      IC=IC+1
      IF(IC.GT.LENC)GO TO 20
      IF(CARD(IC:IC).EQ.BLANK)GO TO 10
      IE=IC+3
      IF(CARD(IC:IE).EQ.PDATA)THEN
	IDATA=1
	IC=IE+1
	GO TO 10
      ELSE IF(CARD(IC:IE).EQ.PHIST)THEN
	IHIST=1
	IC=IE+1
	GO TO 10
      ELSE IF(CARD(IC:IE).EQ.PCORR)THEN
	ICORR=1
	IC=IE+1
	GO TO 10
      ELSE IF(CARD(IC:IE).EQ.PSPLAT)THEN
	ISPLAT=1
        OPEN (UNIT=20,STATUS='NEW')
	IC=IE+1
	GO TO 10
      ELSE
	WRITE(6,9001)CARD
	RETURN
      ENDIF
   20 CONTINUE
      RETURN
 9001 FORMAT('1',5X,'THE FOLLOWING OUTPUT OPTION CARD REQUESTED ',
     1       'AN UNDEFINED OUTPUT OPTION',/,6X,'PLEASE CHECK THE ',
     2       'USER MANUAL FOR THE CORRECT OUTPUT OPTION CARD ',
     3       'SYNTAX',//,3X,'***',A,'***')
      END
      SUBROUTINE OUTDAT(IRK)
C***********************************************************************
C SUBROUTINE OUTDAT OUTPUTS THE SAMPLE AND ITS RANKS
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP

      INTEGER*4 IRK, LOC, I, J, K, I1, I2, ITEMP
      CHARACTER LAB1*2, LAB2*1
      PARAMETER (LAB1='X(',LAB2=')')
      LOC(I,J)=(J-1)*N+I
      K=1
      I1=1
      I2=10
      ITEMP=NV/10
      IF (ITEMP*10 .NE. NV) ITEMP=ITEMP+1
  300 IF (K.EQ.ITEMP) GO TO 350
      IF (K.GT.1) GO TO 320
      WRITE(6,820) TITLE
      IF(IRS.EQ.0)THEN
	IF(IRK.EQ.0)THEN
	  WRITE(6,810)(LAB1,I,LAB2,I=I1,I2)
	ELSE
	  WRITE(6,850)(LAB1,I,LAB2,I=I1,I2)
	ENDIF
      ELSE
	IF(IRK.EQ.0)THEN
	  WRITE(6,815)(LAB1,I,LAB2,I=I1,I2)
	ELSE
	  WRITE(6,855)(LAB1,I,LAB2,I=I1,I2)
	ENDIF
      ENDIF
      DO 310 J=1,N
	IF(IRK.EQ.0)THEN
	  WRITE(6,840)J,(X(LOC(J,I)),I=I1,I2)
	ELSE
	  WRITE(6,860)J,(X(LOC(J,I)),I=I1,I2)
	ENDIF
  310 CONTINUE
      GO TO 340
  320 WRITE(6,820)TITLE
      IF(IRS.EQ.0)THEN
	IF(IRK.EQ.0)THEN
	  WRITE(6,830)(LAB1,I,LAB2,I=I1,I2)
	ELSE
	  WRITE(6,870)(LAB1,I,LAB2,I=I1,I2)
	ENDIF
      ELSE
	IF(IRK.EQ.0)THEN
	  WRITE(6,835)(LAB1,I,LAB2,I=I1,I2)
	ELSE
	  WRITE(6,875)(LAB1,I,LAB2,I=I1,I2)
	ENDIF
      ENDIF
      DO 330 J=1,N
	IF(IRK.EQ.0)THEN
	  WRITE(6,840)J,(X(LOC(J,I)),I=I1,I2)
	ELSE
	  WRITE(6,860)J,(X(LOC(J,I)),I=I1,I2)
	ENDIF
  330 CONTINUE
  340 K=K+1
      I1=I1+10
      I2=I2+10
      GO TO 300
  350 IF (K.GT.1) GO TO 370
      WRITE(6,820) TITLE
      IF(IRS.EQ.0)THEN
	IF(IRK.EQ.0)THEN
	  WRITE(6,810)(LAB1,I,LAB2,I=I1,NV)
	ELSE
	  WRITE(6,850)(LAB1,I,LAB2,I=I1,NV)
	ENDIF
      ELSE
	IF(IRK.EQ.0)THEN
	  WRITE(6,815)(LAB1,I,LAB2,I=I1,NV)
	ELSE
	  WRITE(6,855)(LAB1,I,LAB2,I=I1,NV)
	ENDIF
      ENDIF
      DO 360 J=1,N
	IF(IRK.EQ.0)THEN
	  WRITE(6,840)J,(X(LOC(J,I)),I=I1,NV)
	ELSE
	  WRITE(6,860)J,(X(LOC(J,I)),I=I1,NV)
	ENDIF
  360 CONTINUE
      GO TO 390
  370 WRITE(6,820) TITLE
      IF(IRS.EQ.0)THEN
	IF(IRK.EQ.0)THEN
	  WRITE(6,830)(LAB1,I,LAB2,I=I1,NV)
	ELSE
	  WRITE(6,870)(LAB1,I,LAB2,I=I1,NV)
	ENDIF
      ELSE
	IF(IRK.EQ.0)THEN
	  WRITE(6,835)(LAB1,I,LAB2,I=I1,NV)
	ELSE
	  WRITE(6,875)(LAB1,I,LAB2,I=I1,NV)
	ENDIF
      ENDIF
      DO 380 J=1,N
	IF(IRK.EQ.0)THEN
	  WRITE(6,840)J,(X(LOC(J,I)),I=I1,NV)
	ELSE
	  WRITE(6,860)J,(X(LOC(J,I)),I=I1,NV)
	ENDIF
  380 CONTINUE
  390 CONTINUE
      RETURN
  810 FORMAT('0LATIN HYPERCUBE SAMPLE INPUT VECTORS',//,' RUN NO.',
     1       2X,A,I1,A,8(7X,A,I1,A),6X,A,I2,A)
  815 FORMAT('0RANDOM SAMPLE INPUT VECTORS',//,' RUN NO.',
     1       2X,A,I1,A,8(7X,A,I1,A),6X,A,I2,A)
  820 FORMAT('1',A100)
  830 FORMAT('0LATIN HYPERCUBE SAMPLE INPUT VECTORS',//,' RUN NO.',
     1       3X,A,I2,A,9(6X,A,I2,A))
  835 FORMAT('0RANDOM SAMPLE INPUT VECTORS',//,' RUN NO.',
     1       3X,A,I2,A,9(6X,A,I2,A))
  840 FORMAT('0',I5,1P,10E11.3)
  850 FORMAT('0RANKS OF LATIN HYPERCUBE SAMPLE INPUT VECTORS',//,
     1       ' RUN NO.',5X,A,I1,A,8(7X,A,I1,A),6X,A,I2,A)
  855 FORMAT('0RANKS OF RANDOM SAMPLE INPUT VECTORS',//,' RUN NO.',
     1        5X,A,I1,A,8(7X,A,I1,A),6X,A,I2,A)
  860 FORMAT('0',I5,10F11.0)
  870 FORMAT('0RANKS OF LATIN HYPERCUBE SAMPLE INPUT VECTORS',//,
     1       ' RUN NO.',5X,A,I2,A,9(6X,A,I2,A))
  875 FORMAT('0RANKS OF RANDOM SAMPLE INPUT VECTORS',//,' RUN NO.',
     1        5X,A,I2,A,9(6X,A,I2,A))
      END
      SUBROUTINE PMTRX(NC,IFLAG)
C***********************************************************************
C SUBROUTINE PMTRX PRINTS OUT A CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (NROW=28)
      PARAMETER (NCOL=15)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM

      INTEGER*4 NC, IFLAG, LOC1, I, J, I1, I2, NPROW, NPCOL, NR
      INTEGER*4 J1, J2, IC, K, K1, II
      LOC1(I,J)=J+(I*I-I)/2
      I=1
      I1=1
      I2=NROW
      NPROW=NC/NROW
      IF (NPROW*NROW.NE.NC) NPROW=NPROW+1
      NPCOL=NC/NCOL
      IF (NPCOL*NCOL.NE.NC) NPCOL=NPCOL+1
      DO 100 NR=1,NPROW
      J1=1
      J2=NCOL
      DO 90 IC=1,NPCOL
      IF (J1.GT.I2) GO TO 90
      IF (IFLAG.EQ.2) GO TO 10
      IF (IFLAG.EQ.3) GO TO 20
      IF(IFLAG.EQ.4)GO TO 20
      WRITE(6,180) TITLE
      IF(IRS.EQ.0)WRITE(6,110)I
      IF(IRS.NE.0)WRITE(6,115)I
      GO TO 30
   10 WRITE(6,180) TITLE
      IF(IRS.EQ.0)WRITE(6,120)I
      IF(IRS.NE.0)WRITE(6,125)I
      GO TO 30
   20 WRITE(6,180) TITLE
      IF(IFLAG.EQ.3)WRITE(6,130)I
      IF(IFLAG.EQ.4)WRITE(6,190)I
   30 IF (NC.LT.I2) I2=NC
      IF (NC.LT.J2) J2=NC
      DO 60 K=I1,I2
      K1=K
      IF (K.GE.J1) GO TO 40
      WRITE(6,140)
      GO TO 60
   40 IF (K.GT.J2) GO TO 50
      WRITE(6,150)LCM(K),(CORR(LOC1(K,J)),J=J1,K)
      GO TO 60
   50 WRITE(6,150)LCM(K),(CORR(LOC1(K,J)),J=J1,J2)
   60 CONTINUE
      IF (K1.LT.J2) GO TO 70
      WRITE(6,160)(LCM(II),II=J1,J2)
      WRITE(6,170)
      GO TO 80
   70 WRITE(6,160)(LCM(II),II=J1,K1)
      WRITE(6,170)
   80 CONTINUE
      J1=J1+NCOL
      J2=J2+NCOL
      I=I+1
   90 CONTINUE
      I1=I1+NROW
      I2=I2+NROW
  100 CONTINUE
      RETURN
  110 FORMAT('0CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE LATIN',
     1       ' HYPERCUBE SAMPLE FOR RAW DATA',30X,'PAGE ',I3)
  115 FORMAT('0CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'RANDOM SAMPLE FOR RAW DATA',40X,'PAGE ',I3)
  120 FORMAT('0CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE LATIN',
     1       ' HYPERCUBE SAMPLE FOR RANK DATA',30X,'PAGE ',I3)
  125 FORMAT('0CORRELATIONS AMONG INPUT VARIABLES CREATED BY THE ',
     1       'RANDOM SAMPLE FOR RANK DATA',40X,'PAGE ',I3)
  130 FORMAT('0INPUT RANK CORRELATION MATRIX',95X,'PAGE ',I3)
  140 FORMAT('0')
  150 FORMAT('0',I5,15(F8.4))
  160 FORMAT('0',5X,15I8)
  170 FORMAT('0VARIABLES')
  180 FORMAT('1',A100)
  190 FORMAT('0ADJUSTED RANK CORRELATION MATRIX',92X,'PAGE ',I3)
      END
      SUBROUTINE POSDEF(ITEST)
C***********************************************************************
C SUBROUTINE POSDEF IS USED IN THE POSITIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM
      COMMON/PDMAT/Z, D
      REAL*4       Z(NVAR,NVAR), D(NVAR)

      REAL*4 WK((NVAR*(NVAR+1))/2), EIG
      INTEGER*4 ITEST, NN, M, NP, ICONV, INFO
C
C     NN = THE DIMENSIONS OF Z,  I.E.  Z(NN,NN).
      DATA NN/NVAR/
C
C
C     EIG = THE VALUE THAT THE NEGATIVE EIGENVALUES ARE SET EQUAL TO.
      DATA EIG/0.001/
C
C
C     M = MAXIMUM NUMBER OF ITERATIONS ALLOWED.
      M=20
C
      NP=NCM
      ITEST=0
      ICONV=0
  100 CONTINUE
      ITEST=ITEST+1
      IF(ITEST.GT.M)THEN
	WRITE(6,1000)
	RETURN
      ENDIF
      REWIND 3
      WRITE(3)CORR
      CALL SSPEV(CORR,NP,D,Z,NN,WK,1,INFO)
      CALL FINDIT(NP,NN,EIG,ICONV)
      IF(ICONV.EQ.0)GO TO 100
      REWIND 3
      READ(3)CORR
      RETURN
 1000 FORMAT('1','THE INPUT RANK CORRELATION MATRIX IS NOT POSITIVE ',
     1       'DEFINITE.',/,' AN ITERATIVE PROCEDURE HAS FAILED TO ',
     2       'PRODUCE A POSITIVE DEFINITE MATRIX AFTER 20 ITERATIONS.',
     3       /,' THEREFORE, THE PROGRAM HAS BEEN TERMINATED.',/,' THE',
     4       ' USER NEEDS TO REEVALUATE THE RELATIONSHIP OF THE ',
     5       'CORRELATED VARIABLES IN THE MATRIX.')
      END
      FUNCTION PYTHAG(A,B)
C***********************************************************************
C FUNCTION PYTHAG IS USED IN THE POSTIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  PYTHAG
C***REFER TO  EISDOC
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C***ROUTINES CALLED    (NONE)
C***END PROLOGUE  PYTHAG
C***********************************************************************
      REAL*4 PYTHAG, A, B
C
      REAL*4 P, Q, R, S, T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
	 R = (Q/P)**2
	 T = 4.0E0 + R
	 IF (T .EQ. 4.0E0) GO TO 20
	 S = R/T
	 P = P + 2.0E0*P*S
	 Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      FUNCTION R1MACH(I)
C***********************************************************************
C FUNCTION R1MACH IS USED TO SET FLOATING POINT MACHINE
C DEPENDENT CONSTANTS
C
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  801001   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns single precision machine dependent constants
C***DESCRIPTION
C     R1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          A = R1MACH(I)
C
C     where I=1,...,5.  The (output) value of A above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Single-Precision Machine Constants
C  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  R1MACH(3) = B**(-T), the smallest relative spacing.
C  R1MACH(4) = B**(1-T), the largest relative spacing.
C  R1MACH(5) = LOG10(B)
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***END PROLOGUE  R1MACH
C***********************************************************************
      INTEGER*4 I, SMALL(2), LARGE(2), RIGHT(2), DIVER(2), LOG10(2)
C
      REAL*4 R1MACH, RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C    (EXPRESSED IN INTEGER AND HEXADECIMAL)
C  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
C  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
      DATA SMALL(1) /       128 /
      DATA LARGE(1) /    -32769 /
      DATA RIGHT(1) /     13440 /
      DATA DIVER(1) /     13568 /
      DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
C
      R1MACH = RMACH(I)
      RETURN
C
      END
      SUBROUTINE RANKER
C***********************************************************************
C SUBROUTINE RANKER IS USED TO FIND THE RANKS OF A VECTOR OF DATA
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/RANK/XV, RXV, IWK
      REAL*4      XV(NMAX),RXV(NMAX)
      INTEGER*4   IWK(NMAX)
      REAL*4 AVG, XHOLD, RHOLD
      INTEGER*4 I, NTIES, II, J, K
      DO 10 I=1,N
   10 RXV(I)=FLOAT(I)
      CALL HPSRT
      DO 20 I=1,N
      IWK(I)=IFIX(RXV(I))
   20 RXV(I)=FLOAT(I)
C
C     FIND TIES
      I=0
   30 I=I+1
      IF (I.GE.N) GO TO 80
      IF (XV(I).NE.XV(I+1)) GO TO 30
C
C     COUNT TIES
      NTIES=2
   40 II=I+NTIES
      IF (II.GT.N) GO TO 50
      IF (XV(I).NE.XV(II)) GO TO 50
      NTIES=NTIES+1
      GO TO 40
C
C     AVERAGE TIED RANKS
   50 AVG=0.0
      DO 60 J=1,NTIES
   60 AVG=AVG+RXV(I+J-1)
      AVG=AVG/FLOAT(NTIES)
      I=I-1
      DO 70 J=1,NTIES
      I=I+1
   70 RXV(I)=AVG
      GO TO 30
   80 CONTINUE
C
C     REORDER
      I=1
   90 K=IWK(I)
      IF (K.NE.I) GO TO 100
      I=I+1
      IF (I.LT.N) GO TO 90
      RETURN
  100 XHOLD=XV(I)
      RHOLD=RXV(I)
      XV(I)=XV(K)
      RXV(I)=RXV(K)
      XV(K)=XHOLD
      RXV(K)=RHOLD
      IWK(I)=IWK(K)
      IWK(K)=K
      GO TO 90
      END
      SUBROUTINE RAYLEIG(J)
C***********************************************************************
C SUBROUTINE RAYLEIG GENERATES THE RAYLEIGH DISTRIBUTION
C
C    f(X) = C*(X-A)/(B-A)*EXP(-.5*C*(X-A)**2/(B-A))/(1.-EXP(-.5*C*(B-A))
C                                                           A <= X <= B.
C    F(X) = (1.-EXP(-0.5*C*(X-A)**2/(B-A)))/(1.-EXP(-0.5*C*(B-A)))
C
C   REVISION LOG:
C    2/25/93  HJI  REPLACED C BY C/(B-A) SO THAT C CAN BE THE MAXIMUM
C                  INTRUSION RATE (LAMBDA) INSTEAD OF THE RATE OF CHANGE
C                  OF THE INTRUSION RATE.
C***********************************************************************

      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 LAMBDA, PROBINC, A, B, DENOM, STRTPT, R
      REAL*4 RNUM
      INTEGER*4 J, LOC, I
      LOC(I,J)=(J-1)*N+I
      PROBINC=1./FLOAT(N)
      IF(IRS.EQ.1)PROBINC=1.0
      READ(8)A,B,LAMBDA

C   SELECT A NORMALIZATION CONDITION:
CXC   THE FOLLOWING SETS THE NORMALIZATION TO F(B)=1:
CX      DENOM=1.-EXP(-0.5*LAMBDA*(B-A))
C   THE FOLLOWING SETS THE NORMALIZATION TO F(infinity)=1:
      DENOM=1.

      STRTPT=0.
      DO 10 I=1,N
        RNUM=RAN(ISEED)
	R=PROBINC*RNUM+STRTPT
	X(LOC(I,J))=A+SQRT(-2.*(B-A)*ALOG(1.-R*DENOM)/LAMBDA)
	IF(IRS.EQ.0)STRTPT=STRTPT+PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   10 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE RAYLEXP(J)
C***********************************************************************
C SUBROUTINE RAYLEXP GENERATES THE RAYLEIGH-EXPONENTIAL DISTRIBUTION
C
C FOR A<=X<=B,  THE RAYLEIGH DISTRIBUTION IS USED:
C      f(X) = LAMBDA*(X-A)*EXP(-0.5*LAMBDA*(X-A)**2))
C    AND
C      F(X) = 1 - EXP(-0.5*LAMBDA*(X-A)**2)
C
C FOR B<=X<=C,  THE EXPONENTIAL DISTRIBUTION IS USED:
C      f(X) = LAMBDA*(B-A)*EXP( -0.5*LAMBDA*(B-A)*(X-0.5*(A+B)) )
C    AND
C      F(X) = 1 - EXP( -0.5*LAMBDA*(B-A)*(X-0.5*(A+B)) )
C
C***********************************************************************
      INTEGER*4 I, J, LOC, NMAX, NVAR
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      REAL*4 R, A, B, FB, LAMBDA, PROBINC, STRTPT
      REAL*4 RNUM

      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      LOC(I,J)=(J-1)*N+I

      PROBINC=1./FLOAT(N)
      IF (IRS.EQ.1) PROBINC=1.0
      READ(8) A, B, STRTPT, LAMBDA
      FB=1.-EXP(-0.5*LAMBDA*(B-A)**2)
      STRTPT=0.
      DO 10 I=1,N
        RNUM=RAN(ISEED)
	R=PROBINC*RNUM+STRTPT
	IF (R .LE. FB) THEN
	   X(LOC(I,J))=A+SQRT(-2.*ALOG(1.-R)/LAMBDA)
	ELSE
	   X(LOC(I,J))=0.5*(A+B)-ALOG(1.-R)/(LAMBDA*(B-A))
	ENDIF
	IF(IRS.EQ.0)STRTPT=STRTPT+PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   10 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE RDPAR
C***********************************************************************
C SUBROUTINE RDPAR PROCESSES THE PARAMETER STATEMENTS, DEFINES
C VARIABLES IN COMMON /PARAM/ AND STORES DISTRIBUTION INFORMATION
C ON UNITS 7, 8 AND 9
C***********************************************************************
      INTEGER*4  LENC, LENTC, NCVAR, NINTMX, NMAX, NVAR
      PARAMETER (NSTUD=100)
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
c*******************************
      PARAMETER (NCVAR=1000)
c******************************      
      PARAMETER (NINTMX=50)
      PARAMETER (LENC=80)
      PARAMETER (LENTC=11)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/STAR/NSUBOB, NINT, SUBINT
      INTEGER*4   NSUBOB(NINTMX), NINT
      REAL*4      SUBINT(NINTMX+1)
      COMMON/UICORR/ICVAR, JCVAR, CVAR, NCV
      INTEGER*4     ICVAR(NCVAR), JCVAR(NCVAR), NCV
      REAL*4        CVAR(NCVAR)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP

      REAL*4 AVG, FNORM, STDDEV, XVAL(500), A, B, P, Q, C, D
      INTEGER*4 I, ININE, IP, IRSET, ITEST, IZERO, NDF, NP, K
      INTEGER*4 JSTUD
      COMMON/STU/ STUMAX, STUMIN, STUDMX, STUDMN
      REAL*4      STUMAX, STUMIN, STUDMX(NSTUD), STUDMN(NSTUD)
      LOGICAL  LOGDIST
      CHARACTER*1  CH1
      CHARACTER CARD*(LENC), TMPCRD*(LENTC), CDUM*1
      CHARACTER PTITLE*6,PRSAMP*14,PRSEED*12,PNOBS*5,LOBS*12,PNREP*6,
     1          PC*4,PCMINP*18,LCMINP*29,POUT*7,PBETA*5,PNOR*7,PLNOR*10,
     2          PUNIS*9,PLUNIS*12,PUSDST*18,BLANK*1,
     3          PERIOD*1,MINUS*1,PIRP*15,PEXPONT*12,PRAYLEI*9
      CHARACTER*8  PRAYLEX, PSTUDENT, PUNI
      CHARACTER*11 PLSTUD, PLUNI, PTRIAG
      CHARACTER*19 IRSAMP
      PARAMETER (PTITLE='TITLE ',PRSAMP='RANDOM SAMPLE ',PC='CORR',
     1           PRSEED='RANDOM SEED ',PNOBS='NOBS ',
     2           PNREP='NREPS ',PCMINP='CORRELATION MATRIX',
     3           LCMINP='PAIRS OF CORRELATED VARIABLES',
     4           PBETA='BETA ',PNOR='NORMAL ',PLNOR='LOGNORMAL ',
     5           PUNI='UNIFORM ',PLUNI='LOGUNIFORM ',
     6           PLUNIS='LOGUNIFORM* ',PUSDST='USER DISTRIBUTION ',
     7           BLANK=' ',PERIOD='.',MINUS='-',PTRIAG='TRIANGULAR ',
     8           LOBS='OBSERVATIONS',DUM1=9999.0,DUM2=9999.0,
     9           POUT='OUTPUT ',PUNIS='UNIFORM* ',
     A           PIRP='RANDOM PAIRING ', PEXPONT='EXPONENTIAL ',
     B           PRAYLEI='RAYLEIGH ', PRAYLEX='RAYLEXP ',
     C           PSTUDENT='STUDENT ', PLSTUD='LOGSTUDENT ')
      DATA IRSAMP/'***RANDOM SAMPLE***'/
      DATA CDUM / ' ' /
      IZERO=ICHAR('0')
      ININE=ICHAR('9')
      IRSET=0
      JSTUD=1
C SET DEFAULT VALUES FOR PARAMETERS
      CALL SETDEF
C READ PARAMETER CARD
   10 CONTINUE
      READ(5,9001,END=8000)CARD
C TITLE CARD
   20 CONTINUE
      IF(CARD(1:6).EQ.PTITLE)THEN
	READ(CARD,9002) TITLE
	GO TO 10
C RANDOM SAMPLE CARD
      ELSE IF(CARD(1:14).EQ.PRSAMP)THEN
	READ(IRSAMP,9001) TITLE(82:100)
	IRS=1
	GO TO 10
C RANDOM SEED CARD
      ELSE IF(CARD(1:12).EQ.PRSEED)THEN
	CALL DATSQZ(CARD(13:80),PRSEED,TMPCRD)
	READ(TMPCRD,9003)ISEED
C SET RANDOM NUMBER GENERATOR HERE IF NECESSARY
	IRSET=1
	GO TO 10
C RANDOM PAIRING CARD
      ELSE IF(CARD(1:15).EQ.PIRP)THEN
	IRP=1
	GO TO 10
C NUMBER OF OBSERVATIONS CARD
      ELSE IF(CARD(1:5).EQ.PNOBS)THEN
	CALL DATSQZ(CARD(6:80),PNOBS,TMPCRD)
	READ(TMPCRD,9003)N
	IF(N.LT.1)THEN
	  WRITE(6,9006)N
	  RETURN
	ENDIF
	CALL CHKDIM(1,N,NMAX,PNOBS,LOBS)
	GO TO 10
C NUMBER OF REPETITIONS CARD
      ELSE IF(CARD(1:6).EQ.PNREP)THEN
	CALL DATSQZ(CARD(7:80),PNREP,TMPCRD)
	READ(TMPCRD,9003)NREP
	IF(NREP.LT.1)THEN
	  WRITE(6,9007)NREP
	  RETURN
	ENDIF
	GO TO 10
C CORRELATION MATRIX CARD
      ELSE IF(CARD(1:4).EQ.PC)THEN
	READ(5,*,ERR=9000)NCV,
     1                    (ICVAR(I),JCVAR(I),CVAR(I),I=1,NCV)
	CALL CHKDIM(1,NCV,NCVAR,PCMINP,LCMINP)
	ICM=1
	GO TO 10
C OUTPUT OPTIONS CARD
      ELSE IF(CARD(1:7).EQ.POUT)THEN
	CALL OUTCRD(CARD)
	GO TO 10
C BETA DISTRIBUTION CARD
      ELSE IF(CARD(1:5).EQ.PBETA)THEN
	READ(5,*,ERR=9000)A,B,P,Q
	CALL CHKDAT(PBETA,A,B,P,Q)
	CALL WRTCRD(1,CARD(6:80))
	GO TO 10
C NORMAL DISTRIBUTION CARD
      ELSE IF(CARD(1:7).EQ.PNOR)THEN
	READ(5,*,ERR=9000)A,B
	CALL CHKDAT(PNOR,A,B,DUM1,DUM2)
	CALL WRTCRD(2,CARD(8:80))
	GO TO 10
C LOGNORMAL DISTRIBUTION CARD
      ELSE IF(CARD(1:10).EQ.PLNOR)THEN
	READ(5,*,ERR=9000)A,B
	CALL CHKDAT(PLNOR,A,B,A,B)
	CALL WRTCRD(3,CARD(11:80))
	GO TO 10
C UNIFORM DISTRIBUTION CARD
      ELSE IF(CARD(1:8).EQ.PUNI)THEN
	READ(5,*,ERR=9000)A,B
	CALL CHKDAT(PUNI,A,B,DUM1,DUM2)
	CALL WRTCRD(4,CARD(9:80))
	GO TO 10
C LOGUNIFORM DISTRIBUTION CARD
      ELSE IF(CARD(1:11).EQ.PLUNI)THEN
	READ(5,*,ERR=9000)A,B
	CALL CHKDAT(PLUNI,A,B,A,B)
	CALL WRTCRD(5,CARD(12:80))
	GO TO 10
C UNIFORM* DISTRIBUTION CARD
      ELSE IF(CARD(1:9).EQ.PUNIS)THEN
	READ(5,*,ERR=9000)NINT,(NSUBOB(I),I=1,NINT),
     1                    (SUBINT(I),I=1,NINT+1)
	CALL CHKSTR(PUNIS,CARD)
	CALL WRTCRD(6,CARD(10:80))
	GO TO 10
C LOGUNIFORM* DISTRIBUTION CARD
      ELSE IF(CARD(1:12).EQ.PLUNIS)THEN
	READ(5,*,ERR=9000)NINT,(NSUBOB(I),I=1,NINT),
     1                    (SUBINT(I),I=1,NINT+1)
	CALL CHKSTR(PLUNIS,CARD)
	CALL WRTCRD(7,CARD(13:80))
	GO TO 10
C TRIANGULAR DISTRIBUTION CARD
      ELSE IF(CARD(1:11).EQ.PTRIAG)THEN
	READ(5,*,ERR=9000)A,B,C
	CALL CHKTRI(A,B,C)
	CALL WRTCRD(8,CARD(12:80))
	GO TO 10
C USER DISTRIBUTION/DATA CARD(S)
      ELSE IF(CARD(1:18).EQ.PUSDST)THEN
	CALL WRTCRD(9,CARD(19:80))
   30   CONTINUE
	READ(5,9001,END=8000)CARD
	CH1=CARD(1:1)
	ITEST=ICHAR(CH1)
	IF (CH1 .EQ. BLANK .OR. CH1 .EQ. PERIOD. OR. CH1 .EQ. MINUS
     1    .OR. (ITEST .GE. IZERO .AND. ITEST .LE. ININE)) THEN
	  WRITE(7,9001)CARD
	  GO TO 30
	ELSE
	  GO TO 20
	ENDIF
C EXPONENTIAL DISTRIBUTION CARD
      ELSE IF(CARD(1:12).EQ.PEXPONT)THEN
	READ(5,*,ERR=9000)A,B,C
	CALL CHKEXP(A,B,C)
	CALL WRTCRD(10,CARD(12:80))
	GO TO 10
C RAYLEIGH DISTRIBUTION CARD
      ELSE IF(CARD(1:9).EQ.PRAYLEI) THEN
	READ(5,*,ERR=9000)A,B,C
	CALL CHKRAY(A,B,C)
	CALL WRTCRD(11,CARD(10:80))
	GO TO 10
C RAYLEIGH-EXPONENTIAL DISTRIBUTION CARD
      ELSE IF(CARD(1:8).EQ.PRAYLEX) THEN
	READ (5,*,ERR=9000) A,B,C,D
	CALL CHKREX(A,B,C,D)
	CALL WRTCRD(12,CARD(9:80))
	GO TO 10
C STUDENT DISTRIBUTION & DATA CARD(S)
      ELSE IF (CARD(1:8) .EQ. PSTUDENT .OR. CARD(1:11) .EQ. PLSTUD) THEN
	 IF (CARD(1:8) .EQ. PSTUDENT) THEN
	    CALL WRTCRD(13,CARD(9:80))
	    LOGDIST=.FALSE.
	 ELSE
	    CALL WRTCRD(14,CARD(12:80))
	    LOGDIST=.TRUE.
	 ENDIF
	 READ (5,9001,END=8000) CARD
C   Read the number of data points (NP).  This is 1 larger than the
C   number of degrees of freedom.  The average and standard deviation,
C   if available, should on the same line as the number of points.
	 IP=0
	 STDDEV=0.
	 READ (CARD,*,END=40) NP,AVG,STDDEV
   40    CONTINUE
	 READ (5,9001,END=8000) CARD
	 CH1=CARD(1:1)
	 ITEST=ICHAR(CH1)
	 IF (CH1 .EQ. BLANK .OR. CH1 .EQ. PERIOD. OR. CH1 .EQ. MINUS
     &       .OR. (ITEST .GE. IZERO .AND. ITEST .LE. ININE)) THEN
C   Read the sample values from the input line into the XVAL array.
C   IP is the count of values so far and is updated by READVALS.
	    CALL READVALS(CARD,XVAL,IP)
	    GO TO 40
	 ELSE
	    IF (STDDEV .EQ. 0.) THEN
C   Compute the average and standard deviation.
	       AVG=0.0
	       FNORM=0.
               STUMAX=-1.0E+35
               STUMIN=1.0E+35
	       DO 50 K=1,NP
	       IF (LOGDIST) XVAL(K)=ALOG(XVAL(K))
               IF (XVAL(K) .GT. STUMAX) STUMAX=XVAL(K)
               IF (XVAL(K) .LT. STUMIN) STUMIN=XVAL(K)
	       IF (ABS(XVAL(K)) .GT. FNORM) FNORM=ABS(XVAL(K))
   50          AVG=AVG+XVAL(K)
	       AVG=AVG/NP
               STUDMN(JSTUD)=STUMIN
               STUDMX(JSTUD)=STUMAX
C   Scaling by FNORM is used to prevent underflow/overflow.
	       DO 60 K=1,NP
   60          STDDEV=STDDEV+((XVAL(K)-AVG)/FNORM)**2
	       STDDEV=FNORM*SQRT(STDDEV/(NP*(NP-1.)))
	    ENDIF
	    NDF=NP-1
            JSTUD=JSTUD+1
	    WRITE (8) NDF,AVG,STDDEV
	    GO TO 20
	 ENDIF
C UNDEFINED PARAMETER/DATA CARD
      ELSE
	WRITE(6,9004)CARD
	RETURN
      ENDIF
C CHECK FOR INPUT ERRORS
 8000 CONTINUE
      CALL CHKZRO(N,NV,IRSET)
      CALL CHKDIM(2,NV,NVAR,CDUM,CDUM)
      IF(ICM.EQ.1)CALL CMCRD
      RETURN
 9000 CONTINUE
      WRITE(6,9005)CARD
      RETURN
 9001 FORMAT(A)
 9002 FORMAT(A80)
 9003 FORMAT(I11)
 9004 FORMAT('1',5X,'THE FOLLOWING CARD (POSSIBLY BLANK) IS NOT A ',
     1       'VALID PARAMETER/DATA CARD',/,6X,'PLEASE CONSULT THE ',
     2       'USER MANUAL FOR THE CORRECT PARAMETER/DATA CARD SYNTAX',
     3        //,3X,'***',A,'***')
 9005 FORMAT('1',5X,'THE FOLLOWING PARAMETER CARD DID NOT HAVE THE ',
     1       'CORRECT DATA CARD ASSOCIATED WITH IT',/,6X,'PLEASE ',
     2       'CONSULT THE USER MANUAL FOR THE CORRECT DATA ',
     3       'CARD SYNTAX',/,3X,'***',A,'***')
 9006 FORMAT('1',5X,'THE NUMBER OF OBSERVATIONS REQUESTED IS LESS ',
     1       'THAN ONE',I5)
 9007 FORMAT('1',5X,'THE NUMBER OF REPETITIONS REQUESTED IS LESS ',
     1       'THAN ONE',I5)
      END
      FUNCTION RIERFC1 (Y)
C***********************************************************************
C FUNCTION RIERFC1 IS USED IN GENERATING THE NORMAL AND
C LOGNORMAL DISTRIBUTIONS
C
C     THIS IS THE SAME AS THE RIERFC ROUTINE IN THE DEAMOS LIBRARY
C     THE NAME MODIFICATION WAS TO PREVENT OUR LIBRARY FROM CALLING
C     ANOTHER LIBRARY.                 ELF   OCTOBER 1980
C
C     WRITTEN BY D.E. AMOS AND S.L. DANIEL, SEPTEMBER, 1972.
C
C     REFERENCES
C         HASTINGS, C.JR., APPROXIMATIONS FOR DIGITAL COMPUTERS,
C         PRINCETON UNIVERSITY PRESS, PRINCETON, N.J., 1955
C
C         COMMUNICATION FROM L.F. SHAMPINE FOR CHEBYSHEV COEFFICIENTS.
C
C     ABSTRACT
C         RIERFC EVALUATES THE INVERSE COERROR FUNCTION DEFINED BY
C
C                          Y= ERFC(X)     0.LE.X.LT.INFINITY
C
C         WHERE 0.LT.Y.LE.1. CHEBYSHEV APPROXIMATIONS ON
C
C              EXP(-81).LE.Y.LT.0.1, 0.1.LE.Y.LT.0.5, 0.5.LE.Y.LE.1.
C
C         ARE USED WITH A CHANGE OF VARIABLES
C
C         YY=C1*W+C2, W=SQRT(-LN(Y)), YY=5.*Y-1.5, YY=2.*(1.-Y)
C
C         RESPECTIVELY. THE INVERSE OF THE NORMAL DISTRIBUTION IS GIVEN
C         BY
C                          SQRT(2)*RIERFC(2.*(1.-RN))  0.5.LE.RN.LT.1.0,
C                      X=
C                         -SQRT(2)*RIERFC(2.*RN)       0.0.LT.RN.LT.0.5.
C
C         THE RELATIVE ERROR IN RIERFC DECREASES FROM 1.E-10 TO 3.E-13
C         AS Y INCREASES FROM EXP(-81) TO 1.0.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT
C
C           Y      - Y, EXP(-81).LE.Y.LE.1.
C
C         OUTPUT
C
C           RIERFC - VALUE FOR THE INVERSE COERROR FUNCTION
C
C     ERROR CONDITIONS
C         Y.LT.EXP(-81) OR Y.GT.1 ARE FATAL ERRORS
C
C
C     RIERFC USES SUBROUTINES ERRCHK, ERRGET, ERRPRT, ERXSET, ERSTGT
C     COMPILE DECKS RIERFC, ERRCHK
C
C***********************************************************************
      REAL*4 A(22,3), A1(22), A2(22), A3(22)
      EQUIVALENCE (A(1,1),A1(1))
      EQUIVALENCE (A(1,2),A2(1))
      EQUIVALENCE (A(1,3),A3(1))
      REAL*4 RIERFC1, Y, C1, C2, W, D, TD, VNP1, VN, TEMP
      INTEGER*4 I, J, L, K
C
      DATA(A1(I),I=1,22)/9.18725611735013E-01,0.,1.68792878000327E-02
     1,0.,6.60337139058300E-04,0.,3.20203849839380E-05,0.,1.720
     260607522481E-06,0.,9.81965971588191E-08,0.,5.83049613537653E
     3-09,0.,3.56019351836136E-10,0.,2.21968915783128E-11,0.,1.
     440639693109741E-12,0.,9.02597345404862E-14,0./
C
      DATA(A2(I),I=1,22)/1.54701109458613E+00,-3.31460331083896E-01,4.33
     1001124090060E-02,-1.06564004165532E-02,2.90613542304156E-03,-8.618
     272838022491E-04,2.67933751795053E-04,-8.60838893942933E-05,2.83232
     3058814598E-05,-9.48870819734494E-06,3.22422655069385E-06,-1.108157
     478472076E-06,3.84464770797987E-07,-1.34439275565208E-07,4.73255976
     5052393E-08,-1.67556011100019E-08,5.96199003969093E-09,-2.130705032
     691886E-09,7.64427040920545E-10,-2.75198005584737E-10,9.93792246090
     7789E-11,-3.59877382902119E-11/
C
      DATA(A3(I),I=1,22)/1.10642888011036E+01,4.34299147561447E+00,-2.33
     1781774969295E-02,4.23345215362947E-03,8.68757084192089E-06,-5.9826
     21113270881E-04,4.50490139240298E-04,-2.54858131942102E-04,1.278241
     389261340E-04,-5.97873878043957E-05,2.66474012012582E-05,-1.1438183
     46209267E-05,4.75393030377615E-06,-1.91759589929610E-06,7.508064655
     594834E-07,-2.84791180387123E-07,1.04187791696225E-07,-3.6456724368
     69145E-08,1.20129296139030E-08,-3.61030126779729E-09,9.123561400817
     759E-10,-1.36851363400914E-10/
C
      DATA C1,C2/2.35777520630369E-01,1.35777520630369E+00/
C
      IF (Y.LT.6.63967719958073E-36.OR.Y.GT.1.0) GO TO 50
      IF (Y.GE.0.5) GO TO 10
      IF (Y.GE.0.1) GO TO 20
      J=3
      W=SQRT(-LOG(Y))
      D=C1*W-C2
      GO TO 30
   10 J=1
      D=1.-Y
      D=D+D
      GO TO 30
   20 J=2
      D=5.*Y-1.5
   30 TD=D+D
      VNP1=0.
      VN=0.
      DO 40 L=1,21
      K=22-L+1
      TEMP=VN
      VN=TD*VN-VNP1+A(K,J)
   40 VNP1=TEMP
      RIERFC1=D*VN-VNP1+.5*A(1,J)
      IF (J.EQ.1) RIERFC1=D*RIERFC1
      RETURN
   50 WRITE(6,60)
C      RIERFC1=0.
      RETURN
C      RETURN
C
   60 FORMAT(' Y LESS THAN EXP(-81.) OR Y GREATER THAN 1.0')
      END
      SUBROUTINE READVALS(INPLINE,XVAL,IP)
C   READ NUMERIC FIELDS FROM TEXT STRING INPLINE.
C   STORE THE RESULTS STARTING AT XVAL(IP+1), AND INCREMENT IP AFTER
C   EACH FIELD.
C   REVISION LOG:
C    5/19/93  HJI  ORIGINAL VERSION.

      IMPLICIT NONE
      INTEGER*4  LENTEXT, NEXTBL, NEXTNB
      INTEGER*4  IP, KI, KF, LASTNB
      REAL*4     XVAL(*)
      CHARACTER*(*) INPLINE

      LASTNB=LENTEXT(INPLINE)
      KI=1
      KF=0
   10 KI=NEXTNB(INPLINE,KI,LASTNB)
      KF=NEXTBL(INPLINE,KI,LASTNB)
      IF (KF .EQ. 0) THEN
	 KF=LASTNB
      ELSE
	 KF=KF-1
      ENDIF
      IP=IP+1
      READ(INPLINE(KI:KF),*) XVAL(IP)
      KI=KF+1
      IF (KF .LT. LASTNB) GO TO 10

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SETDEF
C***********************************************************************
C SUBROUTINE SETDEF SETS THE DEFAULT VALUES OF THE PARAMETERS
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/OBSTR/NSTR, NOBSTR
      INTEGER*4    NSTR, NOBSTR(NVAR)

      INTEGER*4 I
      TITLE=' '
      N=0
      NV=0
      IRS=0
      ICM=0
      NREP=1
      IRP=0
      IDATA=0
      IHIST=0
      ICORR=0
      ISPLAT=0
      NSTR=0
      DO 20 I=1,NVAR
      IDIST(I)=0
      NOBSTR(I)=0
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SIFT (XV,N)
C***********************************************************************
C SUBROUTINE SIFT WILL SORT A VECTOR OF DATA IN INCREASING ORDER
C***********************************************************************
      REAL*4 XV(N), A
      INTEGER*4 N, M, K, J, I, L
      M=N
   10 M=M/2
      IF (M) 30,20,30
   20 RETURN
   30 K=N-M
      J=1
   40 I=J
   50 L=I+M
      IF (XV(I)-XV(L)) 70,70,60
   60 A=XV(I)
      XV(I)=XV(L)
      XV(L)=A
      I=I-M
      IF (I) 70,70,50
   70 J=J+1
      IF (J-K) 40,40,10
      END
      SUBROUTINE SSPEV(A,N,E,V,LDV,WORK,JOB,INFO)
C***********************************************************************
C SUBROUTINE SSPEV IS USED IN THE POSITIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  SSPEV
C***DATE WRITTEN   800808   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D4A1
C***KEYWORDS  EISPACK,EIGENVALUES,EIGENVECTORS,SYMMETRIC,REAL,PACKED
C***AUTHOR  KAHANER, K. K., (NBS)
C           MOLER, C. B., (U. OF NEW MEXICO)
C           STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  TO COMPUTE THE EIGENVALUES AND, OPTIONALLY, THE EIGEN-
C            VECTORS OF A REAL SYMMETRIC MATRIX STORED IN PACKED FORM.
C***DESCRIPTION
C     LICEPACK.  THIS VERSION DATED 08/08/80.
C     DAVID KAHANER, CLEVE MOLER, PETE STEWART
C          N.B.S.       U.N.M.     N.B.S./U.MD.
C
C     ABSTRACT
C      SSPEV COMPUTES THE EIGENVALUES AND, OPTIONALLY, THE EIGENVECTORS
C      OF A REAL SYMMETRIC MATRIX STORED IN PACKED FORM.
C
C     CALL SEQUENCE PARAMETERS-
C       (THE VALUES OF PARAMETERS MARKED WITH * (STAR) WILL BE  CHANGED
C         BY SSPEV.)
C
C        A*      REAL(N*(N+1)/2)
C                REAL SYMMETRIC PACKED INPUT MATRIX.  CONTAINS UPPER
C                TRIANGLE AND DIAGONAL OF A, BY COLUMN (ELEMENTS
C                11, 12, 22, 13, 23, 33, ...).
C
C        N       INTEGER
C                SET BY THE USER TO
C                THE ORDER OF THE MATRIX A.
C
C        E*      REAL(N)
C                ON RETURN FROM SSPEV, E CONTAINS THE EIGENVALUES OF A.
C                SEE ALSO INFO BELOW.
C
C        V*      REAL(LDV,N)
C                ON RETURN FROM SSPEV, IF THE USER HAS SET JOB
C                = 0        V IS NOT REFERENCED.
C                = NONZERO  THE N EIGENVECTORS OF A ARE STORED IN THE
C                FIRST N COLUMNS OF V.  SEE ALSO INFO BELOW.
C
C        LDV     INTEGER
C                SET BY THE USER TO
C                THE LEADING DIMENSION OF THE ARRAY V IF JOB IS ALSO
C                SET NONZERO.  IN THAT CASE, N MUST BE .LE. LDV.
C                IF JOB IS SET TO ZERO, LDV IS NOT REFERENCED.
C
C        WORK*   REAL(2N)
C                TEMPORARY STORAGE VECTOR.  CONTENTS CHANGED BY SSPEV.
C
C        JOB     INTEGER
C                SET BY THE USER TO
C                = 0        EIGENVALUES ONLY TO BE CALCULATED BY SSPEV.
C                           NEITHER V NOR LDV ARE REFERENCED.
C                = NONZERO  EIGENVALUES AND VECTORS TO BE CALCULATED.
C                           IN THIS CASE, A & V MUST BE DISTINCT ARRAYS.
C                           ALSO, IF LDA .GT. LDV, SSPEV CHANGES ALL THE
C                           ELEMENTS OF A THRU COLUMN N.  IF LDA < LDV,
C                           SSPEV CHANGES ALL THE ELEMENTS OF V THROUGH
C                           COLUMN N.  IF LDA=LDV, ONLY A(I,J) AND V(I,
C                           J) FOR I,J = 1,...,N ARE CHANGED BY SSPEV.
C
C       INFO*   INTEGER
C               ON RETURN FROM SSPEV, THE VALUE OF INFO IS
C               = 0 FOR NORMAL RETURN.
C               = K IF THE EIGENVALUE ITERATION FAILS TO CONVERGE.
C                   EIGENVALUES AND VECTORS 1 THROUGH K-1 ARE CORRECT.
C
C
C     ERROR MESSAGES-
C          NO. 1   RECOVERABLE  N IS GREATER THAN LDV AND JOB IS NONZERO
C          NO. 2   RECOVERABLE  N IS LESS THAN ONE
C
C     SUBROUTINES USED
C
C      EISPACK- IMTQL2, TQLRAT, TRBAK3, TRED3
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IMTQL2,TQLRAT,TRBAK3,TRED3,XERRWV
C***END PROLOGUE  SSPEV
C***********************************************************************
      INTEGER*4 I, INFO, J, LDV, M, N, JOB, NV
      REAL*4 A((N*(N+1))/2), E(N), V(LDV,N), WORK(1)
C***FIRST EXECUTABLE STATEMENT  SSPEV
       NV=(N*(N+1))/2
       IF(N .GT. LDV) CALL XERRWV('SSPEV-N .GT. LDV.',1,1,0,0,0,0,0.,0.)
       IF(N .GT. LDV) RETURN
       IF(N .LT. 1) CALL XERRWV('SSPEV-N .LT. 1',2,1,0,0,0,0,0.,0.)
       IF(N .LT. 1) RETURN
C
C       CHECK N=1 CASE
C
      E(1) = A(1)
      INFO = 0
      IF(N .EQ. 1) RETURN
C
      IF(JOB.NE.0) GO TO 20
C
C     EIGENVALUES ONLY
C
      CALL TRED3(N,NV,A,E,WORK(1),WORK(N+1))
      CALL TQLRAT(N,E,WORK(N+1),INFO)
      RETURN
C
C     EIGENVALUES AND EIGENVECTORS
C
   20 CALL TRED3(N,NV,A,E,WORK(1),WORK(1))
      DO 30 I = 1, N
	DO 25 J = 1, N
   25     V(I,J) = 0.
   30   V(I,I) = 1.
      CALL IMTQL2(LDV,N,E,WORK,V,INFO)
      M = N
      IF(INFO .NE. 0) M = INFO - 1
      CALL TRBAK3(LDV,N,NV,A,M,V)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE STUDENT(J,IDT,ISTUD)
C***********************************************************************
C SUBROUTINE STUDENT GENERATES THE STUDENT T DISTRIBUTION
C                X
C               _
C      F(K,X) = \  1/(BETA(1/2,K/2)*(1+T**2/K)**((K+1)/2))  dT
C               -
C            - INFINITY
C***********************************************************************
      IMPLICIT NONE
C   Declarations for external functions:
      REAL*4     INVCST
C   Declarations for arguments and local variables:
      INTEGER*4  I, IDT, J, NDF
      REAL*4   AVG, PROBINC, R, STDDEV, STRTPT

C   Reset the starting point to the beginning of the next subinterval
C   unless a RANDOM SAMPLE has been specified
      LOGICAL  LOGDIST
C***********************************************************************
C THE FOLLOWING NINE LINES OF CODE ARE REQUIRED
C***********************************************************************
      INTEGER*4  NMAX, NVAR
      INTEGER*4  ISTUD, NSTUD
      PARAMETER (NSTUD=100)
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      COMMON/STU/ STUMAX, STUMIN, STUDMX, STUDMN
      REAL*4      STUMAX, STUMIN, STUDMX(NSTUD), STUDMN(NSTUD)
      REAL*4      RNUM

C   The following function definition is required
      INTEGER*4  LOC
      LOC(I,J)=(J-1)*N+I

      IF (IDT .EQ. 13) THEN
	 LOGDIST=.FALSE.
      ELSE IF (IDT .EQ. 14) THEN
	 LOGDIST=.TRUE.
      ENDIF

      READ (8) NDF, AVG, STDDEV
C
C   Set the starting point (STRTPT) equal to zero and the probability
C   increment (PROBINC) equal to 1/N for a LHS where N is the sample size
      STRTPT=0.0
      PROBINC=1.0/FLOAT(N)
C
C   If a RANDOM SAMPLE has been specified in the parameter list then the
C   argument IRS has been set equal to 1 in the main program, hence the
C   probability increment is set equal to 1 so that all observations are
C   selected by using the interval (0,1)
      IF (IRS .EQ. 1) PROBINC=1.0
C
C   This loop will obtain the N sample values:
      DO 70 I=1,N
C   R is a randomly selected point in the current subinterval obtained
C   by using the random number generator RAN
         RNUM=RAN(ISEED)
	 R=STRTPT+PROBINC*RNUM

C   This loop will select the specific value of the random variable
C   corresponding to R through the inverse cumulative distribution
C   function.  These values are stored in the vector X through the
C   use of the LOC function.
	 X(LOC(I,J))=STDDEV*INVCST(NDF,R)+AVG
         IF (X(LOC(I,J)) .LT. STUDMN(ISTUD)) X(LOC(I,J))=STUDMN(ISTUD)
         IF (X(LOC(I,J)) .GT. STUDMX(ISTUD)) X(LOC(I,J))=STUDMX(ISTUD)
	 IF (LOGDIST) X(LOC(I,J))=EXP(X(LOC(I,J)))

C   Reset the starting point to the beginning of the next subinterval
C   unless a RANDOM SAMPLE has been specified
	 IF (IRS .NE. 1) STRTPT=STRTPT+PROBINC
         IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   70 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      ISTUD=ISTUD+1
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C***********************************************************************
C SUBROUTINE TQLRAT IS USED IN THE POSITIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  TQLRAT
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EISPACK,EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH ET AL
C***PURPOSE  COMPUTES EIGENVALUES OF SYMMETRIC TRIDIAGONAL MATRIX
C            A RATIONAL VARIANT OF THE QL METHOD.
C***DESCRIPTION
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG(A,B) FOR SQRT(A**2 + B**2).
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  *MATRIX EIGENSYSTEM ROUTINES-EISPACKGUIDE*,
C                 B.T.SMITH,J.M.BOYLE,J.J.DONGARRA,B.S.GARBOW,
C                 Y.I.KEBE,V.C.KLEMA,C.B.MOLER,SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQLRAT
C
C***********************************************************************
      INTEGER*4 I,J,L,M,N,II,L1,MML,IERR
      REAL*4 D(N),E2(N)
      REAL*4 B,C,F,G,H,P,R,S,MACHEP
      REAL*4 PYTHAG
C
      DATA MACHEP/1.0E0/
C***FIRST EXECUTABLE STATEMENT  TQLRAT
      IF (MACHEP .NE. 1.0E0) GO TO 10
   05 MACHEP = 0.5E0*MACHEP
      IF (1.0E0 + MACHEP .GT. 1.0E0) GO TO 05
      MACHEP = 2.0E0*MACHEP
C
   10 IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0E0
      B = 0.0E0
      E2(N) = 0.0E0
C
      DO 290 L = 1, N
	 J = 0
	 H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
	 IF (B .GT. H) GO TO 105
	 B = H
	 C = B * B
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
	    IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
	 J = J + 1
C     .......... FORM SHIFT ..........
	 L1 = L + 1
	 S = SQRT(E2(L))
	 G = D(L)
	 P = (D(L1) - G) / (2.0E0 * S)
	 R = PYTHAG(P,1.0E0)
	 D(L) = S / (P + SIGN(R,P))
	 H = G - D(L)
C
	 DO 140 I = L1, N
  140    D(I) = D(I) - H
C
	 F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
	 G = D(M)
	 IF (G .EQ. 0.0E0) G = B
	 H = G
	 S = 0.0E0
	 MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
	 DO 200 II = 1, MML
	    I = M - II
	    P = G * H
	    R = P + E2(I)
	    E2(I+1) = S * R
	    S = E2(I) / R
	    D(I+1) = H + S * (H + D(I))
	    G = D(I) - E2(I) / G
	    IF (G .EQ. 0.0E0) G = B
	    H = G * P / R
  200    CONTINUE
C
	 E2(L) = S * G
	 D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
	 IF (H .EQ. 0.0E0) GO TO 210
	 IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
	 E2(L) = H * E2(L)
	 IF (E2(L) .NE. 0.0E0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
	 IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
	 DO 230 II = 2, L
	    I = L + 2 - II
	    IF (P .GE. D(I-1)) GO TO 270
	    D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C***********************************************************************
C SUBROUTINE TRBAK3 IS USED IN THE POSITIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  TRBAK3
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C4
C***KEYWORDS  EISPACK,EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH ET AL
C***PURPOSE  FORMS EIGENVECTORS OF REAL SYMMETRIC MATRIX FROM THE
C            EIGENVECTORS OF SYMMETRIC TRIDIAGONAL MATRIX FORMED
C            TRED3.
C***DESCRIPTION
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  *MATRIX EIGENSYSTEM ROUTINES-EISPACKGUIDE*,
C                 B.T.SMITH,J.M.BOYLE,J.J.DONGARRA,B.S.GARBOW,
C                 Y.I.KEBE,V.C.KLEMA,C.B.MOLER,SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRBAK3
C
C***********************************************************************
      INTEGER*4 I,J,K,L,M,N,IK,IZ,NM,NV
      REAL*4 A(NV),Z(NM,M)
      REAL*4 H,S
C
C***FIRST EXECUTABLE STATEMENT  TRBAK3
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
	 L = I - 1
	 IZ = (I * L) / 2
	 IK = IZ + I
	 H = A(IK)
	 IF (H .EQ. 0.0E0) GO TO 140
C
	 DO 130 J = 1, M
	    S = 0.0E0
	    IK = IZ
C
	    DO 110 K = 1, L
	       IK = IK + 1
	       S = S + A(IK) * Z(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
	    S = (S / H) / H
	    IK = IZ
C
	    DO 120 K = 1, L
	       IK = IK + 1
	       Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C***********************************************************************
C SUBROUTINE TRED3 IS USED IN THE POSITIVE DEFINITE CHECK OF THE
C CORRELATION MATRIX
C
C***BEGIN PROLOGUE  TRED3
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B1
C***KEYWORDS  EISPACK,EIGENVALUES,EIGENVECTORS
C***AUTHOR  SMITH ET AL
C***PURPOSE  REDUCE REAL SYMMETRIC MATRIX STORED IN PACKED FORM TO
C            SYMMETRIC TRIDIAGONAL MATRIX USING ORTHOGONAL
C            TRANSFORMATIONS.
C***DESCRIPTION
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT.
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  *MATRIX EIGENSYSTEM ROUTINES-EISPACKGUIDE*,
C                 B.T.SMITH,J.M.BOYLE,J.J.DONGARRA,B.S.GARBOW,
C                 Y.I.KEBE,V.C.KLEMA,C.B.MOLER,SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRED3
C
C***********************************************************************
      INTEGER*4 I,J,K,L,N,II,IZ,JK,NV
      REAL*4 A(NV),D(N),E(N),E2(N)
      REAL*4 F,G,H,HH,SCALE
C
C***FIRST EXECUTABLE STATEMENT  TRED3
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO  300 II = 1, N
	 I = N + 1 - II
	 L = I - 1
	 IZ = (I * L) / 2
	 H = 0.0E0
	 SCALE = 0.0E0
	 IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
	 DO 120 K = 1, L
	    IZ = IZ + 1
	    D(K) = A(IZ)
	    SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
	 IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = 0.0E0
	 E2(I) = 0.0E0
	 GO TO 290
C
  140    DO 150 K = 1, L
	    D(K) = D(K) / SCALE
	    H = H + D(K) * D(K)
  150    CONTINUE
C
	 E2(I) = SCALE * SCALE * H
	 F = D(L)
	 G = -SIGN(SQRT(H),F)
	 E(I) = SCALE * G
	 H = H - F * G
	 D(L) = F - G
	 A(IZ) = SCALE * D(L)
	 IF (L .EQ. 1) GO TO 290
	 F = 0.0E0
C
	 DO 240 J = 1, L
	    G = 0.0E0
	    JK = (J * (J-1)) / 2
C     .......... FORM ELEMENT OF A*U ..........
	    DO 180 K = 1, L
	       JK = JK + 1
	       IF (K .GT. J) JK = JK + K - 2
	       G = G + A(JK) * D(K)
  180       CONTINUE
C     .......... FORM ELEMENT OF P ..........
	    E(J) = G / H
	    F = F + E(J) * D(J)
  240    CONTINUE
C
	 HH = F / (H + H)
	 JK = 0
C     .......... FORM REDUCED A ..........
	 DO 260 J = 1, L
	    F = D(J)
	    G = E(J) - HH * F
	    E(J) = G
C
	    DO 260 K = 1, J
	       JK = JK + 1
	       A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
C
  290    D(I) = A(IZ+1)
	 A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRIANG(J)
C***********************************************************************
C SUBROUTINE TRIANG GENERATES THE TRIANGULAR DISTRIBUTION
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 PROBINC, A, B, C, C1, C2, STRTPT, R
      REAL*4 RNUM
      INTEGER*4 J, LOC, I
      LOC(I,J)=(J-1)*N+I
      PROBINC=1./FLOAT(N)
      IF(IRS.EQ.1)PROBINC=1.0
      READ(8)A,B,C
      C1=C-A
      C2=(B-A)/C1
      STRTPT=0.
      DO 10 I=1,N
        RNUM=RAN(ISEED)
	R=PROBINC*RNUM+STRTPT
	IF(R.LE.C2)THEN
	  X(LOC(I,J))=A+SQRT(R*C1*(B-A))
	ELSE
	  X(LOC(I,J))=C-SQRT((1.-R)*C1*(C-B))
	ENDIF
	IF(IRS.EQ.0)STRTPT=STRTPT+PROBINC
        IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   10 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE UNIFRM(J,IDT)
C***********************************************************************
C SUBROUTINE UNIFRM GENERATES THE UNIFORM, LOGUNIFORM, UNIFORM*,
C OR THE LOGUNIFORM* DISTRIBUTIONS
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)

      REAL*4 A, B, PROBINC, STRTPT
      REAL*4 PROB, RN, STRTY, YPROB, YPROBT
      REAL*4 RNUM
      INTEGER*4 J, IDT, LOC, I, NINT, NPT, K, ID
      LOC(I,J)=(J-1)*N+I
      I=0
      PROB=0.0
      RN=0.0
      STRTY=0.0
      YPROB=0.0
      YPROBT=0.0
      IF(IDT.EQ.6.OR.IDT.EQ.7)THEN
	READ(8)NINT
      ELSE
	NINT=1
	NPT=N
      ENDIF
      PROB=1./FLOAT(N)
      DO 40 K=1,NINT
	IF(IDT.EQ.4.OR.IDT.EQ.5)THEN
	   READ(8)A,B
	ELSE
	   READ(8)NPT,A,B
	   IF(NPT.EQ.0)GO TO 40
	ENDIF
	IF(IDT.EQ.5.OR.IDT.EQ.7)THEN
	  A=LOG10(A)
	  B=LOG10(B)
	ENDIF
	PROBINC=(B-A)/FLOAT(NPT)
	IF(IRS.NE.0)PROBINC=B-A
	DO 30 ID=1,NPT
	  I=I+1
	  IF(IRS.EQ.0)THEN
	    STRTPT=(ID-1)*PROBINC
	    STRTY=(ID-1)*PROB
	  ELSE
	    STRTPT=0.0
	  ENDIF
          RNUM=RAN(ISEED)
	  X(LOC(I,J))=A+STRTPT+PROBINC*RNUM
          YPROB=YPROBT+STRTY+PROB*RNUM
	  IF(IDT.EQ.5.OR.IDT.EQ.7)X(LOC(I,J))=10.**(X(LOC(I,J)))
          IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), YPROB, RNUM
   30   CONTINUE
      YPROBT=YPROB
   40 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
      SUBROUTINE VIF
C***********************************************************************
C SUBROUTINE VIF COMPUTES THE VARIANCE INFLATION FACTOR OF A
C CORRELATION MATRIX
C***********************************************************************
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/CMATR/CORR, LCM, NCM
      REAL*4       CORR((NVAR*(NVAR+1))/2)
      INTEGER*4    LCM(NVAR), NCM

      REAL*4 CRKMX, CRK
      INTEGER*4 LOC1, I, J
      LOC1(I,J)=J+(I*I-I)/2
      CALL DSINV(NV)
      CRKMX=0.0
      DO 630 I=1,NV
	CRK=CORR(LOC1(I,I))
	IF(CRK.GT.CRKMX)CRKMX=CRK
  630 CONTINUE
      WRITE(6,9001)CRKMX
      RETURN
 9001 FORMAT('0THE VARIANCE INFLATION FACTOR FOR THIS MATRIX IS',
     1       F6.2)
      END
      SUBROUTINE WRTCRD(ITYPE,LABEL)
C***********************************************************************
C SUBROUTINE WRTCRD PROCESSES THE DISTRIBUTION PARAMETER STATEMENTS
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (LENC=80)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP

      CHARACTER LABEL*(*), LABVAR*(LENC), BLANK*1
      PARAMETER (BLANK=' ')
      INTEGER*1 ITYPE, NC, IC
      NV=NV+1
      IDIST(NV)=ITYPE
      NC=LEN(LABEL)
      IC=0
      LABVAR=BLANK
  100 CONTINUE
      IC=IC+1
      IF(IC.GT.NC)GO TO 200
      IF(LABEL(IC:IC).EQ.BLANK)GO TO 100
      LABVAR=LABEL(IC:NC)
  200 CONTINUE
      WRITE(9)LABVAR
      RETURN
      END
      SUBROUTINE WRTPAR
C***********************************************************************
C SUBROUTINE WRTPAR PRINTS OUT THE DISTRIBUTION PARAMETERS
C***********************************************************************
      PARAMETER (NVAR=100)
      PARAMETER (LENC=60)
      PARAMETER (LEND=14)
      COMMON/HEADNG/TITLE
      CHARACTER*100 TITLE
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP

      REAL*4 A, B, P, Q, XMN, VAR, C, D
      INTEGER*4 I, ID, NINT, J, INT
      CHARACTER LABEL*(LENC)
      CHARACTER*15 DIST(LEND)
      DATA DIST(1)/'BETA'/,DIST(2)/'NORMAL'/,DIST(3)/'LOGNORMAL'/,
     1     DIST(4)/'UNIFORM'/,DIST(5)/'LOGUNIFORM'/,
     2     DIST(6)/'UNIFORM*'/,DIST(7)/'LOGUNIFORM*'/,
     3     DIST(8)/'TRIANGULAR'/,DIST(9)/'USER SUPPLIED'/,
     4     DIST(10)/'EXPONENTIAL'/, DIST(11)/'RAYLEIGH '/,
     5     DIST(12)/'RAYL-EXPON '/, DIST(13)/'STUDENT '/,
     6     DIST(14)/'LOGSTUDENT '/
      REWIND 8
      REWIND 9
      WRITE(6,9001) TITLE
      DO 100 I=1,NV
	ID=IDIST(I)
	READ(9)LABEL
	IF(ID.EQ.1)THEN
	  READ(8)A,B,P,Q
	  XMN=(A*Q+B*P)/(P+Q)
	  VAR=P*Q*(B-A)**2/((P+Q+1)*(P+Q)**2)
	  WRITE(6,9002)I,DIST(ID),A,B,LABEL
	  WRITE(6,9006)P,Q,XMN,VAR
	  GO TO 100
	ELSE IF(ID.EQ.8 .OR. ID .EQ. 10 .OR. ID .EQ. 11) THEN
	  READ(8)A,B,C
	  WRITE(6,9007)I,DIST(ID),LABEL,A,B,C
	  GO TO 100
	ELSE IF(ID.EQ.9)THEN
	  WRITE(6,9003)I,LABEL
	  GO TO 100
	ELSE IF(ID.EQ.6.OR.ID.EQ.7)THEN
	  READ(8)NINT
	  WRITE(6,9004)I,DIST(ID),NINT,LABEL
	  DO 50 J=1,NINT
	    READ(8)INT,A,B
	    WRITE(6,9005)INT,A,B
   50     CONTINUE
	  GO TO 100
	ELSE IF (ID .EQ. 12) THEN
	  READ (8) A,B,C,D
	  WRITE(6,9008) I,DIST(ID),LABEL,A,B,C,D
	  GO TO 100
	ELSE IF (ID .EQ. 13) THEN
	    READ (8) INT,A,B
	    WRITE(6,9013) I,DIST(ID),LABEL, INT,A,B
	    GO TO 100
	ELSE IF (ID .EQ. 14) THEN
	    READ (8) INT,A,B
	    WRITE(6,9014) I,DIST(ID),LABEL, INT,A,B
	    GO TO 100
	ELSE
	  READ(8)A,B
	  WRITE(6,9002)I,DIST(ID),A,B,LABEL
	ENDIF
  100 CONTINUE
      RETURN
 9001 FORMAT('1'//4X,A100//,4X,'VARIABLE  DISTRIBUTION',10X,
     1       'RANGE',12X,'LABEL')
 9002 FORMAT('0',5X,I3,5X,A,2X,1PG12.4,' TO ',1PG12.4,2X,A)
 9003 FORMAT('0',5X,I3,5X,'USER SUPPLIED DISTRIBUTION ',12X,A)
 9004 FORMAT('0',5X,I3,5X,A,2X,'WITH ',I2,' SUBINTERVALS',6X,A)
 9005 FORMAT(' ',13X,I4,' OBS',5X,1PG12.4,' TO ',1PG12.4)
 9006 FORMAT(' ',13X,'WITH PARAMETERS  P = ',F12.4,2X,'Q = ',F12.4/,
     1       14X,'THIS CHOICE OF PARAMETERS GIVES A ',/,14X,
     2       'POPULATION MEAN OF ',1PG12.4,'  AND A',/,14X,
     3       'POPULATION VARIANCE OF ',1PG12.4)
 9007 FORMAT('0   ',I3,5X,A,'  WITH PARAMETERS BELOW   ',A/27X,
     1       'A= ',1PG15.4,/,27X,'B= ',1PG15.4,/,27X,'C= ',1PG15.4)
 9008 FORMAT('0   ',I3,5X,A,'  WITH PARAMETERS BELOW   ',A/
     1       27X,'A= ',1PG15.4/27X,'B= ',1PG15.4/27X,'C= ',1PG15.4/
     2       27X,'D= ',1PG15.4)
 9012 FORMAT('0',5X,I3,5X,A,2X,'WITH PARAMETERS BELOW',5X,A,/
     1        27X,'A= ',1PG15.4/27X,'B= ',1PG15.4/
     2        27X,'C= ',1PG15.4/22X,'LAMBDA= ',1PG15.4)
 9013 FORMAT('0',5X,I3,5X,A,2X,'WITH PARAMETERS BELOW',5X,A,/
     1        27X,I3,' DEGREES OF FREEDOM'/27X,'AVG=',1PG15.4/
     &        27X,'STDDEV= ',1PG12.4)
 9014 FORMAT('0',5X,I3,5X,A,2X,'WITH PARAMETERS BELOW',5X,A,/
     1        27X,I3,' DEGREES OF FREEDOM'/27X,'LOG AVG=',1PG15.4/
     &        27X,'LOG STDDEV= ',1PG12.4)
      END
      SUBROUTINE XERPRT(MESSG)
C***********************************************************************
C SUBROUTINE XERPRT IS AN ERROR HANDLING ROUTINE USED IN THE
C POSITIVE DEFINITE CHECK OF THE CORRELATION MATRIX
C
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  PRINTS ERROR MESSAGES
C***DESCRIPTION
C     ABSTRACT
C        PRINT THE CHARACTER MESSAGE IN MESSG
C        ON EACH FILE INDICATED BY XGETUA.
C***REFERENCES  JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR
C                 HANDLING PACKAGE*, SAND78-1189, SANDIA LABORATORIES,
C                 1978.
C***ROUTINES CALLED  I1MACH,XGETUA
C***END PROLOGUE  XERPRT
C***********************************************************************
      INTEGER*4 LUN(5), I1MACH, FIRST, LAST
      INTEGER*4 KUNIT, IUNIT, NUNIT, I
      CHARACTER F*10, G*14, LA*1, LCOM*1, LBLANK*1, MESSG*(*)
C***FIRST EXECUTABLE STATEMENT  XERPRT
C     PRINT THE MESSAGE
      CALL XGETUA(LUN,NUNIT)
      DO 50 KUNIT = 1,NUNIT
	 IUNIT = LUN(KUNIT)
	 IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         FIRST = 1
         DO WHILE (FIRST .LE. LEN(MESSG))
           LAST = FIRST + 71
           IF (LAST .GT. LEN(MESSG)) LAST = LEN(MESSG)
	   WRITE (IUNIT,'(1X,A)') MESSG(FIRST:LAST)
           FIRST = LAST + 1
         END DO
50    CONTINUE
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***********************************************************************
C SUBROUTINE XERRWV IS AN ERROR HANDLING ROUTINE USED IN THE
C POSITIVE DEFINITE CHECK OF THE CORRELATION MATRIX
C
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  (NONE)
C***AUTHOR  JONES, R. E., (SLA)
C***PURPOSE  PROCESSES ERROR MESSAGE ALLOWING 2 INTEGER AND TWO REAL
C            VALUES TO BE INCLUDED IN THE MESSAGE.
C***DESCRIPTION
C     ABSTRACT
C        XERRWV PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER
C        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE
C        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL.
C        (SEE SUBROUTINE XSETF FOR DETAILS.)
C        IN ADDITION, UP TO TWO INTEGER VALUES AND TWO REAL
C        VALUES MAY BE PRINTED ALONG WITH THE MESSAGE.
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        MESSG - THE CHARACTER MESSAGE TO BE PROCESSED.
C        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE.
C                NERR MUST NOT BE ZERO.
C        LEVEL - ERROR CATEGORY.
C                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR.
C                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS
C                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.)
C                =0 MEANS THIS IS A WARNING MESSAGE ONLY.
C                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE
C                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY
C                   TIMES THIS CALL IS EXECUTED.
C        NI    - NUMBER OF INTEGER VALUES TO BE PRINTED. (0 TO 2)
C        I1    - FIRST INTEGER VALUE.
C        I2    - SECOND INTEGER VALUE.
C        NR    - NUMBER OF REAL VALUES TO BE PRINTED. (0 TO 2)
C        R1    - FIRST REAL VALUE.
C        R2    - SECOND REAL VALUE.
C
C     EXAMPLE
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM'
C    1            // ' (R2).',77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  19 MAR 1980
C***REFERENCES  JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR
C                 HANDLING PACKAGE*, SAND78-1189, SANDIA LABORATORIES,
C                 1978.
C***ROUTINES CALLED  I1MACH,J4SAVE,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
C***********************************************************************
      REAL*4 R1, R2
      INTEGER*4 NERR, LEVEL, NI, I1, I2
      INTEGER*4 NR, LUN(5), LKNTRL, J4SAVE, MAXMES
      INTEGER*4 KOUNT, LERR, LLEVEL, MKNTRL, NUNIT, KUNIT
      INTEGER*4 IUNIT, I1MACH, IFATAL
      CHARACTER MESSG*(*)
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((LEN(MESSG).GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
	 IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...')
	 CALL XERPRT('XERRWV -- INVALID INPUT')
	 IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.')
	 IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,IUNIT)
	 RETURN
   10 CONTINUE
C     RECORD MESSAGE
      LLEVEL = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,LEN(MESSG),NERR,LEVEL,KOUNT)
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
	 IF (LKNTRL.LE.0) GO TO 20
	    CALL XERPRT(' ')
C           INTRODUCTION
	    IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1  ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.')
	    IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...')
	    IF (LLEVEL.EQ.1) CALL XERPRT ('RECOVERABLE ERROR IN...')
	    IF (LLEVEL.EQ.2) CALL XERPRT ('FATAL ERROR IN...')
   20    CONTINUE
C        MESSAGE
	 CALL XERPRT(MESSG)
	 CALL XGETUA(LUN,NUNIT)
	 DO 50 KUNIT=1,NUNIT
	    IUNIT = LUN(KUNIT)
	    IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
	    IF (NI.GE.1) WRITE (IUNIT,22) I1
	    IF (NI.GE.2) WRITE (IUNIT,23) I2
	    IF (NR.GE.1) WRITE (IUNIT,24) R1
	    IF (NR.GE.2) WRITE (IUNIT,25) R2
   22       FORMAT (11X,'IN ABOVE MESSAGE, I1=',I10)
   23       FORMAT (11X,'IN ABOVE MESSAGE, I2=',I10)
   24       FORMAT (11X,'IN ABOVE MESSAGE, R1=',E20.10)
   25       FORMAT (11X,'IN ABOVE MESSAGE, R2=',E20.10)
	    IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
	       WRITE (IUNIT,30) LERR
   30          FORMAT (' ERROR NUMBER =',I10)
   40       CONTINUE
   50    CONTINUE
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
	 IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.')
	 IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.')
C        PRINT ERROR SUMMARY
	 CALL XERSAV(' ',-1,0,0,IUNIT)
  120 CONTINUE
C     ABORT
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***********************************************************************
C SUBROUTINE XERSAV IS AN ERROR HANDLING ROUTINE USED IN THE
C POSITIVE DEFINITE CHECK OF THE CORRELATION MATRIX
C
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  (NONE)
C***AUTHOR  JONES, R. E., (SLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NERR, LEVEL are as in XERRWV,
C       when NMESSG=0 the tables will be dumped and cleared,
C       when NMESSG < zero the tables will be dumped but not cleared,
C       and when NMESSG > zero the tables are not dumped.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR
C                 HANDLING PACKAGE*, SAND78-1189, SANDIA LABORATORIES,
C                1978.
C***ROUTINES CALLED  I1MACH,XGETUA
C***END PROLOGUE  XERSAV
C***********************************************************************
      INTEGER*4 NMESSG, NERR, LEVEL, ICOUNT, LASTCHR
      INTEGER*4 NERTAB(10), LEVTAB(10), KOUNT(10), KOUNTX, NCHAR
      INTEGER*4 NUNIT, KUNIT, IUNIT, I, II, LUN(5), NCOL, I1MACH
      CHARACTER MESSG*(*), F*13, MESTAB(10)*10
      DATA MESTAB / 10*' ' /
      DATA NERTAB / 10*0 /
      DATA LEVTAB / 10*0 /
      DATA KOUNT  / 10*0 /
      DATA KOUNTX/0/
      DATA F / '(1X,A10,3I10)' /
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
	 IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
	 CALL XGETUA(LUN,NUNIT)
	 DO 60 KUNIT=1,NUNIT
	    IUNIT = LUN(KUNIT)
	    IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
	    WRITE (IUNIT,10)
   10       FORMAT ('0          ERROR MESSAGE SUMMARY' /
     1              ' START OF MSG    NERR     LEVEL     COUNT')
C           PRINT BODY OF TABLE
	    DO 20 I=1,10
	       IF (KOUNT(I).EQ.0) GO TO 30
	       WRITE (IUNIT,F) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
	    IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED=',I10)
	    WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
	 IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
	 DO 70 I=1,10
   70       KOUNT(I) = 0
	 KOUNTX = 0
	 RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      LASTCHR = 10
      IF (LEN(MESSG) .LT. 10) LASTCHR = LEN(MESSG)
      DO 90 I=1,10
	 II = I
	 IF (KOUNT(I).EQ.0) GO TO 110
	 IF (MESSG(1:LASTCHR).NE.MESTAB(I)) GO TO 90
	 IF (NERR.NE.NERTAB(I)) GO TO 90
	 IF (LEVEL.NE.LEVTAB(I)) GO TO 90
	 GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
	 KOUNTX = KOUNTX+1
	 ICOUNT = 1
	 RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
	 ICOUNT = KOUNT(II)
	 RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MESSG(1:LASTCHR)
	 NERTAB(II) = NERR
	 LEVTAB(II) = LEVEL
	 KOUNT(II)  = 1
	 ICOUNT = 1
	 RETURN
      END
      SUBROUTINE XGETUA(IUNIT,N)
C***********************************************************************
C SUBROUTINE XGETUA IS AN ERROR HANDLING ROUTINE USED IN THE
C POSITIVE DEFINITE CHECK OF THE CORRELATION MATRIX
C
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  (NONE)
C***AUTHOR  JONES, R. E., (SLA)
C***PURPOSE  RETURNS UNIT NUMBER(S) TO WHICH ERROR MESSAGES ARE BEING
C            SENT
C***DESCRIPTION
C     ABSTRACT
C        XGETUA MAY BE CALLED TO DETERMINE THE UNIT NUMBER OR NUMBERS
C        TO WHICH ERROR MESSAGES ARE BEING SENT.
C        THESE UNIT NUMBERS MAY HAVE BEEN SET BY A CALL TO XSETUN,
C        OR A CALL TO XSETUA, OR MAY BE A DEFAULT VALUE.
C
C     DESCRIPTION OF PARAMETERS
C      --OUTPUT--
C        IUNIT - AN ARRAY OF ONE TO FIVE UNIT NUMBERS, DEPENDING
C                ON THE VALUE OF N.  A VALUE OF ZERO REFERS TO THE
C                DEFAULT UNIT, AS DEFINED BY THE I1MACH MACHINE
C                CONSTANT ROUTINE.  ONLY IUNIT(1),...,IUNIT(N) ARE
C                DEFINED BY XGETUA.  THE VALUES OF IUNIT(N+1),...,
C                IUNIT(5) ARE NOT DEFINED (FOR N .LT. 5) OR ALTERED
C                IN ANY WAY BY XGETUA.
C        N     - THE NUMBER OF UNITS TO WHICH COPIES OF THE
C                ERROR MESSAGES ARE BEING SENT.  N WILL BE IN THE
C                RANGE FROM 1 TO 5.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C***REFERENCES  JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR
C                 HANDLING PACKAGE*, SAND78-1189, SANDIA LABORATORIES,
C                 1978.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
C***********************************************************************
      INTEGER*4 IUNIT(5), N, J4SAVE, I, IND
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
	 IND = I+4
	 IF (I.EQ.1) IND = 3
	 IUNIT(I) = J4SAVE(IND,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE USRDST(J)
      INTEGER*4   I, J, K, LOC, NP
      CHARACTER*80  INPLINE
      LOGICAL   CONTIN, SPECPRB
C***********************************************************************
C THE FOLLOWING SIX LINES OF CODE ARE REQUIRED BY USRDST
C***********************************************************************
      PARAMETER (NMAX=10000)
      PARAMETER (NVAR=100)
      COMMON/PARAM/ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST,IRP
      INTEGER*4    ISEED,N,NV,IRS,ICM,NREP,IDATA,IHIST,
     1             ICORR,ISPLAT,IDIST(NVAR),IRP
      COMMON/SAMP/X
      REAL*4      X(NMAX*NVAR)
      REAL*4 XVAL, FREQ, CDF, STRTPT, PROBINC, R
      REAL*4 RNUM
C
C THE FOLLOWING LINE OF CODE IS SUPPLIED BY THE USER
C XVAL AND FREQ MUST BE DIMENSIONED TO THE NUMBER OF UNIQUE VALUES
C OR NUMBER OF SAMPLE VALUES THAT THE DISCRETE RANDOM VARIABLE TAKES
C ON AND CDF MUST BE DIMENSIONED TO THE NUMBER OF VALUES PLUS 1
C
      DIMENSION XVAL(5000),FREQ(5000),CDF(5001)
C
C
C THE FOLLOWING FUNCTION DEFINITION IS REQUIRED BY USRDST
C
      LOC(I,J)=(J-1)*N+I
C
C READ IN THE NUMBER OF DATA POINTS (NP),
C THE PROBABILITY STATUS: EQUAL OR SPECIFIED(DEFAULT),
C AND THE DISTRIBUTION TYPE: CONTINUOUS OR DISCRETE(DEFAULT).
C
C IF THE KEYWORD 'SPEC' IS INCLUDED OR IF 'EQUAL' IS OMITTED,
C NP PAIRS OF VALUE,PROBABILITY WILL BE READ.
C IF THE KEYWORD 'EQUAL' IS INCLUDED, NP VALUES WILL BE READ.
C
C IF THE KEYWORD 'DISCR' IS INCLUDED OR IF 'CONTIN' IS OMITTED,
C SAMPLING IS FROM ACTUAL VALUES ONLY.
C IF THE KEYWORD 'CONTIN' IS INCLUDED, SAMPLING USES LINEAR
C INTERPOLATION BETWEEN THE DATA.
C
C IF KEYWORDS ARE ABSENT FROM THE FIRST LINE, ASSUME VALUE,PROBABILITY
C PAIRS AND A DISCRETE DISTRIBUTION.  THIS ALLOWS COMPATIBILITY WITH
C SOME EARLIER VERSIONS OF LHS.
C
C NP IS THE NUMBER OF UNIQUE VALUES OF THE RANDOM VARIABLE
      READ (7,'(A)'), INPLINE
      READ (INPLINE,*) NP
      SPECPRB=.TRUE.
      CONTIN=.FALSE.
      IF (INDEX(INPLINE,'EQUAL') .GT. 0 .OR.
     +    INDEX(INPLINE,'equal') .GT. 0) SPECPRB=.FALSE.
      IF (INDEX(INPLINE,'CONTIN') .GT. 0 .OR.
     +    INDEX(INPLINE,'contin') .GT. 0) CONTIN=.TRUE.
C
C XVAL(K) IS THE KTH UNIQUE VALUE OF THE RANDOM VARIABLE
C FREQ(K) IS THE PROBABILITY ASSOCIATED WITH THE KTH UNIQUE VALUE
C (OR THE INTERVAL BETWEEN XVAL(K) AND XVAL(K+1)).
C NOTE THAT THE READ STATEMENT MUST BE OF THE FORM READ(7,*)....
C
      IF (SPECPRB) THEN
	 DO 20 K=1,NP
   20    READ (7,*) XVAL(K),FREQ(K)
      ELSE
C
C READ IN THE SAMPLE VALUES
C NOTE THAT THE READ STATEMENT MUST BE OF THE FORM READ(7,*)....
C
	 READ(7,*) (XVAL(K),K=1,NP)
	 IF (CONTIN) THEN
	    DO 30 K=1,NP
   30       FREQ(K)=1./(NP-1)
	 ELSE
	    DO 40 K=1,NP
   40       FREQ(K)=1./NP
	 ENDIF
      ENDIF
C
C CONSTRUCT THE CUMULATIVE DISTRIBUTION FUNCTION
C
      CDF(1)=0.0
      DO 50 K=1,NP
   50 CDF(K+1)=CDF(K)+FREQ(K)
C
C SET THE STARTING POINT (STRTPT) EQUAL TO ZERO AND THE PROBABILITY
C INCREMENT (PROBINC) EQUAL TO 1/N FOR A LHS WHERE N IS THE SAMPLE SIZE
C
      STRTPT=0.0
      PROBINC=1.0/FLOAT(N)
C
C IF A RANDOM SAMPLE HAS BEEN SPECIFIED IN THE PARAMETER LIST THEN THE
C ARGUMENT IRS HAS BEEN SET EQUAL TO 1 IN THE MAIN PROGRAM, HENCE THE
C PROBABILITY INCREMENT IS SET EQUAL TO 1 SO THAT ALL OBSERVATIONS ARE
C SELECTED BY USING THE INTERVAL (0,1)
C
      IF (IRS.EQ.1) PROBINC=1.0
C
C THIS LOOP WILL OBTAIN THE N SAMPLE VALUES
C
      DO 70 I=1,N
C
C R IS A RANDOMLY SELECTED POINT IN THE CURRENT SUBINTERVAL OBTAINED
C BY USING THE RANDOM NUMBER GENERATOR RAN
C
        RNUM=RAN(ISEED)
	R=STRTPT+PROBINC*RNUM
C
C THIS LOOP WILL SELECT THE SPECIFIC VALUE OF THE RANDOM VARIABLE
C CORRESPONDING TO R THROUGH THE INVERSE CUMULATIVE DISTRIBUTION
C FUNCTION THESE VALUES ARE STORED IN THE VECTOR X THROUGH THE
C USE OF THE LOC FUNCTION
C
	DO 60 K=1,NP
	IF (R .GE. CDF(K) .AND. R .LT. CDF(K+1)) THEN
	   IF (CONTIN) THEN
	      X(LOC(I,J))=XVAL(K)+(R-CDF(K))*(XVAL(K+1)-XVAL(K))/
     +                        FREQ(K)
	   ELSE
	      X(LOC(I,J))=XVAL(K)
	   ENDIF
	ENDIF
   60   CONTINUE
C
C RESET THE STARTING POINT TO THE BEGINNING OF THE NEXT SUBINTERVAL
C UNLESS A RANDOM SAMPLE HAS BEEN SPECIFIED
C
      IF (IRS .NE. 1) STRTPT=STRTPT+PROBINC
      IF(ISPLAT.NE.0) WRITE(20,*) J, I, X(LOC(I,J)), R, RNUM
   70 CONTINUE
      IF(ISPLAT.NE.0) WRITE(20,*)'*break'
      RETURN
      END
C-----------------------------------------------------------------------
C CMS REPLACEMENT HISTORY, Element LHS2_LHS.FOR
C *6     8-AUG-1997 10:53:21 LNSMITH "MODIFICATIONS MADE FOR USE IN INEEL PA."
C *5     6-MAR-1996 13:06:03 LNSMITH "MODIFICATIONS MADE TO STUDENT DISTRIBUTIONS."
C *4     1-FEB-1996 15:19:42 LNSMITH "UPDATED CODE VERSION NUMBER TO REFLECT RELINK REQUIREMENT WITH MODIFIED LIBRARY ROUTINES."
C *3    11-SEP-1995 15:19:54 LNSMITH "UPDATED CODE VERSION NUMBER."
C *2    24-AUG-1995 13:52:53 CMWILLI "PER EJDECOT, REMOVE UNNECESSARY MODULES, FIX CODDING ERRORS, AND REMOVE USE OF HOLLERITH
CCONSTANTS AND HOLLERITH DATA IN INTEGER VARIABLES."
C *1     9-AUG-1995 09:53:10 LNSMITH "INITIAL LOAD."
C CMS REPLACEMENT HISTORY, Element LHS2_LHS.FOR
