      SUBROUTINE SNDPST_FV3S(LM,NSTAT,LENGTH)
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C MAIN PROGRAM: ETA_SNDP
C   PRGMMR: ROGERS           ORG: NP22        DATE: 1999-09-24
C
C ABSTRACT:  THIS ROUTINE POSTS PROFILE DATA AND WRITES
C   OUTPUT IN BUFR FORMAT.  THIS REPLACES CODE THAT USED TO
C   RUN INSIDE OF CHKOUT IN THE ETA MODEL.
C     
C PROGRAM HISTORY LOG:
C   95-07-26  MIKE BALDWIN
C   96-05-07  MIKE BALDWIN - USE SMASK TO SET SOIL VARS TO MISSING
C   96-11-22  MIKE BALDWIN - ADD OPTION OF DOING MONOLITHIC FILE OR
C                            SPLIT OUT FILES OR BOTH
C   97-12-16  MIKE BALDWIN - NEW MULTI-LEVEL PARAMETERS (SUCH
C                            AS SOIL MOISTURE)
C   98-07-23  ERIC ROGERS  - MADE Y2K COMPLIANT
C   98-09-29  MIKE BALDWIN - SET ACC/AVE VARS TO MISSING AT T=0
C   99-04-01  GEOFF MANIKIN - MAJOR CHANGES FEATURING A REMOVAL OF
C                            THE DISTNICTION OF CLASS0 - ALL OUTPUT
C                            IS NOW CLASS1.  ALSO, THE FIELDS OF
C                            VISIBILITY AND CLOUD BASE PRESSURE 
C                            ARE ADDED.
C   99-09-03  JIM TUCCILLO - REDUCED MEMORY REQUIREMENTS BY CHANGING
C                            THE SIZE OF PRODAT AND INTRODUCING LPRO.
C                            ALSO, PRODAT'S STRUCTURE WAS CHANGED TO
C                            PROVIDE STRIDE-1 ACCESS.
C                            NOTE: THIS CODE CAN STILL BE MODIFIED TO
C                            REDUCE MEMORY. THE CHANGES TODAY WILL
C                            NOT AFFECT ITS FUNCTIONALITY BUT WILL
C                            ALLOW IT TO RUN ON A WINTERHAWK NODE
C                            WITHOUT PAGING. THE MEMORY REQUIREMENT
C                            IS NOW AT ABOUT 260 MBs.
C   00-03-10  GEOFF MANIKIN - CODE CHANGED TO BE READY FOR ETA EXTENSION
C                            TO 60 HOURS
C   00-5-15  ERIC ROGERS   - PUT NSOIL AND LM1 IN INCLUDED PARAMETER 
C                             FILE
C  
C   03-8-01  BINBIN ZHOU   -MODIFY TO ALSO WORK ON RSM 
C                           1.PUT NSTP and LCL1SL1 in INCLUDE file       
C                           2.USE SOIL LEVELS(4 for ETA and 2 for RSM TO 
C                             COMPUTE INDEX AND SOME PARAMETERS          
C  06-07-31  BINBIN ZHOU    -MODIFY TO WORK ON NCAR WRF
C                           The PROF output of NCAR-WRF(ARW) has no rain profile 
C                           while NCEP-WRF(NMM) has, so SNDPST.f for ARW is very different 
C                           from SNDPST for NMM.    
C  08-07-01  BINBIN ZHOU    -MODIFY TO WORK ON NCAR WRF  
C
C USAGE:   
C   INPUT ARGUMENT LIST:
C     LM - number of model layers
C     NSTAT - upper limit on number of stations to process
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       CALWXT
C       CALHEL
C       BFRIZE
C     LIBRARY:
C       BUFRLIB     
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C***
C***   7/26/95 M. BALDWIN
C*** 
C***     SOUNDING POST PROCESSOR
C***     PROGRAM TO READ THE PROFILE OUTPUT FILE ON THE CRAY
C***     AND PRODUCE DIAGNOSTIC QUANTITIES AND PACK INTO BUFR
C***
C--------------------------------------------------------------------
C
C    PARMS FOR HOURLY PROFILER OUTPUT
C      LM    - MAX NUMBER OF VERTICAL LEVELS
C      NPNT  - MAX NUMBER OF OUTPUT TIMES    
C      
C          TO HOLD ALL VARIABLES FOR ANY CLASS OF OUTPUT
C          (MAX NO MULTI-LAYER VARIABLES*LM + NO OF SINGLE LAYER VARS)
C      LCL1ML - MAX NUMBER OF MULTI-LAYER VARIABLES FOR CLASS 0 OR 1
C      LCL1SL - MAX NUMBER OF SINGLE LAYER VARIABLES FOR CLASS 0 OR 1
C      LCL1SOIL - MAX NUMBER OF SOIL LAYER VARIABLES FOR CLASS 0 OR 1
C      NSTAT - NUMBER OF STATIONS
C
C        NOTE: THESE NUMBERS WILL BE LARGER THAN THE NUMBERS
C              COMING OUT OF THE MODEL IN THE BINARY FILE SINCE
C              WE ARE COMPUTING SOME ADDITIONAL VARIABLES TO GO INTO
C              BUFR IN THIS PROGRAM.
C
C--------------------------------------------------------------------

      PARAMETER (NSOIL=4)
      PARAMETER (NSTP=88)
      PARAMETER (INCR=1)

      PARAMETER (NPNT=NSTP                
     &, SPVAL=-99999.0,SMISS=1.E10                
     &, LCL1ML=15,LCL1SL=52,LCL1SOIL=2  
     &, LCL1ML1=15,LCL1SL1=58,ROG=287.04/9.8
     &, NWORD=(LCL1ML+1)*60+2*LCL1SL+NSOIL*LCL1SOIL)                               
c20080701     &, NWORD=964)               !14*60+2*58+4*2=964, Eta/RSM/NCAR_WRF use this size

      PARAMETER (SNOCON=1.4594E5,RAINCON=1.1787E4)
      
C
      LOGICAL LVLWSE,SEQFLG(8),NEED
      CHARACTER*16 SEQNM1(8), SBSET
      CHARACTER*80 CLIST1(8),FMTO,ASSIGN
      CHARACTER*8 CISTAT
      DIMENSION FPACK(NWORD)

      DOUBLE PRECISION, allocatable:: PRODAT(:,:,:)
      DOUBLE PRECISION, allocatable:: LPRO(:,:),LLMH(:)

C      REAL(8) LPRO(NSTAT,NPNT), LLMH(NSTAT)
C
C     THE PURPOSE OF LPRO IS TO HOLD THE VALUES OF RISTAT UNTIL
C     THEY ARE COPIED TO FRODAT. THIS ADDITION WILL ALLOW PRODAT
C     TO BE A REAL(4) ARRAY ( AND SAVE A CONSIDERABLE AMOUNT OF 
C     MEMORY) . PRODAT CAN BE FURTHER REDUCED WITH SOME MORE EFFORT.
C                    JIM TUCCILLO
C
      REAL(8) RISTAT 
cwas  REAL(8) PRODAT(NSTAT,NPNT,NWORD),RISTAT 
      REAL(8) FRODAT(NWORD),WORKK(NWORD)


      real, allocatable, dimension(:):: P,T,U,V,Q,PFL,TFL,QFL,
     &                         CWTR,IMXR,PINT,ZINT

C      DIMENSION P(LM),T(LM),U(LM),V(LM),Q(LM),PINT(LM+1),ZINT(LM+1)
C      DIMENSION PFL(LM),TFL(LM),QFL(LM)
C      REAL CWTR(LM),IMXR(LM)
      INTEGER IDATE(3),NP1(8),NLVL(2)
      INTEGER NSTAT_TRUE,NALG
      EQUIVALENCE (CISTAT,RISTAT)
C--------------------------------------------------------------------     
C
C     SET OUTPUT UNITS FOR CLASS 1 PROFILE FILE.
C      LCLAS1 - OUTPUT UNIT FOR CLASS 1 BINARY FILE
C      LTBCL1 - INPUT UNIT FOR CLASS 1 BUFR TABLE FILE
C      LUNCL1 - OUTPUT UNIT FOR CLASS 1 BUFR FILE
C
C--------------------------------------------------------------------     
                            I N T E G E R
     & LCLAS1,LTBCL1,LUNCL1,STDOUT
C--------------------------------------------------------------------     
                            L O G I C A L
     & MONOL,BRKOUT
C--------------------------------------------------------------------     
       NAMELIST /OPTION/ MONOL,BRKOUT
                            D A T A
     & LCLAS1 / 76 /
     &,LTBCL1 / 32 /
     &,LUNCL1 / 78 /
     &,STDOUT / 6 /
     &,SEQNM1 /'HEADR','PROFILE','SURF','FLUX',
     &         'HYDR','D10M','SLYR','XTRA'/
     &,SEQFLG /.FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,
     &   .TRUE.,.FALSE./
     &,LVLWSE /.TRUE./
c      CALL W3TAGB('ETA_SNDP',1999,0267,0084,'NP22')


      PRODAT=0.
      FMTO='("ln -s ${DIRD}",I5.5,".",I4.4,3I2.2,'//
     &     '"  fort.",I2.2)'
C
C   SET MODEL TOP PRESSURE
C
       NFCST=LENGTH/INCR+1

       PTOP=2.0*100.0
        write(0,*) 'PTOP is: ', PTOP

       LM1=LM
       write(0,*) 'LM, LM1: ', LM,LM1

C allocate the needed arrays

        write(0,*) 'allocating arrays using LM: ', LM
        allocate(P(LM),T(LM),U(LM),V(LM),Q(LM),PFL(LM),TFL(LM))
        allocate(PINT(LM+1))
        allocate(ZINT(LM+1))
        allocate(QFL(LM),CWTR(LM),IMXR(LM))


        allocate(PRODAT(NWORD,NSTAT,NPNT))
        allocate(LPRO(NSTAT,NPNT),LLMH(NSTAT))

C
C   READ IN SWITCHES TO CONTROL WHETHER TO DO...
C     MONOL=.TRUE.   DO MONOLITHIC FILE
C     BRKOUT=.TRUE.  DO BREAKOUT FILES
C
       MONOL=.TRUE.
       BRKOUT=.FALSE.
       READ(11,OPTION,END=12322) !read from file eta_sndp.parm.mono
12322  CONTINUE
C
       IFCSTL=-99
       JHR=0
C----------------------------------------------------------------------
C---READ STATION DATA--------------------------------------------------
C----------------------------------------------------------------------
       write(0,*) 'LM1 is : ', LM1
      LUNIT=66   !read from binary data file 
cBZHOU      LRECPR=4*(8+9+LCL1ML1*LM1+LCL1SL1)        !for RSM/ETA
      LRECPR=4*(8+9+LCL1ML*LM1+LCL1SL)                !for NCAR-WRF

        write(0,*) 'LRECPR is: ', LRECPR

      OPEN(UNIT=LUNIT,ACCESS='DIRECT',RECL=LRECPR,IOSTAT=IER)
      NREC=0
 33   CONTINUE
      NREC=NREC+1      

!        write(0,*) 'will read NREC: ', NREC

!      DO 4000 JHR = 1, NFCST
!
!        write(*,*) ' ===== FCST ==== ', JHR
!
!       DO 3000 NST = 1, NSTAT              !Binbin: Station loop, one station on loop
!                                           !Note: the order of fcst time and stations are reversed in WRF and in ETA/RSM
!                                           !In ETA/RSM, reading recode is first one station store all fcst hours, so first read
!                                           !    station, then read fcst time
!                                           !In WRF, all stations saved in one fcst hour (i.e. one files), so first read
!                                           !    fcst hour, then read station
!                                           !So use DO 3000/4000 loop here
!      NREC=(JHR-1)*NSTAT + NST

!      READ(LUNIT,REC=NREC,IOSTAT=IRR,END=999) IHRST,IDATE,IFCST,  
      READ(LUNIT,REC=NREC,IOSTAT=IRR,ERR=999) IHRST,IDATE,IFCST,  
     &   ISTAT,CISTAT,(FPACK(N),N=1,9),(FPACK(N),N=10,FPACK(7))

        write(0,*) 'IFCST, ISTAT,CISTAT : ', IFCST, ISTAT,CISTAT

      IF(IRR.NE.0) THEN
       WRITE(*,*) NREC, '  read error, IRR=',IRR
      END IF

c      if(ISTAT.eq.725064.and.IFCST.eq.14400) then
c        write(*,*)'Just read:',NREC,IHRST,IDATE,IFCST,ISTAT,CISTAT
c        do i=414,427
c         write(*,*)i,FPACK(i)
c        end do
c      end if


        IF (IFCST.GT.IFCSTL) THEN
         INUMS=1
         JHR=JHR+1
         IFCSTL=IFCST
        ELSE
         INUMS=INUMS+1
        ENDIF

!        INUMS=NST

        IYR=IDATE(3)
        IMON=IDATE(1)
        IDAY=IDATE(2)
        RLAT=FPACK(1)
        RLON=FPACK(2)
	write(0,*) 'rlat, rlon: ', rlat, rlon
        ELEV=FPACK(3)
        LLMH(INUMS)=NINT(FPACK(4))
        LMH=NINT(FPACK(4))


        DO 25 L=1,LMH
C   REVERSE ORDER SO THAT P(1) IS THE TOP AND P(LMH) IS THE BOTTOM
         LV=LMH-L+1
         P(LV)=FPACK(L+9)
         T(LV)=FPACK(L+9+LMH)
         U(LV)=FPACK(L+9+LMH*2)
         V(LV)=FPACK(L+9+LMH*3)
         Q(LV)=FPACK(L+9+LMH*4)

!         write(0,*) 'LV, P,T,Q: ',LV, P(LV),T(LV),Q(LV)
C   IF THE CLOUD WATER IS NEGATIVE, IT IS ICE
         CWTR(LV)=FPACK(L+9+LMH*6)
         IF (CWTR(LV).LT.0.) THEN
           IMXR(LV)= -1. * CWTR(LV)
           CWTR(LV)= 0. 
           FPACK(L+9+LMH*6) = 0.
         ELSE
           IMXR(LV) = 0.
         ENDIF
 25     CONTINUE

C  USE SEA MASK TO SET SOIL/SFC VARIABLES TO MISSING VALUES
C    (IF SEA)
C
cbinbin  SM  =FPACK(13*LMH+54)              !this is for Eta
         SM  =FPACK(13*LMH+36+9+2*NSOIL+1)  !this is for both Eta(SM=54), RSM(SM=50), RWF-NCAR(54)
                                            !ie, #45, #41, #45 of sfc variables 
cZhou     IF (SM.GT.0.5) THEN                !Note all var following soil temperature needs modified
C   SMSTAV
CBZhou, WRF-NCAR has such value          FPACK(13*LMH+15)=SMISS            !soil moisture 
C
C   SUBSHX
C          FPACK(13*LMH+21)=SMISS            !sub-sfc heat flux
C   SNOPCX
C          FPACK(13*LMH+22)=SMISS            !flux  of snow phase change
C   ACSNOW, SMSTOT, SNO, ACSNOM, SSROFF, BGROFF, SOILTB
CBZhou, WRF-NCAR has ACSNOW, ACSNOM, SSROFF, BGROFF
C          DO LKJ=20,26
C            FPACK(13*LMH+LKJ+9)=SMISS
C          ENDDO
C   SFCEXC, VEGFRC, CMC, SMC(1:4), STC(1:4), (if RSM, SMC(1:2),STC(1:2) )
c Binbin modify following code:
c          DO LKJ=34,44
c            FPACK(13*LMH+LKJ+9)=SMISS
c          ENDDO
c as:
cZhou         DO LKJ=34,36
cZhou            FPACK(13*LMH+LKJ+9)=SMISS
cZhou         ENDDO
cZhou
cZhou         DO LKJ=37,37+2*NSOIL-1
cZhou            FPACK(13*LMH+LKJ+9)=SMISS
cZhou         ENDDO

ccZhou       ENDIF
C

C  GET PPT FOR CALWXT
C
         PPT  =FPACK(13*LMH+16)    !total rain
C     COMPUTE PINT,ZINT
C
        
C       Flip P, T, and Q
        DO L=1,LMH
        TFL(L)=T(L)
        QFL(L)=Q(L)
        PFL(L)=P(L)
        ENDDO

        PINT(1)=PTOP
        DO L=1,LMH-1
          DP1=PFL(L+1)-PFL(L)
          PINT(L+1)=PFL(L)+0.5*DP1
        ENDDO
        PINT(LMH+1)=PFL(LMH)+0.5*DP1
        ZINT(LMH+1)=FPACK(3)
        DO L=LMH,1,-1
         TV2=TFL(L)*(1.0+0.608*QFL(L))
         ZZ=ROG*TV2*ALOG(PINT(L+1)/PINT(L))
         ZINT(L)=ZINT(L+1)+ZZ
        ENDDO
C
C     CALL PRECIP TYPE SUBROUTINE.
C

      CALL CALWXT(TFL,QFL,PFL,PINT,LMH,LM,PPT,IWX)

      CALL CALWXT_RAMER(TFL,QFL,PFL,PINT,LMH,LM,PPT,IWX2)
      CALL CALWXT_BOURG(TFL,QFL,PINT,LMH,LM,PPT,ZINT,IWX3)
      CALL CALWXT_REVISED(TFL,QFL,PFL,PINT,LMH,LM,PPT,IWX4)

! not possible for ARWs microphysics
C!      CALL CALWXT_EXPLICIT(LMH,TSKIN,PPT,SR,RIME,IWX5)
      IWX5=0

      CALL CALWXT_DOMINANT(PPT,IWX1,IWX2,IWX3,IWX4,IWX5,
     *                     CSNO,CICE,CFZR,CRAI)

C
C
C   COMPUTE HELICITY AND STORM MOTION
C
      CALL CALHEL(U,V,P,ZINT,PINT,LMH,LM,HELI,UST,VST)
C
C   COMPUTE VISIBILITY
C   FIRST, EXTRACT THE SNOW RATIO AND SEA LEVEL PRESSURE
      SR=FPACK(9+13*LMH+49) ! set to 0 in EM 

cBZhou      SLP=FPACK(13*LMH+10)  !WRF-NCAR has no sea level pressure output
       !use sfc pressure and elevation FPACK(3)to derive sea level pressure aproximately:
       !at present, just set it = 101300.0 Pa      
C      SLP=101300.0
       SLP=FPACK(13*LMH+10)

       if (PPT.LT.0.)then
         PPT=0.0
       endif

       SNORATE=(SR/100.)*PPT/3600.
       RAINRATE=(1-(SR/100.))*PPT/3600.
       TERM1=(T(LMH)/SLP)**0.4167
       TERM2=(T(LMH)/(P(LMH)))**0.5833
       TERM3=RAINRATE**0.8333
       QRAIN=RAINCON*TERM1*TERM2*TERM3
       TERM4=(T(LMH)/SLP)**0.47
       TERM5=(T(LMH)/(P(LMH)))**0.53
       TERM6=SNORATE**0.94
       QSNO=SNOCON*TERM4*TERM5*TERM6
       TT=T(LMH)
       QV=Q(LMH)
       QCD=CWTR(LMH)
       QICE=IMXR(LMH)
       PPP=P(LMH)
 
       CALL CALVIS(QV,QCD,QRAIN,QICE,QSNO,TT,PPP,HOVI)  

C   COMPUTE CLOUD BASE PRESSURE
C   FIRST, EXTRACT THE CONVECTIVE CLOUD BASE
       HBOT=FPACK(13*LMH+9+50)
       CLIMIT =1.0E-06
       NEED = .TRUE.
       CDBP = SMISS
       CBOT = 5000                       
 
       DO L=LMH,1,-1
C GSM
C START AT THE FIRST LAYER ABOVE GROUND, AND FIND THE
C   FIRST LAYER WITH A VALUE OF CLOUD WATER GREATER THAN
C   THE SIGNIFICANT LIMIT (VALUE DESIGNATED BY Q. ZHAO).
C   THIS LAYER WILL BE THE CLOUD BOTTOM UNLESS THE BOTTOM
C   OF THE CONVECTIVE CLOUD (HBOT) IS FOUND BELOW IN WHICH
C   CASE HBOT BECOMES THE CLOUD BASE LAYER.

        
        IF ((CWTR(L)+IMXR(L)).GT.CLIMIT.AND.NEED) THEN
            CBOT=L
            IF (HBOT.GT.CBOT) THEN
              CBOT = HBOT
            ENDIF
            NEED=.FALSE.
          ENDIF
       ENDDO


       IF (CBOT.GT.LMH) THEN  !cloud base
          CDBP=SMISS
       ELSE
          CDBP=P(INT(CBOT))
       ENDIF

         
C
C
C   SET ACC/AVERAGED VARIABLES TO MISSING IF IFCST=0
C
      IF (IFCST.EQ.0) THEN
          DO L=1,LMH
           FPACK(L+9+LMH*7)=SMISS
           FPACK(L+9+LMH*8)=SMISS
          ENDDO
          DO JK=16,29
           FPACK(13*LMH+JK)=SMISS
          ENDDO
          DO JK=32,34
           FPACK(13*LMH+JK)=SMISS
          ENDDO
      ENDIF
C

C      ADD 9 SINGLE LEVEL VARIABLES TO THE OUTPUT
C      TACK THEM ON TO THE END;  WE DON'T NEED CONVECTIVE
C      CLOUD BASE, THOUGH, SO WRITE OVER THAT RECORD.
C      WRITE OVER THE DUMMY RECORDS AS WELL, SO GO BACK
C      3 PLACES (NLEN-2) FOR STARTING POINT - GSM
          NLEN = FPACK(7)         
          FPACK(NLEN+5) = CDBP     
          FPACK(NLEN-2) =   CSNO
          FPACK(NLEN-1) = CICE
          FPACK(NLEN) = CFZR
          FPACK(NLEN+1) = CRAI
          FPACK(NLEN+2) = UST
          FPACK(NLEN+3) = VST
          FPACK(NLEN+4) = HELI
          FPACK(NLEN+6) = HOVI
          FPACK(5) = FPACK(5) + 1  !add ice mixing ratio space 
          FPACK(6) = FPACK(6) + 6  !9 new variables but write over 3
          FPACK(7) = 9 + FPACK(5)*FPACK(4) + FPACK(6)


c      if(ISTAT.eq.725064.and.IFCST.eq.14400) then
c        write(*,*)'after process:',NREC,IHRST,IDATE,IFCST,ISTAT,CISTAT
c        do i=414,427
c         write(*,*)i,FPACK(i)
c        end do
c      end if


C
C           PLACE DATA INTO PRODAT IN PROPER LOCATIONS

          PRODAT (1,INUMS,JHR) = FLOAT(IFCST)
          PRODAT (2,INUMS,JHR) = FLOAT(ISTAT)

C         RISTAT is a REAL(8) variable by virtue of the fact that it 
C         is equivalenced to CISTAT. Everything else stored in PRODAT
C         is REAL(4). We have made PRODAT REAL(4) but need a REAL(8)
C         array for storing RISTAT - that is what LPRO is. Farther
C         down in the code, we will pull values out of LPRO and store
C         in FRODAT ( a REAL(8) array ).
cwas      PRODAT (3,INUMS,JHR) = RISTAT 


          LPRO   (  INUMS,JHR) = RISTAT
          PRODAT (4,INUMS,JHR) = FPACK (1)
          PRODAT (5,INUMS,JHR) = FPACK (2)
          PRODAT (6,INUMS,JHR) = FPACK (3)
          PRODAT (7,INUMS,JHR) = 1 

          DO IJ = 10, 13*LMH+9                              !13 profile variables
            PRODAT (IJ-2,INUMS,JHR) = FPACK (IJ)            !PRODAT(8,INUMS,JHR)=FPACK(10) 
          ENDDO                                             !PRODAT(9,INUMS,JHR)=FPACK(11)
C    TACK ON THE ICE WATER TO THE PROFILE SECTION           !..........
C    IT IS CURRENTLY WRITTEN IN REVERSE ORDER.
          DO L=1,LMH
            LV=LMH-L+1
            PRODAT(L+7+LMH*13,INUMS,JHR)=IMXR(LV)           !14th profile variable, ie. ice mixing ratio
          ENDDO
          DO IJ = 13*LMH+10,NLEN+6                          !all other surface variable
            PRODAT (IJ+LMH-2,INUMS,JHR) = FPACK (IJ)
          ENDDO

c      if(ISTAT.eq.725064.and.IFCST.eq.14400) then
c        write(*,*)'In prodat:',INUMS,JHR
c        do i=440,453
c         write(*,*)i,PRODAT(i,INUMS,JHR)
c        end do
c      end if

        GOTO 33

! 3000    CONTINUE
!        write(*,*) 'READ FCST ', JHR, 'done !'
! 4000    CONTINUE

 999    CONTINUE
 
        write(*,*) 'Write all of data into one big Bufr file ...'


        write(0,*) 'NREC-1: ', NREC-1
        write(0,*) 'estimated true NSTA: ', (NREC-1)/NFCST
        NSTAT_TRUE=(NREC-1)/NFCST
        write(0,*) 'NWORD is: ', NWORD
C
C  WRITE OUT INDIVIDUAL FILES FOR EACH STATION
C
        IF (BRKOUT) THEN
        DO I=1,NSTAT
         NLVL(1)=LLMH(I)
         NLVL(2)=NSOIL
C
         DO J=1,NFCST
          DO IJ = 1, NWORD
            FRODAT(IJ) = PRODAT (IJ,I,J)
          ENDDO
          FRODAT(3) = LPRO(I,J)
          ISTAT=NINT(FRODAT(2))

C
          IF (J.EQ.1) THEN
C
C     INITIALIZE BUFR LISTS SO BFRHDR WILL BE CALLED THE FIRST
C     TIME THROUGH.
C
            WRITE(ASSIGN,FMTO) ISTAT,IYR,IMON,IDAY,IHRST,LCLAS1
            CALL SYSTEM(ASSIGN)
            CLIST1(1)=' '
          ENDIF
C
C           CALL BUFR-IZING ROUTINE
C

           NSEQ = 8
           SBSET = 'ETACLS1'


          CALL BFRIZE(LTBCL1,LCLAS1,SBSET,IYR,IMON,IDAY,IHRST
     1,               SEQNM1,SEQFLG,NSEQ,LVLWSE,FRODAT,NLVL,CLIST1,NP1
     2,               WORKK,IER)
          IF(IER.NE.0)WRITE(6,1080)ISTAT,IER,FRODAT(1)
 1080   FORMAT(' SOME SORT OF ERROR ',2I8,F9.1)
C
C
         ENDDO
C
C   FINISHED, CLOSE UP BUFR FILES
C
        NSEQ = 8
        CALL BFRIZE(0,LCLAS1,SBSET,IYR,IMON,IDAY,IHRST
     1,             SEQNM1,SEQFLG,NSEQ,LVLWSE,FRODAT,NLVL,CLIST1,NP1
     2,             WORKK,IER)
         ENDDO
         ENDIF
         IF (MONOL) THEN
C
C  WRITE OUT ONE FILE FOR ALL STATIONS
C
C     INITIALIZE BUFR LISTS SO BFRHDR WILL BE CALLED THE FIRST
C     TIME THROUGH.
C
            CLIST1(1)=' '
        DO I=1,NSTAT_TRUE
         NLVL(1)=LLMH(I)
         NLVL(2)=NSOIL
C
         DO J=1,NFCST


          DO IJ = 1, NWORD
            FRODAT(IJ) = PRODAT (IJ,I,J)
          ENDDO

          FRODAT(3) = LPRO(I,J)
          ISTAT=NINT(FRODAT(2))

C
C           CALL BUFR-IZING ROUTINE
C
           NSEQ = 8
           SBSET = 'ETACLS1'
        
c        IF((I.EQ.1 .AND. J.EQ.1).OR.(I.EQ.1300 .AND. J.EQ.64)) THEN
c          write(*,*) LTBCL1,LUNCL1,SBSET,IYR,IMON,IDAY,IHRST,
c     1             SEQNM1,SEQFLG,NSEQ,LVLWSE,NLVL,NP1
c         DO IJ=1,NWORD
c          write(*,*) IJ, FRODAT(IJ)
c         END DO
c        END IF

c Binbin check
c        write(*,'(I4,3I3,2F10.0)') IYR,IMON,IDAY,IHRST,
c     &                            FRODAT(1),FRODAT(2)


          CALL BFRIZE(LTBCL1,LUNCL1,SBSET,IYR,IMON,IDAY,IHRST
     1,               SEQNM1,SEQFLG,NSEQ,LVLWSE,FRODAT,NLVL,CLIST1,NP1
     2,               WORKK,IER)
          IF(IER.NE.0)WRITE(6,1080)ISTAT,IER,FRODAT(1)

C
         ENDDO
         ENDDO
C
C   FINISHED, CLOSE UP BUFR FILES
C
        NSEQ = 8
        CALL BFRIZE(0,LUNCL1,SBSET,IYR,IMON,IDAY,IHRST
     1,             SEQNM1,SEQFLG,NSEQ,LVLWSE,FRODAT,NLVL,CLIST1,NP1
     2,             WORKK,IER)
         ENDIF
        deallocate(P,T,U,V,Q,PFL,TFL)
        deallocate(PINT)
        deallocate(ZINT)
        deallocate(QFL,CWTR,IMXR)
        deallocate(PRODAT,LPRO,LLMH)
C
        WRITE(STDOUT,*) ' END OF SOUNDING POST '
        END SUBROUTINE SNDPST_FV3S
