      SUBROUTINE MAXMIN(KFILDO,HRLY,RMXMN_6HR,ITIMEZ,CCALL,NAME,RESFLD,
     *                  ND1,IFOUR,MAXHRL,MXMN,IPRT,NSTA)
C
C        SEPTEMBER 1998   WEISS  TDL MOS-2000
C                                BASED ON A MAX/MIN CODE WRITTEN BY
C                                JACK SETTELMAIER, APRIL 1997
C        DECEMBER  2000   WEISS  MODIFIED TO ALLOW FOR THE
C                                CALCULATION OF MAX/MIN TEMPERATURES
C                                FOR PART-TIME STATIONS.
C        OCTOBER   2001   WEISS  IF STATEMENTS CONTAINING REAL VALUES
C                                CHANGED TO INTEGERS.
C        MAY       2003   WEISS  ARGUMENT LIST CHANGES; ELIMINATION 
C                                OF THE LOCAL ARRAY HFLD, CCALL1 NOW
C                                REFERENCED AS CCALL
C        MAY       2003   GLAHN  MINOR REARRANGEMENT OF COMMENTS AND
C                                TYPE STATEMENTS; WHITE SPACE 
C        SEPT      2009   VEENHU ADDED THE VARIABLE NAME TO INPUT
C                                LIST. ADDED CODE TO CHECK THE STATION
C                                STATE AND ADDED CALL TO MAXMIN_AK
C                                TO PROCESS ALASKAN STAIONS.
C
C        PURPOSE
C            THIS SUBROUTINE COMPUTES DAYTIME/NIGHTTIME MAX/MIN
C            TEMPERATURES FOR STATION (VECTOR) DATA FOR DATES 
C            STARTING ON DECEMBER 1, 1996 OR LATER. INPUTS INCLUDE
C            THE 6-HRLY MINIMUM OR MAXIMUM TEMPERATURES AND
C            HOURLY OBSERVATIONS.
C
C        DATA SET USE
C            KFILDO - DEFAULT UNIT NUMBER FOR OUTPUT (PRINT) FILE
C                     (OUTPUT).
C
C        VARIABLES
C              KFILDO = DEFAULT UNIT NUMBER FOR OUTPUT (PRINT) FILE
C                       (INPUT).
C           HRLY(N,J) = THE HOURLY TEMPERATURE OBS FOR EACH STATION
C                       WHERE N=1,ND1 AND J=1,MAXHRL (INPUT).
C      RMXMN_6HR(N,J) = ARRAY OF THE 6-HRLY MINS OR MAXS, BEFORE THEY
C                       ARE PUT IN THE ARRAYS MINS AND MAXS, WHERE
C                       N=1,ND1 AND J=1,4 (INPUT).
C           ITIMEZ(N) = ARRAY STATION'S TIME ZONE, WHERE N=1,ND1,
C                       (INPUT).
C           CCALL(N)  = AN ARRAY OF NSTA NUMBER OF STATION CALL
C                       LETTERS USED FOR DIAGNOSTICS, WHERE N=1,ND1
C                       (INPUT).
C           RESFLD(N) = THE RESULTANT FIELD (MAX OR MIN) FOR EACH
C                       STATION, WHERE N=1,ND1 (OUTPUT).
C                 ND1 = MAXIMUM NUMBER OF STATIONS (INPUT).
C               IFOUR = VALUE OF 4, MAXIMUM NUMBER OF 6-HRLY MINS OR 
C                       MAXS IN A GIVEN DAY (INPUT).
C              MAXHRL = MAXIMUM NUMBER OF HOURS (=25) USED TO STORE
C                       THE HOURLY TEMPERATURES FROM ACROSS ALL
C                       POTENTIAL TIME ZONES (INPUT).
C                MXMN = FLAG TO INDICATE DAYTIME MAX (=2) OR
C                       NIGHTTIME MIN (=1) PROCESSING (INPUT).
C             NAME(K) = NAMES OF STATIONS (K=1,NSTA).  USED FOR PRINTOUT
C                       ONLY.  (CHARACTER*20)  (INPUT)
C                IPRT = PRINT FLAG USED TO (IPRT > 0), PRINT OUT
C                       INCREASING AMOUNTS OF DIAGNOSTICS
C                       1=SOME(RECOMMENDED), 2=ALL DIAGNOSTICS (INPUT).
C                NSTA = NUMBER OF STATIONS THAT THE MAX OR MIN WILL BE 
C                       CALCULATED FOR (INPUT).
C
C        OTHER VARIABLES
C
C                 
C                CMAX = LOGICAL PARAMETER FOR DAYTIME MAXIMUM
C                       PROCESSING (INTERNAL).
C                CMIN = LOGICAL PARAMETER FOR NIGHTTIME MINIMUM
C                       PROCESSING (INTERNAL).
C             HOLD(K) = ARRAY USED TO HOLD THE MIN/MAX VALUES (MIN1/MAX1
C                       VS SEGMENT B) AND (MIN3/MAX3 VS SEGMENT C). 
C                       MIN/MAX PROCESSING USES THESE VALUES ONLY IF 
C                       SIMPLE CHECKING FAILS, WHERE K=1,4 (INTERNAL).
C             IMSG(L) = ARRAY TO HOLD THE NUMBER OF MISSINGS IN A RANGE
C                       OF DATA; USUALLY HOURLY, FOR A GIVEN STATION,
C                       WHERE L=1,3 (OUTPUT). 
C              INUMZN = NUMBER OF HOURS (=14) USED TO SELECT HOURLY
C                       OBS REPRESENTING THE FIRST AND THIRD 6-HOUR
C                       TIME PERIODS (INTERNAL).
C           ITIME1(K) = THE STARTING POINT OF A STATION'S 2 HOUR WINDOW
C                       CHECK (6AM, 6PM: MAX, AND 6PM,7AM: MIN), A - B
C                       AND C - D BOUNDARIES, WHERE K=1,4 (INTERNAL).
C           ITIME2(K) = THE ENDING POINT OF A STATION'S 2 HOUR WINDOW
C                       CHECK (8AM, 8PM: MAX, AND 8PM,9AM: MIN), A - B
C                       AND C - D BOUNDARIES, WHERE K=1,4 (INTERNAL).
C           ITIMEZ_AK = USE TO INPUT TIME ZONE INFORMATION TO MAXMIN_AK
C                       (INTERNAL).
C              IWK(M) = THE ARRAY ELEMENT OF RMXMN_6HR TO CHOOSE IN
C                       FILLING ARRAYS MINS AND MAXS DEPENDING ON 
C                       TIME ZONE, WHERE M=1,3 (INTERNAL).
C                       MAX ATLANTIC: 12Z,18Z,00Z  OTHERS: 18Z,00Z,06Z
C                       MIN ATLANTIC: 00Z,06Z,12Z  OTHERS: 06Z,12Z,18Z
C                  IZ = DETERMINES STARTING POINT IN HRLY ARRAY FOR 
C                       RECONSTRUCTING 6 HOUR MAX/MIN VALUES BASED ON 
C                       HOURLY OBS, DEPENDING ON TIME ZONE (INTERNAL).
C                JEND = THE ENDING POINT OF THE RANGE OF HOURLY OBS
C                       MAKING UP THE TIME PERIOD OF A 6-HR MIN/MAX
C                       ESTIMATE (INTERNAL).
C              JSTART = THE STARTING POINT OF THE RANGE OF HOURLY OBS
C                       MAKING UP THE TIME PERIOD OF A 6-HR MIN/MAX
C                       ESTIMATE (INTERNAL).
C                KEND = ENDING ARRAY ELEMENT ASSOCIATED WITH KSTART 
C                       (SEE KSTART) (INTERNAL). 
C              KSTART = STARTING ARRAY ELEMENT TO START CHECKING RMXABCD
C                       OR RMNABCD FOR MISSING VALUES WHEN NECESSARY.
C                       A STATION'S PART-TIME STATUS WILL DETERMINE THE
C                       VALUE (INTERNAL).
C           MAXHRL_AK = MAXIMUM NUMBER OF HOURS (=19) USED TO STORE
C                       THE HOURLY TEMPERATURES FOR PROCESSING OF 
C                       ALASKA MAX/MIN TEMPERATURE (INTERNAL).
C            NUM_MISS = FOR PART-TIME MAX AND MIN ESTIMATES, THE NUMBER
C                       OF MISSING HOURLY OBS FOR SEGMENTS B AND C
C                       (MAXMIN_PART) FOR A GIVEN STATION (INTERNAL).
C              RECON1 = FOR PART-TIME STATIONS, LOGICAL PARAMETER USED 
C                       TO INDICATE EITHER MAX1/MIN1 HAS BEEN
C                       GENERATED SUCCESSFULLY (INTERNAL).
C              RECON3 = FOR PART-TIME STATIONS, LOGICAL PARAMETER USED
C                       TO INDICATE EITHER MAX3/MIN3 HAS BEEN
C                       GENERATED SUCCESSFULLY (INTERNAL).
C               RFILL = LOGICAL PARAMETER INDICATING THAT SUBROUTINE
C                       RMXMN_FILL HAS BEEN CALLED AT LEAST ONCE
C                       (ALL THAT IS REQUIRED) (INPUT/INTERNAL).
C          RMAXS(N,L) = WORK ARRAY USED TO STORE THE THREE 6-HRLY MAXS 
C                       USED FOR DAYTIME MAX ESTIMATION, WHERE 
C                       N=1,ND1 AND L=1,3 (INTERNAL).
C          RMINS(N,L) = WORK ARRAY USED TO STORE THE THREE 6-HRLY MINS 
C                       USED FOR NIGHTTIME MIN ESTIMATION, WHERE 
C                       N=1,ND1 AND L=1,3 (INTERNAL).
C          RMNABCD(J) = HOLDS THE SEGMENTED MIN TEMPERATURES OF THE 
C                       FIRST AND THIRD 6 HOUR TEMPERATURES A,B AND 
C                       C,D RESPECTIVELY. B AND C ARE THE STARTING 
C                       ENDING SEGMENTS OF THE NIGHTTIME MIN WINDOW 
C                       FOR A GIVEN STATION, WHERE J=1,4 (INTERNAL). 
C          RMNHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST 
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (INTERNAL).  
C          RMXABCD(J) = HOLDS THE SEGMENTED MAX TEMPERATURES OF THE 
C                       FIRST AND THIRD 6 HOUR TEMPERATURES A,B AND 
C                       C,D RESPECTIVELY. B AND C ARE THE STARTING
C                       ENDING SEGMENTS OF THE DAYTIME MAX WINDOW 
C                       FOR A GIVEN STATION, WHERE J=1,4 (INTERNAL).
C          RMXHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MAX TEMPERATURE SEGMENTS A,B,C AND D, 
C                       WHERE K=1,14 (INTERNAL).
C             RMXMN13 = THE MAX1/MIN1 AND/OR MAX3/MIN3 VALUE RETURNED
C                       SUBROUTINE MAXMIN_PART, FOR PART-TIME STATION
C                       MAX/MIN ESTIMATES (INTERNAL).
C     RMXMN_6HR_AK(J) = ARRAY OF THE 6-HRLY MINS OR MAXS PASSED TO
C                       MAXMIN_AK WHERE J=1,4 (INTERNAL).
C           RESFLD_AK = THE RESULTING MAX/MIN FROM THE MAXMIN_AK
C                       SUBROUTINE.
C               STATE = CHARACTER LEN=2 WHICH CONTAINS THE STATE OF 
C                       STATION BEING PROCESSED.
C               ZMISS = MISSING VALUE OF 9999. (INTERNAL). 
C
C        INTERNAL SUBROUTINES
C             CKFRMSG = SUBROUTINE TO CHECK HRLY ARRAY FOR NUMBER OF 
C                       MISSINGS OVER 6 HOUR PERIODS OR SPECIFIC
C                       HOURLY OBSERVATIONS. 
C
C              CKABCD = SUBROUTINE TO CHECK THE NUMBER OF MISSING HOURLY 
C                       OBS FOR SEGMENTS A,B,C AND D.
C
C          RMXMN_FILL = SUBROUTINE TO FILL THE ARRAYS RMXHRLY AND 
C                       RMNHRLY WITH HOURLY TEMPERATURE OBS WHICH
C                       ARE THEN USED LATER TO GENERATE 
C                       SEGMENTS A,B,C AND D.
C
C           MAXMIN_AK = SUBROUTINE TO CALCULATE MAX/MIN FOR ALASKA
C                       STATIONS. ALASKA STATION REQUIRED DIFFERENT
C                       PROCESSING DUE TO THE MAX/MIN TIME WINDOW
C                       CHANGE IMPLIMENTED BY THE ALASKA REGION.
C
C         MAXMIN_PART = SUBROUTINE TO CALCULATE (RECONSTRUCT) THE 
C                       MAX1/MIN1 AND/OR MAX3/MIN3 FOR PART-TIME 
C                       STATIONS.
C
      IMPLICIT NONE
C  
      CHARACTER*8 CCALL(ND1),CCALL_AK
      CHARACTER*20 NAME(ND1),TEMP_NAME
      CHARACTER*2 STATE
C
      LOGICAL CMAX,CMIN,RECON1,RECON3,RFILL
C
      INTEGER, PARAMETER :: INUMZN=14
C
      INTEGER ITIMEZ(NSTA),ITIMEZ_AK
      INTEGER KFILDO,ND1,IFOUR,MAXHRL,MXMN,IPRT,NSTA,JSTART,JEND,
     1        IZ,I,J,NUM_MISS,KSTART,KEND,IZMISS
      INTEGER MAXHRL_AK
      INTEGER IWK(3),IMSG(3),ITIME1(4),ITIME2(4)
C
      REAL RMAXS(ND1,3),RMINS(ND1,3)
C        RMAXS( , ) AND RMINS( , ) ARE AUTOMATIC ARRAYS.
      REAL RMXHRLY(INUMZN),RMNHRLY(INUMZN),RESFLD(ND1)
      REAL RMXABCD(4),RMNABCD(4),HOLD(4),RESFLD_AK
      REAL RMXMN_6HR(ND1,IFOUR),HRLY(ND1,MAXHRL)
      REAL RMXMN13,ZMISS,HRLY_AK(19),RMXMN_6HR_AK(3)
C
      DATA IMSG/3*0/
      MAXHRL_AK=19
C
C******************************************************
C
C        STEP 1A. INITIALIZE RESULTANT DAYTIME/NIGHTTIME MAX/MIN
C
      ZMISS=9999.
      IZMISS=9999
C
      DO 4 I=1,ND1
!         PRINT*,"Loop number 4.  I,ND1,RESFLD=",I,ND1,RESFLD(I)
        RESFLD(I)=-ZMISS
  4   CONTINUE
C
      CMIN=.FALSE.
      CMAX=.FALSE.
      IF(MXMN.EQ.1)CMIN=.TRUE.
      IF(MXMN.EQ.2)CMAX=.TRUE.
      PRINT*, "ENTERING MAXMINT",NSTA
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        *** MAIN DO LOOP ***
C        *** FOR EACH STATION..............
C
      DO 300 I=1,NSTA
C         PRINT*,"LOOKING AT STATION NUMBER:",I
        IF((IPRT.GE.1).AND.(I.EQ.1)) WRITE(KFILDO,*) '*******'
        IF((IPRT.GE.1).AND. (I.EQ.1))WRITE(KFILDO,*) 'NEW DAY'
        IF((IPRT.GE.1).AND. (I.EQ.1))WRITE(KFILDO,*) '*******'
C
C
C       FIRST CHECK WHICH STATE THE STATION IS IN. IN JULY 2009
C       THE CODE WAS CHANGED TO HANDEL A NEW MAX/MIN TIME 
C       WINDOW DEFINITION USED BY THE ALASKAN REGION. SINCE THE
C       CHANGES APPLY ONLY TO ALASKAN STATIONS A SEPARATE 
C       MAXMIN_AK SUBROUTINE WAS WRITTEN.
C
C
        TEMP_NAME=NAME(I)
C        STATE=TEMP_NAME(19:20)
        STATE=TEMP_NAME(19:20)
C       CHANGE BY LEVINE AFTER LOOKING AT MDL METAR DICTIONARY
C 
        IF(STATE == '  ') THEN
C           PRINT *,'NO STATE FOR STATION: '  
           WRITE(KFILDO,*) 'NO STATE FOR STATION:',CCALL(I),NAME(I)
        ELSE
           WRITE(KFILDO,*) 'STATION,STATE=',CCALL(I),STATE
        ENDIF
C
C       IF THE STATE EQUALS AK THEN CALL MAXMIN_AK
C
        IF(STATE == 'AK') THEN
C
          HRLY_AK(1:19)=HRLY(I,7:25)
          RMXMN_6HR_AK(1:3)=RMXMN_6HR(I,2:4)
          ITIMEZ_AK=ITIMEZ(I)
          CCALL_AK=CCALL(I)
          RESFLD_AK=9999.
C
          CALL MAXMIN_AK(KFILDO,HRLY_AK,RMXMN_6HR_AK,
     1                     ITIMEZ_AK,CCALL_AK,RESFLD_AK,
     2                     MAXHRL_AK,MXMN,IPRT)
 
          RESFLD(I)=RESFLD_AK
          GOTO 300
        ENDIF
C
C        STEP 1B.INITIALIZE THE RECON1, RECON3 AND RFILL FOR 
C                EACH STATION
        RECON1=.FALSE.
        RECON3=.FALSE.
        RFILL=.TRUE.
C
C        STEP 2. MAKE SURE THE TIME ZONE IS A VALID ONE 
C                IF NOT, A SPECIAL MISSING VALUE WILL BE ASSIGNED 
C                TO THE MAX/MIN. (FORMERLY AT STATEMENT NUMBER 29)
C
        IF(ITIMEZ(I).LE.-11.OR.ITIMEZ(I).GE.-3) THEN
          WRITE(KFILDO,*) 'INVALID TIME ZONE= ',
     *    ITIMEZ(I),' FOR STATION= ',CCALL(I)
	  RESFLD(I)=ZMISS
          GO TO 300
       ELSE
          WRITE(KFILDO,*) 'STATION,TIMEZONE=',CCALL(I),ITIMEZ(I)
       ENDIF
C
C         STEP 3. PREPARE TO CHECK IF THE 6 HOUR VALUES ARE MISSING AND 
C                 REPLACE IF POSSIBLE (ONLY IF ONE OR LESS MISSING 
C                 HOURLY VALUES).  
C
C         STEP 3A. SET JSTART, JEND AND IZ TO PROCESS DIFFERENT 
C                  TIME ZONES.
C
C                  FOR MAIN US, AK, AND HI TIME ZONES
	IF(ITIMEZ(I).LE.-5) THEN
	  JSTART=7
	  JEND=25
	  IZ=0
	  IWK(1)=2
	  IWK(2)=3
	  IWK(3)=4
C                FOR PUERTO RICO'S TIME ZONE 
        ELSEIF(ITIMEZ(I).EQ.-4)THEN
	  JSTART=1
	  JEND=19
	  IZ=6
	  IWK(1)=1
	  IWK(2)=2
	  IWK(3)=3
        ENDIF 
C
C        STEP 3B. FOR NIGHTIME MIN DO THE FOLLOWING:
C
C              1. FILL UP WORK ARRAY (MINS) USING THE INPUT ARRAY 
C                 RMXMN_6HR BY COPYING THE THREE 6-HRLY MINS TO BE
C                 USED TO CALCULATE THE NIGHTTIME MIN.
        IF(CMIN) THEN
          RMINS(I,1)=RMXMN_6HR(I,IWK(1))
          RMINS(I,2)=RMXMN_6HR(I,IWK(2))
          RMINS(I,3)=RMXMN_6HR(I,IWK(3))
C
C              2. CHECK THE HOURLY OBS WITHIN THE JSTART AND JEND RANGE
C                 FOR MISSINGS BY CALLING CKFRMSG. THE HOURLY OBS MAY
C                 HAVE TO BE USED TO RECREATE A MISSING 6-HRLY MIN. 
          CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,JSTART,JEND,6,ZMISS,IMSG)
C
C              3. CHECK EACH OF THE 6-HRLY MINS FOR MISSING VALUES.
C                 IF ONE 6-HRLY MIN IS MISSING, AND THERE ARE ENOUGH
C                 HOURLY OBS TO REGENERATE THE 6-HRLY MIN, THEN DO SO.
C                 IF THERE IS NO 6-HRLY MIN AND NO WAY TO RECREATE ONE,
C                 THEN SET THE NIGHTTIME MIN TO MISSING AND GO 
C                 TO THE NEXT STATION. 
C
	  DO 210 J=1,3
C
            IF(NINT(RMINS(I,J)).EQ.IZMISS) THEN
C
              IF(IMSG(J).GE.2) THEN
C
C                  STEP 3BB. PART-TIME STATION MIN1 AND/OR MIN3
C                            RECONSTRUCTION
C
                IF((J.EQ.1).OR.(J.EQ.3)) THEN
C                   SPECIAL CALL OF RMXMN_FILL (SEE STEP 5A.) 
                  IF(RFILL) CALL RMXMN_FILL(INUMZN,I,IZ,MAXHRL,ND1,
     *                           HRLY,RMXHRLY,RMNHRLY,ZMISS,
     *                           CMAX,CMIN,RFILL)
                  NUM_MISS=0
                  RMXMN13=ZMISS
                  CALL MAXMIN_PART(INUMZN,ITIMEZ(I),J,
     *                             RMXHRLY,RMNHRLY,RMXMN13,ZMISS,
     *                             CMAX,CMIN,RECON1,RECON3,NUM_MISS,
     *                             IPRT,KFILDO)
C
                  IF(NINT(RMXMN13).EQ.IZMISS) THEN
                    RESFLD(I)=ZMISS
                    IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',
     *                CCALL(I),'TOO FEW OBS TO MAKE A 6HR', 
     *                ' MIN',J,': PART-T MISS=',NUM_MISS,
     *                ' FULL-T MISS=',IMSG(J)
                    GO TO 300
                  ELSE
                    RMINS(I,J)=RMXMN13
                  ENDIF
C
                ELSEIF(J.EQ.2) THEN 
                  RESFLD(I)=ZMISS
                  IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',
     *                          CCALL(I),'NOT ENOUGH OBS TO MAKE',
     *                          ' A 6HR MIN',J,': MISSING#=',IMSG(J)
                  GO TO 300
                ENDIF
C
              ELSE
                RMINS(I,J)=MIN(HRLY(I,J*6+1-IZ),HRLY(I,J*6+2-IZ),
     *                        HRLY(I,J*6+3-IZ),HRLY(I,J*6+4-IZ),
     *                        HRLY(I,J*6+5-IZ),HRLY(I,J*6+6-IZ),
     *                        HRLY(I,J*6+7-IZ))
                IF(IPRT.GE.1) WRITE(KFILDO,*) ' FOR STATION= ',
     *                        CCALL(I),'6HR MIN',J,' HAD TO BE REMADE'
              ENDIF
C
            ENDIF
C
  210     CONTINUE
C
          IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *    'MIN1-3=',RMINS(I,1),' ',RMINS(I,2),' ',RMINS(I,3)
        ENDIF
C
C         STEP 3C. FOR DAYTIME MAX DO THE FOLLOWING:
C
C              1. FILL UP WORK ARRAY (MAXS) USING THE INPUT ARRAY
C                 RMXMN_6HR BY COPYING THE THREE 6-HRLY MAXS TO BE
C                 USED TO CALCULATE THE DAYTIME MAX.
C
        IF(CMAX) THEN
          RMAXS(I,1)=RMXMN_6HR(I,IWK(1))
          RMAXS(I,2)=RMXMN_6HR(I,IWK(2))
          RMAXS(I,3)=RMXMN_6HR(I,IWK(3))
          IF(IPRT.EQ.2) THEN
             WRITE(KFILDO,*) 'STATION, ORIGINAL MAXES=',CCALL(I),
     *       RMAXS(I,1),RMAXS(I,2),RMAXS(I,3)
          ENDIF
C
C              2. CHECK THE HOURLY OBS WITHIN THE JSTART AND JEND RANGE
C                 FOR MISSINGS BY CALLING CKFRMSG. THE HOURLY OBS MAY
C                 HAVE TO BE USED TO RECREATE A MISSING 6-HRLY MAX.
C              ** ALSO, ZMISS IS CHANGED TO -ZMISS IN THE HRLY ARRAY.**
C
          CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,JSTART,JEND,6,ZMISS,IMSG)
C
C              3. CHECK EACH OF THE 6-HRLY MAXS FOR MISSING VALUES.
C                 IF ONE 6-HRLY MAX IS MISSING, AND THERE ARE ENOUGH
C                 HOURLY OBS TO REGENERATE THE 6-HRLY MAX, THEN DO SO.
C                 IF THERE IS NO 6-HRLY MAX AND NO WAY TO RECREATE ONE,
C                 THEN SET THE DAYTIME MAX TO MISSING AND GO
C                 TO THE NEXT STATION.
C
          DO 215 J=1,3
C
            IF(NINT(RMAXS(I,J)).EQ.IZMISS) THEN
C
              IF(IMSG(J).GE.2) THEN
C
C                  STEP 3CC. PART-TIME STATION MAX1 AND/OR MAX3
C                            RECONSTRUCTION
C
                IF((J.EQ.1).OR.(J.EQ.3)) THEN
C                      SPECIAL CALL OF RMXMN_FILL (SEE STEP 5A.) 
                  IF(RFILL) CALL RMXMN_FILL(INUMZN,I,IZ,MAXHRL,ND1,
     *                           HRLY,RMXHRLY,RMNHRLY,ZMISS,
     *                           CMAX,CMIN,RFILL)
                  NUM_MISS=0
                  RMXMN13=-ZMISS
                  CALL MAXMIN_PART(INUMZN,ITIMEZ(I),J,
     *                             RMXHRLY,RMNHRLY,RMXMN13,ZMISS,
     *                             CMAX,CMIN,RECON1,RECON3,NUM_MISS,
     *                             IPRT,KFILDO)
C
                  IF(NINT(RMXMN13).EQ.-IZMISS) THEN
                    RESFLD(I)=ZMISS
                    IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',
     *                CCALL(I),'TOO FEW OBS TO MAKE A 6HR', 
     *                ' MAX',J,': PART-T MISS=',NUM_MISS,
     *                ' FULL-T MISS=',IMSG(J)
                    GO TO 300
                  ELSE
                    RMAXS(I,J)=RMXMN13
                  ENDIF
C
                ELSEIF(J.EQ.2) THEN 
                  RESFLD(I)=ZMISS
                  IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',
     *                          CCALL(I),'NOT ENOUGH OBS TO MAKE',
     *                          ' A 6HR MAX',J,': MISSING#=',IMSG(J)
                  GO TO 300
                ENDIF
C
              ELSE
                RMAXS(I,J)=MAX(HRLY(I,J*6+1-IZ),HRLY(I,J*6+2-IZ),
     *                        HRLY(I,J*6+3-IZ),HRLY(I,J*6+4-IZ),
     *                        HRLY(I,J*6+5-IZ),HRLY(I,J*6+6-IZ),
     *                        HRLY(I,J*6+7-IZ))
                IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',
     *          CCALL(I),'6HR MAX',J,' HAD TO BE REMADE'
              ENDIF
C
            ENDIF
C
  215     CONTINUE
C
          IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *    'MAX1-3=',RMAXS(I,1),' ',RMAXS(I,2),' ',RMAXS(I,3)
        ENDIF
C
C         STEP 4. SIMPLE CHECK OF THE 6-HR VALUES.
C                 TEST IF MIDDLE 6 HOUR MAX/MIN IS (HIGHEST/LOWEST).
C
C         STEP 4A. IF SIMPLE CHECK FOR MIN IS VALID, SET 
C                  NIGHTTIME MIN AND GO ON TO NEXT STATION.
        IF(CMIN) THEN
C
          IF(NINT(RMINS(I,2)).LE.NINT(RMINS(I,1)).AND.
     *      NINT(RMINS(I,2)).LE.NINT(RMINS(I,3))) THEN
            RESFLD(I)=RMINS(I,2)
            IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *                    'MINS EASY WAY= ',RESFLD(I)
            GO TO 300
          ELSE 
C                 IF SIMPLE CHECK IS NOT VALID, SET NIGHTTIME MIN 
C                 TO A UNIQUE MISSING VALUE AND GO TO STEP 5.
            IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *                    'MINS:NEED TO CHECK INDIVIDUAL HOURS'
             RESFLD(I)=-888.
          ENDIF
C
        ENDIF
C
C        STEP 4B. IF SIMPLE CHECK FOR MAX IS VALID, SET 
C                 DAYTIME MAX AND GO ON TO NEXT STATION.
C
        IF(CMAX) THEN
C
          IF(NINT(RMAXS(I,2)).GE.NINT(RMAXS(I,1)).AND.
     *      NINT(RMAXS(I,2)).GE.NINT(RMAXS(I,3)))THEN
            RESFLD(I)=RMAXS(I,2)
            IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *                    ' MAXS EASY WAY=',RESFLD(I)
            GO TO 300
          ELSE 
C                 IF SIMPLE CHECK IS NOT VALID, SET DAYTIME MAX
C                 TO A UNIQUE MISSING VALUE AND GO TO STEP 5.
            IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *                    'MAXS:NEED TO CHECK INDIVIDUAL HOURS'
             RESFLD(I)=-777.
          ENDIF
C
        ENDIF
C
C         STEP 5. THE SIMPLE CHECK WAS INVALID, THEREFORE MUST CHECK 
C                 6 HOUR PERIODS ONE AND THREE. THE 6 HOUR PERIODS ARE
C                 SEGMENTED AS A & B FOR PERIOD ONE AND C & D FOR 
C                 PERIOD 3. B AND C CORRESPOND TO BEGINNING AND 
C                 END OF THE "WINDOW" HOURLY PORTION OF THE 
C                 DAYTIME MAX OR NIGHTTIME MIN. 
C
C         STEP 5A. IF THE SIMPLE CHECK OF THE 6-HRLY VALUES
C                  WAS NOT CONCLUSIVE, THEN RMXHRLY AND RMNHRLY
C                  ARE FILLED WITH THE HOURLY TEMPERATURE OBS USED TO
C                  GENERATE A,B,C AND D FOR EITHER MAX OR MIN.
C
        IF (RFILL) CALL RMXMN_FILL(INUMZN,I,IZ,MAXHRL,ND1, 
     *                    HRLY,RMXHRLY,RMNHRLY,ZMISS,
     *                    CMAX,CMIN,RFILL)
C
C         STEP 5B. SET UP THE THE MAX OR MIN VALUES OF SEGMENTS
C                   A, B, C, AND D FOR EACH TIME ZONE.
C                   MIN WINDOW = 7PM - 8AM 
C                   MAX WINDOW = 7AM - 7PM
C               1. ALEUTIANS/HAWAII
C
        IF(ITIMEZ(I).EQ.-10) THEN
          RMNABCD(1)=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3),
     *                   RMNHRLY(4),RMNHRLY(5))
          RMNABCD(2)=MIN(RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),RMNHRLY(11),
     *                   RMNHRLY(12),RMNHRLY(13),RMNHRLY(14))
          RMNABCD(4)=9988.
          RMXABCD(1)=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3),
     *                   RMXHRLY(4),RMXHRLY(5))
          RMXABCD(2)=MAX(RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                   RMXHRLY(11),RMXHRLY(12),RMXHRLY(13))
          RMXABCD(4)=RMXHRLY(14)
C              2. ALASKA
        ELSEIF(ITIMEZ(I).EQ.-9) THEN
          RMNABCD(1)=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3),
     *                   RMNHRLY(4))
          RMNABCD(2)=MIN(RMNHRLY(5),RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *                   RMNHRLY(11),RMNHRLY(12),RMNHRLY(13))
          RMNABCD(4)=RMNHRLY(14)
          RMXABCD(1)=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3),
     *                   RMXHRLY(4))
          RMXABCD(2)=MAX(RMXHRLY(5),RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                   RMXHRLY(11),RMXHRLY(12))
          RMXABCD(4)=MAX(RMXHRLY(13),RMXHRLY(14))
C              3.             PACIFIC
        ELSEIF(ITIMEZ(I).EQ.-8) THEN
          RMNABCD(1)=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3))
          RMNABCD(2)=MIN(RMNHRLY(4),RMNHRLY(5),RMNHRLY(6),
     *                   RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *                   RMNHRLY(11),RMNHRLY(12))
          RMNABCD(4)=MIN(RMNHRLY(13),RMNHRLY(14))
          RMXABCD(1)=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3))
          RMXABCD(2)=MAX(RMXHRLY(4),RMXHRLY(5),RMXHRLY(6),
     *                   RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                   RMXHRLY(11))
          RMXABCD(4)=MAX(RMXHRLY(12),RMXHRLY(13),RMXHRLY(14))
C              4. MOUNTAIN
        ELSEIF(ITIMEZ(I).EQ.-7) THEN
          RMNABCD(1)=MIN(RMNHRLY(1),RMNHRLY(2))
          RMNABCD(2)=MIN(RMNHRLY(3),RMNHRLY(4),RMNHRLY(5),
     *                   RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *                   RMNHRLY(11))
          RMNABCD(4)=MIN(RMNHRLY(12),RMNHRLY(13),RMNHRLY(14))
          RMXABCD(1)=MAX(RMXHRLY(1),RMXHRLY(2))
          RMXABCD(2)=MAX(RMXHRLY(3),RMXHRLY(4),RMXHRLY(5),
     *                   RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10))
          RMXABCD(4)=MAX(RMXHRLY(11),RMXHRLY(12),RMXHRLY(13),
     *                   RMXHRLY(14))
C              5. CENTRAL
        ELSEIF(ITIMEZ(I).EQ.-6) THEN
          RMNABCD(1)=RMNHRLY(1)
          RMNABCD(2)=MIN(RMNHRLY(2),RMNHRLY(3),RMNHRLY(4),RMNHRLY(5),
     *                   RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10))
          RMNABCD(4)=MIN(RMNHRLY(11),RMNHRLY(12),RMNHRLY(13),
     *                   RMNHRLY(14))
          RMXABCD(1)=RMXHRLY(1)
          RMXABCD(2)=MAX(RMXHRLY(2),RMXHRLY(3),RMXHRLY(4),RMXHRLY(5),
     *                   RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9))
          RMXABCD(4)=MAX(RMXHRLY(10),RMXHRLY(11),RMXHRLY(12),
     *                   RMXHRLY(13),RMXHRLY(14))
C              6. EASTERN
        ELSEIF(ITIMEZ(I).EQ.-5) THEN
          RMNABCD(1)=9988.
          RMNABCD(2)=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3),RMNHRLY(4),
     *                   RMNHRLY(5),RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9))
          RMNABCD(4)=MIN(RMNHRLY(10),RMNHRLY(11),RMNHRLY(12),
     *                   RMNHRLY(13),RMNHRLY(14))
          RMXABCD(1)=-9988.
          RMXABCD(2)=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3),RMXHRLY(4),
     *                   RMXHRLY(5),RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=RMXHRLY(8)
          RMXABCD(4)=MAX(RMXHRLY(9),RMXHRLY(10),RMXHRLY(11),
     *                   RMXHRLY(12),RMXHRLY(13),RMXHRLY(14))
C              7. ATLANTIC
        ELSEIF(ITIMEZ(I).EQ.-4) THEN
          RMNABCD(1)=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3),RMNHRLY(4),
     *                   RMNHRLY(5))
          RMNABCD(2)=MIN(RMNHRLY(6),RMNHRLY(7))
          RMNABCD(3)=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),RMNHRLY(11),
     *                   RMNHRLY(12),RMNHRLY(13),RMNHRLY(14))
          RMNABCD(4)=9988.
          RMXABCD(1)=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3),RMXHRLY(4),
     *                   RMXHRLY(5))
          RMXABCD(2)=MAX(RMXHRLY(6),RMXHRLY(7))
          RMXABCD(3)=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),RMXHRLY(11),
     *                   RMXHRLY(12),RMXHRLY(13))
          RMXABCD(4)=RMXHRLY(14)
        ENDIF
C
C        STEP 5C. SET UP 2 HOUR BREAKPOINTS AROUND EDGES OF SEGMENTS
C                 A,B,C AND D FOR MIN/MAX "TREND". COMPARISONS
C                 USING THESE HOURLY TEMPERATURES WORKS ESSENTIALLY 
C                 LIKE A "TIE BREAKER" IF A=B AND/OR C=D.
C                 FOR EXAMPLE, A=B AND/OR C=D.
C                 MIN WINDOW = 7PM - 8AM
C                 MAX WINDOW = 7AM - 7PM

C        FOR USE WITH AB PART OF MIN (6PM VS 8PM)
                ITIME1(1)=11+ITIMEZ(I)
                ITIME2(1)=13+ITIMEZ(I)
C
C        FOR USE WITH CD PART OF MIN (7AM VS 9AM)
                ITIME1(2)=24+ITIMEZ(I)
                ITIME2(2)=26+ITIMEZ(I)
C
C        FOR USE WITH AB PART OF MAX (6AM VS 8AM)
                ITIME1(3)=11+ITIMEZ(I)
                ITIME2(3)=13+ITIMEZ(I)
C
C        FOR USE WITH CD PART OF MAX (6PM VS 8PM)
                ITIME1(4)=23+ITIMEZ(I)
                ITIME2(4)=25+ITIMEZ(I)
C
C        STEP 5D. CALL SUBROUTINE CKABCD AND DETERMINE IF THE NUMBER 
C                 OF MISSING PER SEGMENT A,B,C, OR D IS LESS THAN 2.
C
        CALL CKABCD(RMNHRLY,RMXHRLY,RMNABCD,RMXABCD,ITIMEZ(I),
     *              ZMISS,INUMZN)
C
        IF(CMAX) THEN
C
          IF((.NOT. RECON1).AND.(.NOT. RECON3))THEN
	    KSTART=1
	    KEND=4
          ELSEIF((RECON1).AND.(.NOT. RECON3))THEN
	    KSTART=3
	    KEND=4
	    HOLD(3)=RMAXS(I,1) 
	  ELSEIF((RECON3).AND.(.NOT. RECON1))THEN
	    KSTART=1
	    KEND=2
	    HOLD(4)=RMAXS(I,3)
          ELSEIF((RECON1).AND.(RECON3))THEN
	    HOLD(3)=RMAXS(I,1)
	    HOLD(4)=RMAXS(I,3)
            KSTART=1
            KEND=0
          ENDIF
C
          IF(IPRT.GE.1) WRITE(KFILDO,*) ' KSTART=',KSTART,' KEND=', KEND
C
          DO 245 J=KSTART,KEND
C
            IF(NINT(RMXABCD(J)).EQ.-IZMISS) THEN
	      RESFLD(I)=ZMISS
              IF(IPRT.GE.1) WRITE(KFILDO,*) 'RMXABCD(',J,') IS MISSING'
            ENDIF
C
 245      CONTINUE
C
          IF(NINT(RESFLD(I)).EQ.IZMISS)GO TO 300
C
        ELSEIF(CMIN) THEN
C
	  IF((.NOT. RECON1).AND.(.NOT. RECON3))THEN
	    KSTART=1
	    KEND=4
          ELSEIF((RECON1).AND.(.NOT. RECON3))THEN
	    KSTART=3
	    KEND=4
	    HOLD(1)=RMINS(I,1) 
	  ELSEIF((RECON3).AND.(.NOT. RECON1))THEN
	    KSTART=1
	    KEND=2
	    HOLD(2)=RMINS(I,3)
          ELSEIF((RECON1).AND.(RECON3))THEN
	    HOLD(1)=RMINS(I,1)
	    HOLD(2)=RMINS(I,3)
            KSTART=1
            KEND=0
          ENDIF
C
          IF(IPRT.GE.1) WRITE(KFILDO,*) ' KSTART=',KSTART,' KEND=', KEND
C
          DO 250 J=KSTART,KEND
C
            IF(NINT(RMNABCD(J)).EQ.IZMISS) THEN
	      RESFLD(I)=ZMISS
	      IF(IPRT.GE.1) WRITE(KFILDO,*) 'RMNABCD(',J,') IS MISSING'
            ENDIF
C
 250      CONTINUE
C
          IF(NINT(RESFLD(I)).EQ.IZMISS) GO TO 300
        ENDIF
C
C         STEP 6. SINCE SIMPLE CHECK WAS NOT SUFFICIENT,
C         CONTINUE ON WITH A SERIES OF IF TESTS TO DETERMINE WHETHER
C         THE MAX/MIN OCCURRED AT THE BEGINNING OR THE END OF THE
C         THE "WINDOW" PERIOD.
C
C         STEP 6A. DO THE FOLLOWING FOR MIN
C
C         1. FIRST CHECK THE A AND B PORTION 
C
        IF(CMIN) THEN
C
	  IF(.NOT. RECON1) THEN
C
            IF(IPRT.GE.2) WRITE(KFILDO,*) 'MIN CHECKING A VS B'
C
            IF(NINT(RMNABCD(1)).LT.NINT(RMNABCD(2))) THEN
              HOLD(1)=RMNABCD(2)
            ELSEIF(NINT(RMNABCD(1)).GT.NINT(RMNABCD(2))) THEN
              HOLD(1)=RMINS(I,1)
            ELSEIF(NINT(RMNABCD(1)).EQ.NINT(RMNABCD(2))) THEN
C
              IF(NINT(RMNABCD(1)).EQ.NINT(RMINS(I,1)).AND.
     *          NINT(RMNABCD(2)).EQ.NINT(RMINS(I,1))) THEN
                HOLD(1)=RMINS(I,1)
              ELSE
                CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,ITIME1(1),ITIME2(1),
     *                       1,ZMISS,IMSG)
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'NEED TO CHK 6 & 8PM OBS'
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR MIN A VS B '
C
                IF(IMSG(1).GE.1) THEN
                  IF(IPRT.GE.1) WRITE(KFILDO,*) '6 OR 8PM OBS ARE', 
     *            ' MISSING'
                  RESFLD(I)=ZMISS
                  GO TO 300
                ENDIF
C
                IF(NINT(HRLY(I,ITIME1(1))).LT.
     *            NINT(HRLY(I,ITIME2(1)))) THEN
                  HOLD(1)=RMNABCD(2)
                ELSEIF(NINT(HRLY(I,ITIME1(1))).GT.
     *            NINT(HRLY(I,ITIME2(1)))) THEN
                  HOLD(1)=RMINS(I,1)
                ELSEIF(NINT(HRLY(I,ITIME1(1))).EQ.
     *            NINT(HRLY(I,ITIME2(1)))) THEN
C
                  IF(NINT(RMXABCD(1)).NE.-IZMISS.AND.
     *               NINT(RMXABCD(2)).NE.-IZMISS) THEN
C
                    IF(NINT(RMXABCD(1)).LT.NINT(RMXABCD(2))) THEN
                      HOLD(1)=RMNABCD(2)
                    ELSE
                      HOLD(1)=RMINS(I,1)
                    ENDIF
C
                  ELSE
                    HOLD(1)=RMINS(I,1)
                  ENDIF
C
                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
C            2. THEN CHECK THE C AND D PORTION 
C
          IF(.NOT. RECON3) THEN
            IF(IPRT.GE.2) WRITE(KFILDO,*) 'MIN CHECKING C VS D'
C
            IF(NINT(RMNABCD(3)).LT.NINT(RMNABCD(4))) THEN
              HOLD(2)=RMINS(I,3)
            ELSEIF(NINT(RMNABCD(3)).GT.NINT(RMNABCD(4))) THEN
              HOLD(2)=RMNABCD(3)
            ELSEIF(NINT(RMNABCD(3)).EQ.NINT(RMNABCD(4))) THEN
C
              IF(NINT(RMNABCD(3)).EQ.NINT(RMINS(I,3)).AND.
     *          NINT(RMNABCD(4)).EQ.NINT(RMINS(I,3))) THEN
                HOLD(2)=RMINS(I,3)
              ELSE
                CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,ITIME1(2),ITIME2(2),
     *                       1,ZMISS,IMSG)
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'NEED TO CHK 7 & 9AM OBS'
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR MIN C VS D '
C
                IF(IMSG(1).GE.1) THEN
                  IF(IPRT.GE.1) WRITE(KFILDO,*) '7 OR 9AM OBS ARE', 
     *            ' MISSING'
                  RESFLD(I)=ZMISS
                  GO TO 300
                ENDIF
C
                IF(NINT(HRLY(I,ITIME1(2))).LT.
     *            NINT(HRLY(I,ITIME2(2)))) THEN
                  HOLD(2)=RMINS(I,3)
                ELSEIF(NINT(HRLY(I,ITIME1(2))).GT.
     *            NINT(HRLY(I,ITIME2(2)))) THEN
                  HOLD(2)=RMNABCD(3)
                ELSEIF(NINT(HRLY(I,ITIME1(2))).EQ.
     *            NINT(HRLY(I,ITIME2(2)))) THEN
C
                  IF(NINT(RMXABCD(3)).NE.-IZMISS.AND.
     *              NINT(RMXABCD(4)).NE.-IZMISS) THEN
C
                    IF(NINT(RMXABCD(4)).LT.NINT(RMXABCD(3))) THEN
                      HOLD(2)=RMNABCD(3)
                    ELSE
                      HOLD(2)=RMINS(I,3)
                    ENDIF
C
                  ELSE
                    HOLD(2)=RMINS(I,3)
                  ENDIF
C
                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
          RESFLD(I)=MIN(HOLD(1),RMINS(I,2),HOLD(2))
C
          IF(IPRT.GE.1) WRITE(KFILDO,*) 'FOR STATION= ',CCALL(I),
     *                  'MIN FINALLY SET= ',RESFLD(I)
        ENDIF
C
C         STEP 6B. DO THE FOLLOWING FOR MAX
C
C         1. FIRST CHECK THE A AND B PORTION
C
	IF(CMAX) THEN
C
          IF(.NOT. RECON1) THEN
          IF(IPRT.GE.2) WRITE(KFILDO,*) 'MAX CHECKING A VS B'
C
            IF(NINT(RMXABCD(1)).GT.NINT(RMXABCD(2))) THEN
              HOLD(3)=RMXABCD(2)
            ELSEIF(NINT(RMXABCD(1)).LT.NINT(RMXABCD(2))) THEN
              HOLD(3)=RMAXS(I,1)
            ELSEIF(NINT(RMXABCD(1)).EQ.NINT(RMXABCD(2))) THEN
C
              IF(NINT(RMXABCD(1)).EQ.NINT(RMAXS(I,1)).AND.
     *          NINT(RMXABCD(2)).EQ.NINT(RMAXS(I,1))) THEN
                HOLD(3)=RMAXS(I,1)
              ELSE
                CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,ITIME1(3),ITIME2(3),
     *                       1,ZMISS,IMSG)
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'NEED TO CHK 6 & 8AM OBS'
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR MAX A VS B'
C
                IF(IMSG(1).GE.1) THEN
                  IF(IPRT.GE.1) WRITE(KFILDO,*) '6 OR 8AM OBS ARE', 
     *            ' MISSING'
                  RESFLD(I)=ZMISS
                  GO TO 300
                ENDIF
C
                IF(NINT(HRLY(I,ITIME1(3))).GT.
     *            NINT(HRLY(I,ITIME2(3)))) THEN
                  HOLD(3)=RMXABCD(2)
                ELSEIF(NINT(HRLY(I,ITIME1(3))).LT.
     *            NINT(HRLY(I,ITIME2(3)))) THEN
                  HOLD(3)=RMAXS(I,1)
                ELSEIF(NINT(HRLY(I,ITIME1(3))).EQ.
     *            NINT(HRLY(I,ITIME2(3)))) THEN
C
                  IF(NINT(RMNABCD(1)).NE.IZMISS.AND.
     *              NINT(RMNABCD(2)).NE.IZMISS) THEN
C
                    IF(NINT(RMNABCD(1)).GT.NINT(RMNABCD(2))) THEN
                      HOLD(3)=RMXABCD(2)
                    ELSE
                      HOLD(3)=RMAXS(I,1)
                    ENDIF
C
                  ELSE
                    HOLD(3)=RMAXS(I,1)
                  ENDIF
C
                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
C           2. THEN CHECK THE C AND D PORTION
C
          IF(.NOT. RECON3) THEN 
            IF(IPRT.GE.2) WRITE(KFILDO,*) 'MAX CHECKING C VS D'
C
            IF(NINT(RMXABCD(3)).GT.NINT(RMXABCD(4))) THEN
              HOLD(4)=RMAXS(I,3)
            ELSEIF(NINT(RMXABCD(3)).LT.NINT(RMXABCD(4))) THEN
              HOLD(4)=RMXABCD(3)
            ELSEIF(NINT(RMXABCD(3)).EQ.NINT(RMXABCD(4))) THEN
C
              IF(NINT(RMXABCD(3)).EQ.NINT(RMAXS(I,3)).AND.
     *          NINT(RMXABCD(4)).EQ.NINT(RMAXS(I,3))) THEN
                HOLD(4)=RMAXS(I,3)
              ELSE
                CALL CKFRMSG(ND1,MAXHRL,MXMN,I,HRLY,ITIME1(4),ITIME2(4),
     *                       1,ZMISS,IMSG)
                IF(IPRT.GE.2) WRITE(KFILDO,*) 'NEED TO CHK 6 & 8PM OBS'
	        IF(IPRT.GE.2) WRITE(KFILDO,*) 'FOR MAX C VS D'
C
                IF(IMSG(1).GE.1) THEN
                  IF(IPRT.GE.1) WRITE(KFILDO,*) '6 OR 8PM OBS ARE', 
     *            ' MISSING'
                  RESFLD(I)=ZMISS
                  GO TO 300
                ENDIF
C
                IF(NINT(HRLY(I,ITIME1(4))).GT.
     *            NINT(HRLY(I,ITIME2(4)))) THEN
                  HOLD(4)=RMAXS(I,3)
                ELSEIF(NINT(HRLY(I,ITIME1(4))).LT.
     *            NINT(HRLY(I,ITIME2(4)))) THEN
                  HOLD(4)=RMXABCD(3)
                ELSEIF(NINT(HRLY(I,ITIME1(4))).EQ.
     *            NINT(HRLY(I,ITIME2(4)))) THEN
C
                  IF(NINT(RMNABCD(3)).NE.IZMISS.AND.
     *              NINT(RMNABCD(4)).NE.IZMISS) THEN
C
                    IF(NINT(RMNABCD(4)).GT.NINT(RMNABCD(3))) THEN
                      HOLD(4)=RMXABCD(3)
                    ELSE
                      HOLD(4)=RMAXS(I,3)
                    ENDIF
C
                  ELSE
                    HOLD(4)=RMAXS(I,3)
                  ENDIF
C
                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
          RESFLD(I)=MAX(HOLD(3),RMAXS(I,2),HOLD(4))
C
          IF(IPRT.GE.1) WRITE(KFILDO,*)'FOR STATION= ',CCALL(I),
     *                  'MAX FINALLY SET= ',RESFLD(I)
        ENDIF
C
 300  CONTINUE
C
      RETURN
      END
C********************************************************
      SUBROUTINE CKABCD(RMNHRLY,RMXHRLY,RMNABCD,RMXABCD,ITZONE,
     *                  ZMISS,INUMZN)
C
C        AUGUST  1998   WEISS   TDL   MOS-2000
C        OCTOBER 2001   WEISS   IF STATEMENTS CONTAINING REAL VALUES
C                               CHANGED TO INTEGERS
C        MAY     2003   GLAHN   INSERTED WHITE SPACE
C
C        PURPOSE
C             FOR A GIVEN STATION, THIS ROUTINE WILL DETERMINE THE
C             WHETHER SEGMENTS A,B,C AND D CONTAIN A SUFFICIENT 
C             NUMBER OF HOURLY TEMPERATURE VALUES TO BE USED FOR 
C             DAYTIME/MAX AND NIGHTTIME/MIN ESTIMATES.  
C
C        VARIABLES
C
C          RMNHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (INPUT).
C          RMXHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MAX TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (INPUT).
C          RMNABCD(J) = HOLDS THE SEGMENTED MIN TEMPERATURES OF THE
C                       FIRST AND THIRD 6 HOUR TEMPERATURES A,B AND
C                       C,D RESPECTIVELY. B AND C ARE THE STARTING
C                       ENDING SEGMENTS OF THE NIGHTTIME MIN WINDOW
C                       FOR A GIVEN STATION, WHERE J=1,4 (INPUT/OUTPUT).
C          RMXABCD(J) = HOLDS THE SEGMENTED MAX TEMPERATURES OF THE
C                       FIRST AND THIRD 6 HOUR TEMPERATURES A,B AND
C                       C,D RESPECTIVELY. B AND C ARE THE STARTING
C                       ENDING SEGMENTS OF THE DAYTIME MAX WINDOW
C                       FOR A GIVEN STATION, WHERE J=1,4 (INPUT/OUTPUT).
C              ITZONE = ARRAY STATION'S TIME ZONE (INPUT).
C               ZMISS = MISSING VALUE OF 9999. (INPUT).
C              INUMZN = NUMBER OF HOURS (=14) USED TO SELECT HOURLY
C                       OBS REPRESENTING THE FIRST AND THIRD 6-HOUR
C                       TIME PERIODS (INPUT).
C
C        OTHER VARIABLES
C
C                IEND = END COUNTER OF NUMBER OF HOURS PER
C                       SEGMENT (INTERNAL).
C              ISEGHR = NUMBER OF HOURS PER SEGMENT (INTERNAL).
C              ISTART = START COUNTER OF NUMBER OF HOURS PER
C                       SEGMENT (INTERNAL).
C      MAX_COUNT(J,M) = FOR MAX ESTIMATES, THE TOTAL NUMBER OF HOURLY 
C                       OBS USED TO GENERATE SEGMENTS A,B,C, AND D FOR
C                       EACH TIME ZONE, WHERE M=1,7 AND J=1,4 
C                       (INTERNAL). NOTE: M=1,7 = ZONES -10 TO -4.
C             MAXH(L) = FOR MAX ESTIMATES, NUMBER OF HOURLY OBS
C                       PER SEGMENT (A - D) FOR ALL 7 TIME ZONES
C                       (WEST TO EAST) ,WHERE L=1,28 (INTERNAL).
C       MAX_MISSNG(J) = FOR MAX ESTIMATES, THE NUMBER OF MISSING 
C                       HOURLY OBS FOR SEGMENTS A,B,C, AND D FOR EACH
C                       STATION, WHERE J=1,4 (INTERNAL).
C      MIN_COUNT(J,M) = FOR MIN ESTIMATES, THE TOTAL NUMBER OF HOURLY
C                       OBS USED TO GENERATE SEGMENTS A,B,C, AND D FOR
C                       EACH TIME ZONE, WHERE M=1,7 AND J=1,4 
C             MINH(L) = FOR MIN ESTIMATES, NUMBER OF HOURLY OBS
C                       PER SEGMENT (A - D) FOR ALL 7 TIME ZONES
C                       (WEST TO EAST) ,WHERE L=1,28 (INTERNAL).
C                       (INTERNAL). NOTE: M=1,7 = ZONES -10 TO -4.
C       MIN_MISSNG(J) = FOR MIN ESTIMATES, THE NUMBER OF MISSING 
C                       HOURLY OBS FOR SEGMENTS A,B,C, AND D FOR EACH
C                       STATION, WHERE J=1,4 (INTERNAL).
C
      IMPLICIT NONE
C 
      INTEGER ITZONE,INUMZN
      INTEGER MAX_COUNT(4,7),MIN_COUNT(4,7),MAX_MISSNG(4),
     *        MIN_MISSNG(4),MAXH(28),MINH(28) 
      INTEGER ISTART,IEND,ISEGHR,I,II,J,N,NN,IZMISS
C
      REAL RMNABCD(4),RMXABCD(4),RMNHRLY(INUMZN),RMXHRLY(INUMZN),ZMISS
C
      DATA MINH/5,2,7,0,4,3,6,1,3,4,5,2,2,5,4,3,1,6,3,4,0,7,2,5,
     *          5,2,7,0/
      DATA MAXH/5,2,6,1,4,3,5,2,3,4,4,3,2,5,3,4,1,6,2,5,0,7,1,6,
     *          5,2,6,1/
C
C        STEP 1. INITIALIZE
C
      ISTART=0
      IEND=0
      ISEGHR=0
      IZMISS=9999
C
      DO 10 I=1,4
        MAX_MISSNG(I)=0
        MIN_MISSNG(I)=0
 10   CONTINUE
C
      II=0
C
      DO 15 I=1,7
C
        DO 12 J=1,4
          II=II+1
          MAX_COUNT(J,I)=MAXH(II)
          MIN_COUNT(J,I)=MINH(II)
 12     CONTINUE
C
 15   CONTINUE 
C
C        STEP 2A. FOR MAX, COUNT NUMBER OF MISSING PER SEGMENT
C
      ISTART=1
C
      DO 30 N=1,4
        ISEGHR=MAX_COUNT(N,ITZONE+11)
C
        IF(ISEGHR.GT.0) THEN
          IEND=ISTART+ISEGHR-1
C
          DO 25 NN=ISTART,IEND
C
            IF(NINT(RMXHRLY(NN)).EQ.-IZMISS) THEN
              MAX_MISSNG(N)=MAX_MISSNG(N)+1
            ENDIF
C
 25       CONTINUE
C
          ISTART=IEND+1
          IF((ISEGHR.EQ.1).AND.(MAX_MISSNG(N).EQ.1))RMXABCD(N)=-9988.
          IF((ISEGHR.GT.1).AND.(MAX_MISSNG(N).GE.2))RMXABCD(N)=-ZMISS
        ENDIF
C
 30   CONTINUE
C
C        STEP 2B. FOR MIN, COUNT NUMBER OF MISSING PER SEGMENT
C
      ISTART=1
C
      DO 40 N=1,4
        ISEGHR=MIN_COUNT(N,ITZONE+11)
C
        IF(ISEGHR.GT.0) THEN
          IEND=ISTART+ISEGHR-1
C
          DO 35 NN=ISTART,IEND
            IF(NINT(RMNHRLY(NN)).EQ.IZMISS) THEN
              MIN_MISSNG(N)=MIN_MISSNG(N)+1
            ENDIF
 35       CONTINUE
C
          ISTART=IEND+1
          IF((ISEGHR.EQ.1).AND.(MIN_MISSNG(N).EQ.1))RMNABCD(N)=9988.
          IF((ISEGHR.GT.1).AND.(MIN_MISSNG(N).GE.2))RMNABCD(N)=ZMISS
        ENDIF
C
 40   CONTINUE
C
      RETURN    
      END
C*******************************************************************
      SUBROUTINE CKFRMSG(ND1,MAXHRL,MXMN,ISTA,HRLY,JSTART,JEND,IWHICH,
     *                   ZMISS,IMSG)
C
C        AUGUST  1997    SETTELMAIER   TDL   MOS-2000
C        AUGUST  1998    WEISS         MODIFIED FOR REVISED MAXMIN SUBROUTINE
C        OCTOBER 2001    WEISS         IF STATEMENTS CONTAINING REAL VALUES
C                                      CHANGED TO INTEGERS.
C        MAY     2003    GLAHN         INSERTED WHITE SPACE
C
C        PURPOSE
C            THIS SUBROUTINE CHECKS FOR MISSING VALUES AMONGST ELEMENTS
C            OF THE ARRAY HRLY.
C
C        VARIABLES
C
C                 ND1 = MAXIMUM NUMBER OF STATIONS (INPUT).
C              MAXHRL = MAXIMUM NUMBER OF HOURS (=25) USED TO STORE
C                       THE HOURLY TEMPERATURES FROM ACROSS ALL
C                       POTENTIAL TIME ZONES (INPUT).
C                MXMN = FLAG TO INDICATE DAYTIME MAX (=2) OR
C                       NIGHTTIME MIN (=1) PROCESSING (INPUT).
C                ISTA = STATION COUNTER FOR ARRAY HRLY (INPUT). 
C           HRLY(N,J) = THE HOURLY TEMPERATURE OBS FOR EACH STATION
C                       WHERE N=1,ND1 AND J=1,MAXHRL (INTERNAL).
C              JSTART = THE STARTING POINT OF THE RANGE OF HOURLY OBS
C                       MAKING UP THE TIME PERIOD OF A 6-HR MIN/MAX
C                       ESTIMATE (INPUT).
C                JEND = THE ENDING POINT OF THE RANGE OF HOURLY OBS
C                       MAKING UP THE TIME PERIOD OF A 6-HR MIN/MAX
C                       ESTIMATE (INPUT).
C              IWHICH = FLAG FOR 6 HOUR MAX/MIN MISSING CHECK (=6) OR
C                       FLAG FOR 1 HOUR MAX/MIN MISSING CHECK (=1),
C                       (INPUT).
C               ZMISS = MISSING VALUE OF ZMISS (INPUT).
C             IMSG(K) = ARRAY TO HOLD THE NUMBER OF MISSING IN A 
C                       THE HRLY ARRAY, WHERE K=1,3 (INPUT).
C
C        OTHER VARIABLES
C            
C               JFROM = STARTING HOUR FOR A 6 HOUR PERIOD (INTERNAL).
C              JFROM2 = ENDING HOUR FOR A 6 HOUR PERIOD (INTERNAL).
C
      IMPLICIT NONE 
C
      INTEGER ND1,MAXHRL,MXMN,ISTA,JSTART,JEND,IWHICH,IMSG(3)
      INTEGER JFROM,JFROM2
      INTEGER I,J,K,IZMISS
C
      REAL HRLY(ND1,MAXHRL),ZMISS
C
C        STEP 1. INITIALIZE THE MISSING VALUE COUNTER ARRAY 
C
      IZMISS=9999
C
      DO 2 J=1,3
        IMSG(J)=0
  2   CONTINUE
C
C        STEP 2A. FOR IWHICH=6, CHECK A RANGE OF HOURLY OBSERVATIONS
C                 FOR MISSING VALUES OVER A 6-HR PERIOD. THE PERIOD 
C                 IS DETERMINED BY JSTART.
C
      IF(IWHICH.EQ.6) THEN
C
        DO 7 K=1,3
          JFROM=JSTART+(IWHICH*K)-IWHICH
          JFROM2=JSTART+(IWHICH*K)
C
          DO 3 I=JFROM,JFROM2
            IF(NINT(HRLY(ISTA,I)).EQ.IZMISS) THEN
              IMSG(K)=IMSG(K)+1
              IF(MXMN.EQ.2) HRLY(ISTA,I)=-ZMISS
            ENDIF
C
  3       CONTINUE
C
  7     CONTINUE
C
      ENDIF
C
C        STEP 2B. FOR IWHICH=1, CHECK TWO HOURLY OBSERVATIONS
C                 FOR MISSING VALUES. THE PERIOD IS JUST TWO 
C                 SPECIFIC HOURS DETERMINED BY THE JSTART AND JEND.
C
      IF(IWHICH.EQ.1) THEN
        K=1
C
        IF(NINT(HRLY(ISTA,JSTART)).EQ.IZMISS.OR.
     *    NINT(HRLY(ISTA,JSTART)).EQ.-IZMISS) THEN
          IMSG(K)=IMSG(K)+1
        ENDIF
C
        IF(NINT(HRLY(ISTA,JEND)).EQ.IZMISS.OR.
     *    NINT(HRLY(ISTA,JEND)).EQ.-IZMISS) THEN
          IMSG(K)=IMSG(K)+1
        ENDIF
C
      ENDIF
C
      RETURN
      END
C*******************************************************************
      SUBROUTINE RMXMN_FILL(INUMZN,ISTA,IZ,MAXHRL,ND1,
     *                       HRLY,RMXHRLY,RMNHRLY,ZMISS,
     *                       CMAX,CMIN,RFILL)
C
C
C        DECEMBER 2000   WEISS   MDL   MOS-2000
C        OCTOBER  2001   WEISS   IF STATEMENTS CONTAINING REAL VALUES
C                                CHANGED TO INTEGERS
C
C        PURPOSE
C            THIS SUBROUTINE FILLS THE ARRAYS RMXHRLY AND RMNHRLY WITH
C            HOURLY TEMPERATURE OBS. THESE ARRAY ARE THEN USED TO 
C            GENERATE SEGMENTS A,B,C OR D FOR EITHER MAX OR MIN.
C
C        VARIABLES
C
C              INUMZN = NUMBER OF HOURS (=14) USED TO SELECT HOURLY
C                       OBS REPRESENTING THE FIRST AND THIRD 6-HOUR
C                       TIME PERIODS (INPUT).
C                ISTA = STATION COUNTER FOR ARRAY HRLY (INPUT).
C                  IZ = DETERMINES STARTING POINT IN HRLY ARRAY FOR
C                       RECONSTRUCTING 6 HOUR MAX/MIN VALUES BASED ON
C                       HOURLY OBS, DEPENDING ON TIME ZONE (INPUT).
C              MAXHRL = MAXIMUM NUMBER OF HOURS (=25) USED TO STORE
C                       THE HOURLY TEMPERATURES FROM ACROSS ALL
C                       POTENTIAL TIME ZONES (INPUT).
C                 ND1 = MAXIMUM NUMBER OF STATIONS (INPUT).
C           HRLY(N,J) = THE HOURLY TEMPERATURE OBS FOR EACH STATION
C                       WHERE N=1,ND1 AND J=1,MAXHRL (INPUT).
C          RMNHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (OUTPUT).
C          RMXHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (OUTPUT).
C               ZMISS = MISSING VALUE OF 9999. (INPUT).
C                CMAX = LOGICAL PARAMETER FOR DAYTIME MAXIMUM
C                       PROCESSING (INPUT).
C                CMIN = LOGICAL PARAMETER FOR NIGHTTIME MINIMUM
C                       PROCESSING (INPUT).
C               RFILL = LOGICAL PARAMETER INDICATING THAT SUBROUTINE
C                       RMXMN_FILL HAS BEEN CALLED AT LEAST ONCE 
C                       (ALL THAT IS REQUIRED) (INPUT/OUTPUT).
C
      IMPLICIT NONE
C
      LOGICAL CMAX,CMIN,RFILL
C
      INTEGER INUMZN,ISTA,IZ,MAXHRL,ND1
      INTEGER II,J,IZMISS
C
      REAL HRLY(ND1,MAXHRL)
      REAL RMNHRLY(INUMZN),RMXHRLY(INUMZN),ZMISS
C
C        STEP 1. INITIALIZE THE RMXHRLY AND RMNHRLY ARRAYS
C
      IZMISS=9999
C
      DO 2 II=1,INUMZN
          RMXHRLY(II)=-ZMISS
          RMNHRLY(II)=ZMISS
 2    CONTINUE
C
      RFILL=.FALSE.
C
C        STEP 2. FILL ARRAYS RMXHRLY AND RMNHRLY WITH HRLY OBS
C
      IF(CMAX) THEN 
C
        DO 220 J=1,7 
          RMXHRLY(J)=HRLY(ISTA,J+6-IZ) 
          RMNHRLY(J)=HRLY(ISTA,J+6-IZ)
          IF(NINT(HRLY(ISTA,J+6-IZ)).EQ.-IZMISS) RMNHRLY(J)=ZMISS
 220    CONTINUE 
C
        DO 225 J=8,14 
          RMXHRLY(J)=HRLY(ISTA,J+11-IZ) 
          RMNHRLY(J)=HRLY(ISTA,J+11-IZ)
          IF(NINT(HRLY(ISTA,J+11-IZ)).EQ.-IZMISS) RMNHRLY(J)=ZMISS
 225    CONTINUE
C
      ELSEIF(CMIN) THEN
C
        DO 230 J=1,7 
          RMXHRLY(J)=HRLY(ISTA,J+6-IZ)
          IF(NINT(HRLY(ISTA,J+6-IZ)).EQ.IZMISS) RMXHRLY(J)=-ZMISS
          RMNHRLY(J)=HRLY(ISTA,J+6-IZ)
 230    CONTINUE 
C
        DO 235 J=8,14
          RMXHRLY(J)=HRLY(ISTA,J+11-IZ)
          IF(NINT(HRLY(ISTA,J+11-IZ)).EQ.IZMISS) RMXHRLY(J)=-ZMISS
          RMNHRLY(J)=HRLY(ISTA,J+11-IZ) 
 235    CONTINUE
C
      ENDIF
C
      RETURN
      END
C*******************************************************************
      SUBROUTINE MAXMIN_PART(INUMZN,ITZON,JCOUNT,
     *                       RMXHRLY,RMNHRLY,RMXMN13,ZMISS,
     *                       CMAX,CMIN,RECON1,RECON3,NUM_MISS,
     *                       IPRT,KFILDO)
C
C        DECEMBER 2000   WEISS   MDL   MOS-2000
C        OCTOBER  2001   WEISS   IF STATEMENTS CONTAINING REAL VALUES
C                                CHANGED TO INTEGERS.
C
C        PURPOSE
C            THIS SUBROUTINE CALCULATES (RECONSTRUCTS) THE 
C            MAX1/MIN1 AND MAX3/MIN3 VALUES FOR PART-TIME
C            STATIONS. THE MAX1/MIN1 ARE EQUIVALENT TO SEGMENT B
C            AND MAX3/MIN3 ARE EQUIVALENT TO SEGMENT C. THE CHECK
C            OF MISSING HOURLY OBS FOR THE RECONSTRUCTED MAX1/MIN1
C            AND MAX3/MIN3 VALUES IS ALSO CONDUCTED IN THIS SUBROUTINE.
C
C        VARIABLES
C
C              INUMZN = NUMBER OF HOURS (=14) USED TO SELECT HOURLY
C                       OBS REPRESENTING THE FIRST AND THIRD 6-HOUR
C                       TIME PERIODS (INPUT).
C               ITZON = STATION'S TIME ZONE (INPUT).
C              JCOUNT = COUNTER OF THE THREE MAX'S OR MIN'S (INPUT).  
C          RMNHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (INPUT).
C          RMXHRLY(K) = HOLDS THE HOURLY TEMPERATURES FOR THE FIRST
C                       AND THIRD SIX HOUR PERIODS AND IS USED TO
C                       GENERATE MIN TEMPERATURE SEGMENTS A,B,C AND D,
C                       WHERE K=1,14 (INPUT).
C             RMXMN13 = THE RECONSTRUCTED MAX1/MIN1 AND/OR MAX3/MIN3
C                       VALUE RETURNED BY SUBROUTINE MAXMIN_PART, 
C                       FOR PART-TIME STATION MAX/MIN ESTIMATES 
C                       (OUTPUT).
C               ZMISS = MISSING VALUE OF 9999. (INPUT).
C                CMAX = LOGICAL PARAMETER FOR DAYTIME MAXIMUM
C                       PROCESSING (INPUT).
C                CMIN = LOGICAL PARAMETER FOR NIGHTTIME MINIMUM
C                       PROCESSING (INPUT).
C              RECON1 = LOGICAL PARAMETER USED TO INDICATE EITHER
C                       MAX1/MIN1 HAS BEEN GENERATED FOR A 
C                       PART-TIME STATION (INPUT/OUTPUT).
C              RECON3 = LOGICAL PARAMETER USED TO INDICATE EITHER
C                       MAX3/MIN3 HAS BEEN GENERATED FOR A 
C                       PART-TIME STATION (INPUT/OUTPUT).
C                IPRT = PRINT FLAG USED TO (IPRT > 0), PRINT OUT
C                       INCREASING AMOUNTS OF DIAGNOSTICS
C                       1=SOME(RECOMMENDED), 2=ALL DIAGNOSTICS (INPUT).
C              KFILDO = DEFAULT UNIT NUMBER FOR OUTPUT (PRINT) FILE
C                       (INPUT).
C
C        OTHER VARIABLES
C                IEND = END COUNTER OF NUMBER OF HOURS PER
C                       SEGMENT (INTERNAL).
C              ISEGHR = NUMBER OF HOURS PER SEGMENT (INTERNAL).
C              ISTART = START COUNTER OF NUMBER OF HOURS PER
C                       SEGMENT (INTERNAL).
C      MAX_COUNT(J,M) = FOR MAX ESTIMATES, THE TOTAL NUMBER OF HOURLY
C                       OBS USED TO GENERATE SEGMENTS B AND C FOR
C                       EACH TIME ZONE, WHERE M=1,7 AND J=1,2
C             MAXH(L) = FOR MAX ESTIMATES, NUMBER OF HOURLY OBS
C                       PER SEGMENT (B - C) FOR ALL 7 TIME ZONES
C                       (WEST TO EAST) ,WHERE L=1,14 (INTERNAL).
C            MAXST(L) = FOR MAX ESTIMATES, THE STARTING POINT FOR
C                       TESTING MISSING VALUES OF SEGMENTS B AND C
C                       WHERE L=1,7 (INTERNAL).
C      MIN_COUNT(J,M) = FOR MIN ESTIMATES, THE TOTAL NUMBER OF HOURLY
C                       OBS USED TO GENERATE SEGMENTS B AND C FOR
C                       EACH TIME ZONE, WHERE M=1,7 AND J=1,2
C             MINH(L) = FOR MIN ESTIMATES, NUMBER OF HOURLY OBS
C                       PER SEGMENT (B - C) FOR ALL 7 TIME ZONES
C                       (WEST TO EAST) ,WHERE L=1,14 (INTERNAL).
C            MINST(L) = FOR MIN ESTIMATES, THE STARTING POINT FOR
C                       TESTING MISSING VALUES OF SEGMENTS B AND C
C                       WHERE L=1,7 (INTERNAL).
C            NUM_MISS = FOR MAX AND MIN ESTIMATES, THE NUMBER OF
C                       MISSING HOURLY OBS FOR SEGMENTS B AND C
C                       FOR A GIVEN STATION (INTERNAL/OUTPUT).
C              RMAX1B = THE RECONSTRUCTED MAX1 VALUE (INTERNAL).
C              RMAX3C = THE RECONSTRUCTED MAX3 VALUE (INTERNAL).
C              RMIN1B = THE RECONSTRUCTED MIN1 VALUE (INTERNAL).
C              RMIN3C = THE RECONSTRUCTED MIN3 VALUE (INTERNAL).
C
      IMPLICIT NONE
C
      LOGICAL CMAX,CMIN,RECON1,RECON3
C
      INTEGER JCOUNT,INUMZN,ITZON
      INTEGER ISEGHR,ISTART,IEND,NUM_MISS,IPRT,KFILDO
      INTEGER MAXH(14),MAXST(7),MINH(14),MINST(7)
      INTEGER MAX_COUNT(2,7),MIN_COUNT(2,7)
      INTEGER I,II,J,JJ,IZMISS
C
      REAL RMNHRLY(INUMZN),RMXHRLY(INUMZN),ZMISS
      REAL RMAX1B,RMAX3C,RMIN1B,RMIN3C,RMXMN13
C
      DATA MAXH/2,6,3,5,4,4,5,3,6,2,7,1,2,6/
      DATA MAXST/5,4,3,2,1,0,5/
      DATA MINH/2,7,3,6,4,5,5,4,6,3,7,2,2,7/
      DATA MINST/5,4,3,2,1,0,5/
C
C        STEP 1. INITIALIZE MAX1B, MAX3C, MIN1B, MIN3C.
C
          RMAX1B=-ZMISS
	  RMAX3C=-ZMISS
          RMIN1B=ZMISS
	  RMIN3C=ZMISS
          NUM_MISS=0
          RMXMN13=ZMISS
          IZMISS=9999
C
C*********************************************************
C
C        STEP 2A.SET UP THE THE MAX VALUES SEGMENT B
C                AND C FOR EACH TIME ZONE.
C                  MAX WINDOW = 7AM - 7PM
C
      IF(CMAX) THEN
C
        IF(JCOUNT.EQ.1) THEN
C
          SELECT CASE (ITZON)
            CASE (-10)
C               1. ALEUTIANS/HAWAII
              RMAX1B=MAX(RMXHRLY(6),RMXHRLY(7))
            CASE (-9)
C               2. ALASKA
              RMAX1B=MAX(RMXHRLY(5),RMXHRLY(6),RMXHRLY(7))
            CASE (-8)          
C               3. PACIFIC
              RMAX1B=MAX(RMXHRLY(4),RMXHRLY(5),RMXHRLY(6),
     *                   RMXHRLY(7))
            CASE (-7)
C               4. MOUNTAIN
              RMAX1B=MAX(RMXHRLY(3),RMXHRLY(4),RMXHRLY(5),
     *                   RMXHRLY(6),RMXHRLY(7))
            CASE (-6)
C               5. CENTRAL
              RMAX1B=MAX(RMXHRLY(2),RMXHRLY(3),RMXHRLY(4),RMXHRLY(5),
     *                   RMXHRLY(6),RMXHRLY(7))
            CASE (-5)
C               6. EASTERN
              RMAX1B=MAX(RMXHRLY(1),RMXHRLY(2),RMXHRLY(3),RMXHRLY(4),
     *                   RMXHRLY(5),RMXHRLY(6),RMXHRLY(7))
            CASE (-4)
C               7. ATLANTIC
              RMAX1B=MAX(RMXHRLY(6),RMXHRLY(7))
          END SELECT
C
          RECON1=.TRUE.
          RMXMN13=RMAX1B
        ELSEIF(JCOUNT.EQ.3) THEN
C
          SELECT CASE (ITZON)
          CASE(-10) 
C             1. ALEUTIANS/HAWAII
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                 RMXHRLY(11),RMXHRLY(12),RMXHRLY(13))
          CASE(-9)
C              2. ALASKA
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                 RMXHRLY(11),RMXHRLY(12))
          CASE(-8)
C              3. PACIFIC
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),
     *                 RMXHRLY(11))
          CASE(-7)
C              4. MOUNTAIN
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10))
          CASE(-6)
C              5. CENTRAL
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9))
          CASE(-5)
C              6. EASTERN
            RMAX3C=RMXHRLY(8)
          CASE(-4)
C              7. ATLANTIC
            RMAX3C=MAX(RMXHRLY(8),RMXHRLY(9),RMXHRLY(10),RMXHRLY(11),
     *                 RMXHRLY(12),RMXHRLY(13))
          END SELECT
C
          RECON3=.TRUE.
          RMXMN13=RMAX3C
        ENDIF
C
      ENDIF
C
C        STEP 2B. SET UP THE THE MIN VALUES SEGMENT B
C                 AND C FOR EACH TIME ZONE
C                 MIN WINDOW = 7PM - 8AM
C
      IF(CMIN) THEN
C
        IF(JCOUNT.EQ.1) THEN
C
          SELECT CASE (ITZON)
          CASE (-10)
C              1. ALEUTIANS/HAWAII
            RMIN1B=MIN(RMNHRLY(6),RMNHRLY(7))
          CASE (-9)
C              2. ALASKA
            RMIN1B=MIN(RMNHRLY(5),RMNHRLY(6),RMNHRLY(7))
          CASE (-8)
C              3. PACIFIC
            RMIN1B=MIN(RMNHRLY(4),RMNHRLY(5),RMNHRLY(6),
     *                 RMNHRLY(7))
          CASE (-7)
C              4. MOUNTAIN
            RMIN1B=MIN(RMNHRLY(3),RMNHRLY(4),RMNHRLY(5),
     *                 RMNHRLY(6),RMNHRLY(7))
          CASE (-6)
C              5. CENTRAL
            RMIN1B=MIN(RMNHRLY(2),RMNHRLY(3),RMNHRLY(4),RMNHRLY(5),
     *                 RMNHRLY(6),RMNHRLY(7))
          CASE (-5)
C              6. EASTERN
            RMIN1B=MIN(RMNHRLY(1),RMNHRLY(2),RMNHRLY(3),RMNHRLY(4),
     *                 RMNHRLY(5),RMNHRLY(6),RMNHRLY(7))
          CASE (-4)
C              7. ATLANTIC
            RMIN1B=MIN(RMNHRLY(6),RMNHRLY(7))
          END SELECT
C
          RECON1=.TRUE.
          RMXMN13=RMIN1B
        ELSEIF(JCOUNT.EQ.3) THEN
C
          SELECT CASE (ITZON)
          CASE (-10)
C             1. ALEUTIANS/HAWAII
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),RMNHRLY(11),
     *                 RMNHRLY(12),RMNHRLY(13),RMNHRLY(14))
          CASE (-9)
C             2. ALASKA
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *                 RMNHRLY(11),RMNHRLY(12),RMNHRLY(13))
          CASE (-8)
C             3. PACIFIC
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *               RMNHRLY(11),RMNHRLY(12))
          CASE (-7)
C             4. MOUNTAIN
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),
     *                 RMNHRLY(11))
          CASE (-6)
C             5. CENTRAL
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10))
          CASE (-5)
C             6. EASTERN
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9))
          CASE (-4)
C             7. ATLANTIC
            RMIN3C=MIN(RMNHRLY(8),RMNHRLY(9),RMNHRLY(10),RMNHRLY(11),
     *                 RMNHRLY(12),RMNHRLY(13),RMNHRLY(14))
          END SELECT
C
          RECON3=.TRUE.
          RMXMN13=RMIN3C
        ENDIF
C
      ENDIF
C
C*****************************************************************
C
C        STEP 3. CHECK THE RECONSTRUCTED MAX1/MIN1 AND MAX3/MIN3
C                FOR MISSING DATA. THE CHECK WILL BE CONSISTENT
C                THE CHECKS FOUND IN SUBROUTINE CKABCD.
C
      II=0
C
      DO 200 I=1,7
C
        DO 190 J=1,2
          II=II+1
          IF(CMAX) MAX_COUNT(J,I)=MAXH(II)
          IF(CMIN) MIN_COUNT(J,I)=MINH(II)
 190    CONTINUE
C
 200    CONTINUE
C
C
C        STEP 3A. FOR MAX, COUNT NUMBER OF MISSING PER SEGMENT
C
      IF(CMAX) THEN
C
        IF(RECON1) THEN
          ISTART=MAXST(ITZON+11)+1
          J=1
        ENDIF
C
        IF(RECON3) THEN
          ISTART=MAXST(ITZON+11)+MAX_COUNT(1,ITZON+11)+1
          J=2
        ENDIF
C
        ISEGHR=MAX_COUNT(J,ITZON+11)
C
        IF(ISEGHR.GT.0) THEN
          IEND=ISTART+ISEGHR-1
            IF(IPRT.GE.1) WRITE(KFILDO,*) 'PART-TIME MAX'
C
          DO 250 JJ=ISTART,IEND
            IF(NINT(RMXHRLY(JJ)).EQ.-IZMISS) NUM_MISS=NUM_MISS+1
            IF(IPRT.GE.1) WRITE(KFILDO,*) 'MAX HOURLY=',RMXHRLY(JJ),
     *      ' ISEGHR=',ISEGHR,' NUM_MISS=',NUM_MISS
 250      CONTINUE
C
          ISTART=IEND+1
          IF((ISEGHR.EQ.1).AND.(NUM_MISS.EQ.1))RMXMN13=-9998.
          IF((ISEGHR.GT.1).AND.(NUM_MISS.GE.2))RMXMN13=-ZMISS
        ENDIF
CCC
        IF((RECON1).AND.(IPRT.GE.1)) THEN
          IF(.NOT. RECON3) WRITE(KFILDO,*) 'RECONSTRUCTED',
     *    ' MAX1=',RMXMN13
        ENDIF
C
        IF((RECON3).AND.(IPRT.GE.1)) WRITE(KFILDO,*) 'RECONSTRUCTED',
     *  ' MAX3=',RMXMN13
CCC
      ENDIF
C
C        STEP 3B. FOR MIN, COUNT NUMBER OF MISSING PER SEGMENT
C
      IF(CMIN) THEN
C
        IF(RECON1) THEN
          ISTART=MINST(ITZON+11)+1
          J=1
        ENDIF
C
        IF(RECON3) THEN
          ISTART=MINST(ITZON+11)+MIN_COUNT(1,ITZON+11)+1
          J=2
        ENDIF
C
        ISEGHR=MIN_COUNT(J,ITZON+11)
C
        IF(ISEGHR.GT.0) THEN
          IEND=ISTART+ISEGHR-1
	    IF(IPRT.GE.1) WRITE(KFILDO,*) 'PART-TIME MIN'
C
          DO 260 JJ=ISTART,IEND
            IF(NINT(RMNHRLY(JJ)).EQ.IZMISS) NUM_MISS=NUM_MISS+1
	    IF(IPRT.GE.1) WRITE(KFILDO,*) 'MIN HOURLY=',RMXHRLY(JJ),
     *      ' ISEGHR=',ISEGHR,' NUM_MISS=',NUM_MISS
 260      CONTINUE
C
          ISTART=IEND+1
          IF((ISEGHR.EQ.1).AND.(NUM_MISS.EQ.1))RMXMN13=9998.
          IF((ISEGHR.GT.1).AND.(NUM_MISS.GE.2))RMXMN13=ZMISS
        ENDIF
C
        IF((RECON1).AND.(IPRT.GE.1)) THEN
          IF(.NOT. RECON3) WRITE(KFILDO,*) 'RECONSTRUCTED',
     *    ' MIN1=',RMXMN13
        ENDIF
C
        IF((RECON3).AND.(IPRT.GE.1)) WRITE(KFILDO,*) 'RECONSTRUCTED',
     *  ' MIN3=',RMXMN13
C
      ENDIF
C
      RETURN
      END
