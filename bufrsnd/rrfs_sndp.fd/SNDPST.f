      PROGRAM SNDPST
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
C   00--5-15  ERIC ROGERS   - PUT NSOIL AND LM1 IN INCLUDED PARAMETER 
C                             FILE
C   04-11-18  BRAD FERRIER  - Added Cu precip rate; separated cloud water,
C             GEOFF MANIKIN   rain, & ice from CWM array => feed into visibility
C              
C   2017-03-30 Matthew Pyle - Creating a single main program that then will
C                              make calls to
C                             the previous SNDPST_EM.f type routines.
C
C   2022-07-12 Matthew Pyle - Shift to pure FV3 system without kludgy
C                             level/dynamical core blocks.
C
C
C     
C USAGE:   
C   INPUT ARGUMENT LIST:
C     NONE    
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       SNDPST_EM
C       SNDPST_NMM
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

      character(len=3):: core
      integer :: nlev, nstat,length

      read(5,*) nlev,nstat,length

      call SNDPST_FV3S(nlev,nstat,length)

      end program sndpst
