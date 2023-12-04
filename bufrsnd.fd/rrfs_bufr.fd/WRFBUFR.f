      PROGRAM WRFBUFR
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C MAIN PROGRAM: WRFBUFR
C   PRGMMR: PYLE             ORG: EMC/MMB    DATE: 2004-11-25
C     
C ABSTRACT:  
C     THIS PROGRAM DRIVES THE EXTERNAL WRF BUFR POST PROCESSOR.
C     
C PROGRAM HISTORY LOG (FOR THE PROF CODE):
C   99-04-22  T BLACK - ORIGINATOR
C   02-07-01  G MANIKIN - FIXED PROBLEM WITH DHCNVC AND DHRAIN
C                          COMPUTATIONS - SEE COMMENTS BELOW
C   03-04-01  M PYLE - BEGAN CONVERTING FOR WRF
C   04-05-26  M PYLE - MADE CHANGES FOR WRF-NMM
C   04-11-24  M PYLE - ELIMINATED USE OF PARMETA FILE, DIMENSIONS
C                      NOW READ IN FROM WRF OUTPUT FILE, WITH WORKING
C                      ARRAYS ALLOCATED TO NEEDED DIMENSIONS.  UNIFIED
C                      WRF-EM AND WRF-NMM VERSIONS INTO A SINGLE CODE
C                      THAT READS EITHER NETCDF OR BINARY OUTPUT FROM
C                      THE WRF MODEL.
C   05-08-29  M PYLE - ELIMINATE THE NEED TO RETAIN ALL WRF HISTORY FILES.
C   07-08-06  J Du & B Zhou - A new prefilename was defined to correctly
C                      calculate precip rate during INCHOUR interval                   
C     
C USAGE:    WRFPOST
C   INPUT ARGUMENT LIST:
C     NONE     
C
C   OUTPUT ARGUMENT LIST: 
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       NONE
C     LIBRARY:
C       COMMON - CTLBLK
C                RQSTFLD
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : IBM RS/6000 SP
C$$$  
C
C
C
C     INCLUDE ARRAY DIMENSIONS.
C      INCLUDE "parmeta"
C      INCLUDE "mpif.h"
C
C     DECLARE VARIABLES.
C     
C     
C     INCLUDE COMMON BLOCKS.
!tst      INCLUDE "CTLBLK.comm"
C     
C     SET HEADER WRITER FLAGS TO TRUE.
c
C
      real rinc(5)
      integer jdate(8),idate(8)
          integer iii
      character(len=6) :: IOFORM,model
      character(len=98) :: newname
      character(len=256) :: fileName,fileNamedyn,fileNamephys
      character(len=256) :: prefileName,prefileNamedyn,
     &                      prefileNamephys
      character(len=19) :: DateStr
      integer :: DataHandle, IHR, INCR, NSTAT
!      integer, parameter:: INCR=3

C
C**************************************************************************
C
C     START PROGRAM WRFBUFR.
C
!	write(0,*) 'to read statements'
       read(11,111) fileNameDyn
       read(11,111) fileNamePhys
!	write(0,*) 'initial filename= ', filename
       read(11,113) model
!	write(0,*) 'model type= ', model
       read(11,113) IOFORM
!	write(0,*) 'ioform= ', ioform
       read(11,112) DateStr
!	write(0,*) 'datestr= ', datestr
       read(11,*) NFILES
       read(11,*) INCR
       read(11,*) IHR
       read(11,*) NSTAT
       read(11,111) prefileNameDyn
       read(11,111) prefileNamephys
!        write(0,*) 'previous filename= ', prefilename

!!!! CHANGE THIS ASSUMPTION???

! assume for now that the first date in the stdin file is the start date

       read(DateStr,300) iyear,imn,iday,ihrst

C      write(*,*) 'in WRFPOST iyear,imn,iday,ihrst',iyear,imn,iday,ihrst
 300  format(i4,1x,i2,1x,i2,1x,i2)
	
	IDATE=0

         IDATE(2)=imn
         IDATE(3)=iday
         IDATE(1)=iyear
         IDATE(5)=ihrst

 111  format(a256)
 112  format(a19)
 113  format(a6)
C

	do N=1,NFILES

	len=index(filename,' ')-1

!	IHR=(N-1)*INCR

!	add forecast hour to start time
	
	RINC(1)=0.
	RINC(2)=float(IHR)
	RINC(3)=0.
	RINC(4)=0.
	RINC(5)=0.

	call w3movdat(rinc,idate,jdate)

	if (model(1:4) .eq. 'NCEP') then
	write(DateStr,302) JDATE(1),JDATE(2),JDATE(3),JDATE(5)
	elseif (model(1:4) .eq. 'NCAR') then
	write(DateStr,302) JDATE(1),JDATE(2),JDATE(3),JDATE(5)
	endif

c20080707	filename=filename(1:len-19)//DateStr

!	if (model(1:4) .eq. 'NCAR') then
!	filename=filename(1:len-19)//DateStr
!	filename(len-2:len-2)='_'	
!	filename(len-5:len-5)='_'	
!	endif

 301  format(i4,'-',i2.2,'-',i2.2,'T',i2.2,':00:00')
 302  format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':00:00')

!	write(0,*) 'calling prof.... '
!	write(0,*) 'datestr: ', datestr
!	write(0,*) 'fileName ', fileName
!	write(0,*) 'IHR: ', IHR
!	write(0,*) '--------------------------------'

	if (ioform(1:6) .eq. 'netcdf') then 

	if (model(1:4) .eq. 'FV3S') then
        CALL PROF_FV3SAR_NET(fileNamedyn,filenamephys,
     &                       DateStr,IHR,INCR,NSTAT)
	endif

	endif

	write(0,*) 'back from prof'
	END DO
	write(0,*) 'back end do, next line is STOP0'
      STOP0
      END

