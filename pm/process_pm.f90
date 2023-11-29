program pmbufr

 implicit none

 integer, parameter :: mxmn=35, mxlv=255, maxcnt=5000, pm_limit=-15
 integer :: ireadmg,ireadsb,idate,nmsg,ntb,nsubset,nlvl
 integer :: i,j,k,NARG,iret_hdr,iret_ob,interval
 integer :: dhr,iret_hdr1,iret_ob1,cnt,ios
 integer :: unit_table,unit_out,unit_site,unit_var,nt,pm_len
 integer(8) :: z_sid
 character(80):: hdstr='SID XOB YOB DHR TYP T29 SQN PROCN RPT CAT TYPO TSIG'
 character(80):: obstr='TPHR QCIND COPOPM'
 character(8) :: subset,c_sid
 character(10):: analysis_time
 character(20):: str,infile,outfile
 character(20),allocatable,dimension(:)::sta_id, lat_str, lon_str, elv_str, pm_aqi_str, pm_measured_str, pm_str

 real*8,parameter:: bmiss=10e10_8
 real(8) :: hdr(mxmn),obs(mxmn,mxlv)
 real(8) :: rstation_id, lat, lon, elv, pm

 logical :: ifexist

 namelist/setup/analysis_time,infile,outfile,cnt

 ! only c_sid*8 is encoded to bufr file since SID is 64 bit in bufrtable
 equivalence(rstation_id,c_sid)

! read namelist
 inquire(file='namelist.pm', EXIST=ifexist )
 if(ifexist) then
   open(15, file='namelist.pm')
     read(15,setup)
   close(15)
   write(*,'(a45,a20,3x,2a20,i8)') 'analysis_time, infile, outfile, cnt: ' ,analysis_time, infile, outfile, cnt    
   write(*,*)
 else
   write(*,*) "Usage: analysis_time, infile, outfile, pm_count"
 stop 123

 endif

 if (cnt==0) then
   write(*,*) 'No observation available, exit'
   stop
 end if

! open and read pm2.5 ascii file
 unit_var=11
 allocate(sta_id(cnt))
 allocate(lat_str(cnt))
 allocate(lon_str(cnt))
 allocate(elv_str(cnt))
 allocate(pm_aqi_str(cnt))
 allocate(pm_measured_str(cnt))
 allocate(pm_str(cnt))

 open (unit = unit_var, file = infile,status = 'old')
 do i = 1,cnt
   read(unit_var,'(7a16)')            sta_id(i), lat_str(i), lon_str(i), elv_str(i), pm_aqi_str(i), pm_measured_str(i), pm_str(i)
 end do


! open output files and write bufr table to output files

 unit_out=12
 unit_table=13
 open(unit_table,file='pm.bufrtable')
 open(unit_out,file=outfile,action='write',form='unformatted')
 call openbf(unit_out,'OUT',unit_table)

 !call maxout(20000) ! 20k byte any message size
 call datelen(10)

 ! unchaged variables for each profile at current time
 subset="ANOWPM"
 read(analysis_time,'(I10)')  idate    ! string to integer


 do j=1, cnt

  hdr=bmiss
  obs=bmiss;

  ! write out pm which is not empty as single character ", but has valule e.g. "2.0
  pm_len=len(trim(adjustl(pm_str(j))))
  if (pm_len > 1) then

    str=adjustl(lat_str(j))
    read(str(2:),'(f15.6)') lat

    str=adjustl(lon_str(j))
    read(str(2:),'(f15.6)') lon

    str=adjustl(elv_str(j))
    read(str(2:),'(f15.6)') elv

    str=adjustl(pm_str(j))
    read(str(2:),'(f15.6)') pm

    str=adjustl(sta_id(j))
    c_sid=str(2:9)
    read(c_sid,'(z8)',iostat=ios) z_sid
    if (ios/=0) then 
      c_sid="9999999"
      read(c_sid,'(z8)',iostat=ios) z_sid
      if (ios/=0) write (*,'(a,a25)') 'wrong station_id ', str
    end if  
    
    hdr(1)=rstation_id; hdr(2)=lon; hdr(3)=lat; hdr(4)=0.0_8; hdr(5)=102; hdr(10)=6 ! Single level report


    nlvl=1
    obs(1,nlvl)=-1.0_8

    if (pm > 0.0) then
      obs(2,nlvl)=0.0_8
      obs(3,nlvl)=pm*1e-9_8  ! "UG/M3" to KG/M3
      write(*,'(a25,i7,3x,2a14,3f10.2,f15.12)') 'goodObs: pm>0',j,str,c_sid,lat,lon,obs(2,nlvl),obs(3,nlvl)
    else if (pm > pm_limit) then
      obs(2,nlvl)=0.0_8
      obs(3,nlvl)=0.0_8
      write(*,'(a25,i7,3x,2a14,3f10.2,f15.12)') 'limitObs: -15<pm<0',j,str,c_sid,lat,lon,obs(2,nlvl),obs(3,nlvl)
    else
      obs(2,nlvl)=bmiss
      obs(3,nlvl)=bmiss
      write(*,'(a25,i7,3x,2a14,3f10.2,f15.12)') 'badObs: pm<-15',j,str,c_sid,lat,lon,obs(2,nlvl),obs(3,nlvl)
    endif

    ! only encode hdr and obs, no err and qc since these are not defined in pm.bufrtable 
    call openmb(unit_out,subset,idate)
    call ufbint(unit_out,hdr,mxmn,1,iret_hdr,hdstr)
    call ufbint(unit_out,obs,mxmn,nlvl,iret_ob,obstr)
    call writsb(unit_out)

  endif
 enddo

 write(*,*) "Process PM is successful"

 call closbf(unit_out)
 close(unit_var)

end program

