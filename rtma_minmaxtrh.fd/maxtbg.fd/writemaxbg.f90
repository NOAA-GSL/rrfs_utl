subroutine readmaxb(nx,ny,hours,threehours,dname,fnum,bnum)

  !abstract: Read in background values from binary files for an URMA
  !calculate max from them and write to an output file

  !need explicit interface in main program

  implicit none

  include 'param.incl'

  !input variables
  integer(4),intent(in)::nx,ny,fnum,bnum!,nhours,max !,singledim calc singledim on the fly
  integer(4),intent(in)::hours(nhours),threehours(8)
  character(len=7),intent(in)::dname

  !internal variables
  integer(4)::i,j,k,l,lun,ptnum,singledim !counters
  real(4),allocatable::tzsec(:,:),tzhour(:,:)
  real(4),allocatable::temps(:,:)
  real(4),allocatable::hourlytempbg(:,:,:),hourlytempanl(:,:,:)
  real(4),allocatable::threehour(:,:,:)
  real(4),allocatable::sixhour(:,:,:)
  real(4),allocatable::allhrlybg(:,:),allhrlyanl(:,:),allsix(:,:)
  real(4),allocatable::finalmax(:,:)
  integer(4),allocatable::alltz(:)
  real(4),allocatable::maxtonebg(:),maxttwobg(:,:),maxttwoanl(:,:),maxtoneanl(:)
  character(len=20)::desc(max)
  character(len=8)::names(max)
  logical::fexist
  character(len=2)::hchar

  !allocate internal arrays, calculate singledim

  allocate(hourlytempbg(nx,ny,nhours))
  allocate(hourlytempanl(nx,ny,nhours))
  allocate(threehour(nx,ny,8))
  allocate(sixhour(nx,ny,4))
  allocate(tzsec(nx,ny))
  allocate(tzhour(nx,ny))
  singledim=nx*ny
  print*, "maximum singledim=",singledim
  allocate(alltz(singledim))
  allocate(allhrlybg(singledim,nhours))
  allocate(allhrlyanl(singledim,nhours))
  allocate(allsix(singledim,4))
  allocate(maxtonebg(singledim))
  allocate(maxtoneanl(singledim))
  allocate(maxttwobg(nx,ny))
  allocate(maxttwoanl(nx,ny))
  allocate(finalmax(nx,ny))

  !initialize arrays

  tzhour(:,:)=misreal
  tzsec(:,:)=misreal
  sixhour(:,:,:)=misreal
  threehour(:,:,:)=misreal
  hourlytempbg(:,:,:)=misreal
  hourlytempanl(:,:,:)=misreal
  allsix(:,:)=misreal
  allhrlybg(:,:)=misreal
  allhrlyanl(:,:)=misreal
  alltz(:)=misreal
  maxtonebg(:)=misreal
  maxtoneanl(:)=misreal
  maxttwobg(:,:)=misreal
  maxttwoanl(:,:)=misreal
  finalmax(nx,ny)=misreal

  !configure time zone

  if (dname.eq."rtma3d".and.nx.eq.2145.and.ny.eq.1597) then
     open(11,file="conustz.bin",form="unformatted") !time zone file
     read(11) tzsec
     close(11)
  elseif (dname.eq."rtma3d".and.nx.eq.2345.and.ny.eq.1597) then
     open(11,file="conusexttz.bin",form="unformatted") !time zone file
     read(11) tzsec
     close(11)
  elseif (dname.eq."rtma3d".and.nx.eq.2145.and.ny.eq.1377) then
     open(11,file="conustz_ndfdonly.bin",form="unformatted") !time zone file
     read(11) tzsec
     close(11)
  elseif (dname.eq."rtma3d".and.nx.eq.1799.and.ny.eq.1059) then
     open(11,file="conustz_ndfdonly.bin",form="unformatted") !time zone file
     read(11) tzsec
     close(11)
  elseif (dname.eq."akurma") then
     open(11,file="aktz.bin",form="unformatted")
     read(11) tzsec
     close(11)
  elseif (dname.eq."prurma") then
     tzsec(:,:)=-14400
     tzhour(:,:)=-4
  elseif (dname.eq."hiurma") then
     tzsec(:,:)=-32400
     tzhour(:,:)=-9
  elseif (dname.eq."guurma") then
     tzsec(:,:)=14400
     tzhour(:,:)=5
  else
     print*, "FATAL ERROR: Invalid domain name,nx,ny:",dname,nx,ny
     stop
  endif
  tzhour=tzsec/3600
  if (tzsec(1,1).lt.-50000) then
     print*, "FATAL ERROR: Error reading time zone file, value at 1,1=",tzsec(1,1)
     stop
  endif

  allocate(temps(nx,ny))
  temps(:,:)=misreal

  !apply fix for AK region

  if (dname.eq."akurma") then
     do i=1,nx
        do j=1,ny
           !fix time zone to -9 across domain at AR request
           if ((tzhour(i,j).lt.-9).or.(tzhour(i,j).gt.0)) then
              tzhour(i,j)=-9
           endif
        enddo
     enddo
     desc(:)="DESC GOES HERE    AK"
  else
     desc(:)="STN DESC GOES HERE  "
  endif
  names(:)="12345678"
  
  !read in first hourly temps
  if (bnum.eq.1.or.bnum.eq.3) then
     inquire(file=trim(dname)//"_ges_valid06_prevday.bin",exist=fexist)
     if (fexist) then
        open(12,file=trim(dname)//"_ges_valid06_prevday.bin",form="unformatted")
        read(12) temps
        close(12)
        !convert 2d into 1d array, match with hour number (1) to but back into 2d array needed for MDL routine
        do i=1,nx
           do j=1,ny
              if ((temps(i,j).gt.10000).or.(temps(i,j).lt.-10000)) then
                 hourlytempbg(i,j,1)=misreal
              else
                 hourlytempbg(i,j,1) = temps(i,j)
              endif
              ptnum=(ny*(i-1))+j
              if (ptnum.gt.singledim) then
                 print*, "FATAL ERROR: ptnum calc error (line 145 - 1st hour) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,1
                 stop
              endif
              allhrlybg(ptnum,1)=hourlytempbg(i,j,1)
           enddo
        enddo
     else
        print*, "WARNING: HOURLY GUESS FILE MISSING:",dname,1
        allhrlybg(:,1)=misreal
     endif
  endif

  !now the analysis file
  if (bnum.eq.2.or.bnum.eq.3) then
     inquire(file=trim(dname)//"_anl_valid06_prevday.bin",exist=fexist)
     if (fexist) then
        open(13,file=trim(dname)//"_anl_valid06_prevday.bin",form="unformatted")
        read(13) temps
        close(13)
        do i=1,nx
           do j=1,ny
              if ((temps(i,j).gt.10000).or.(temps(i,j).lt.-10000)) then
                 hourlytempanl(i,j,1)=misreal
              else
                 hourlytempanl(i,j,1) = temps(i,j)
              endif
              ptnum=(ny*(i-1))+j
              if (ptnum.gt.singledim) then
                 print*, "FATAL ERROR: ptnum calcerror (line 173 - 1st hour) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,1
                 stop
              endif
              allhrlyanl(ptnum,1)=hourlytempanl(i,j,1)
           enddo
        enddo
     else
        print*, "WARNING: HOURLY ANALYSIS FILE MISSING:",dname,1
        allhrlyanl(:,1)=misreal
     endif
  endif

  !now do the same thing for the other 24 hours
  do k=2,nhours
     lun=(k*2)+12
     write(hchar,"(I2.2)") hours(k)
     if (bnum.eq.1.or.bnum.eq.3) then
        inquire(file=trim(dname)//"_ges_valid"//hchar//".bin",exist=fexist)
        if (fexist) then
           open(lun,file=trim(dname)//"_ges_valid"//hchar//".bin",form="unformatted")
           read(lun) temps
           close(lun)
           do i=1,nx
              do j=1,ny
                 ptnum=(ny*(i-1))+j
                 if (ptnum.gt.singledim) then
                    print*, "FATAL ERROR: ptnum calc error (line 199 - hourly) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,k
                    stop
                 endif
                 if ((temps(i,j).gt.10000).or.(temps(i,j).lt.-10000)) then
                    hourlytempbg(i,j,k)=misreal
                    allhrlybg(ptnum,k)=misreal
                 else
                    hourlytempbg(i,j,k)=temps(i,j)
                    allhrlybg(ptnum,k)=temps(i,j)
                 endif
              enddo
           enddo
        else
           print*, "WARNING: HOURLY GUESS FILE MISSING:",dname,k
           allhrlybg(:,k)=misreal
        endif
     endif
     if (bnum.eq.2.or.bnum.eq.3) then
        inquire(file=trim(dname)//"_anl_valid"//hchar//".bin",exist=fexist)
        if (fexist) then
           lun=lun+1
           open(lun,file=trim(dname)//"_anl_valid"//hchar//".bin",form="unformatted")
           read(lun) temps
           close(lun)
           do i=1,nx
              do j=1,ny
                 ptnum=(ny*(i-1))+j
                 if (ptnum.gt.singledim) then
                    print*, "FATAL ERROR: ptnum calc error (line 227 - hourly) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,k
                    stop
                 endif
                 if ((temps(i,j).gt.10000).or.(temps(i,j).lt.-10000)) then
                    hourlytempanl(i,j,k)=misreal
                    allhrlyanl(ptnum,k)=misreal
                 else
                    hourlytempanl(i,j,k)=temps(i,j)
                    allhrlyanl(ptnum,k)=temps(i,j)
                 endif
              enddo
           enddo
        else
           print*, "WARNING: HOURLY ANALYSIS FILE MISSING:",dname,k
           allhrlyanl(ptnum,k)=misreal
        endif
     endif
  enddo

  !we are not using 3 or 6 hourly arrays, so default to misreal
  threehour(:,:,:)=misreal
  sixhour(:,:,:)=misreal
  allsix(:,:)=misreal
  
  !convert time zones for all points into 1D array (needed for MDL routine)
  do i=1,nx
     do j=1,ny
        ptnum=(ny*(i-1))+j
        if (ptnum.gt.singledim) then
           print*, "FATAL ERROR: ptnum calc error (line 256 - timezone) at i,j,ptnum,singledim,dname=",i,j,ptnum,singledim,dname
           stop
        endif
        !print*, "Looking at time zone for ptnum,i,j=",ptnum,i,j,tzhour(i,j)
        if ((tzhour(i,j).gt.100).or.tzhour(i,j).lt.-100) then
           alltz(ptnum)=alltz(ptnum-1)
        else
           alltz(ptnum)=nint(tzhour(i,j))
        endif
     enddo
  enddo
  
  print*, "Running max calculation..."
  if (bnum.eq.1.or.bnum.eq.3) then
     call maxmin(fnum,allhrlybg,allsix,alltz,names,desc,maxtonebg,singledim,4,25,2,2,singledim)
  endif
  if (bnum.eq.2.or.bnum.eq.3) then
     call maxmin(fnum+1,allhrlyanl,allsix,alltz,names,desc,maxtoneanl,singledim,4,25,2,2,singledim)
  endif
  do l=1,singledim
     j=mod(l,ny)
     if (j.eq.0) then
        j=ny
     endif
     i=((l-j)/ny)+1
     if (bnum.eq.1.or.bnum.eq.3) then
        if (maxtonebg(l).eq.misreal) then
           maxttwobg(i,j)=missing
        elseif (nint(maxtonebg(l)).eq.0) then
           maxttwobg(i,j)=missing
        else
           maxttwobg(i,j)=maxtonebg(l)
        endif
     endif
     if (bnum.eq.2.or.bnum.eq.3) then
        if (maxtoneanl(l).eq.misreal) then
           maxttwoanl(i,j)=missing
        elseif (nint(maxtoneanl(l)).eq.0) then
           maxttwoanl(i,j)=missing
        else
           maxttwoanl(i,j)=maxtoneanl(l)
        endif
     endif
     !assign max value based on max between anl and ges at i,j
     if (bnum.eq.3) then
        !        if (maxttwoanl(i,j).ge.maxttwobg(i,j)) then
        !           finalmax(i,j)=maxttwoanl(i,j)
        !        else
        !           finalmax(i,j)=maxttwobg(i,j)
        !        endif
        !     elseif (bnum.eq.2) then
        !        finalmax(i,j)=maxttwoanl(i,j)
        !     elseif (bnum.eq.1) then
        !        finalmax(i,j)=maxttwobg(i,j)
        if ((maxttwoanl(i,j).eq.missing).and.(maxttwobg(i,j).eq.missing)) then
           finalmax(i,j)=missing
        elseif ((maxttwoanl(i,j).eq.missing).and.(maxttwobg(i,j).ne.missing)) then
           finalmax(i,j)=maxttwobg(i,j)
        elseif ((maxttwoanl(i,j).ne.missing).and.(maxttwobg(i,j).eq.missing)) then
           finalmax(i,j)=maxttwoanl(i,j)
        elseif (maxttwoanl(i,j).ge.maxttwobg(i,j)) then
           finalmax(i,j)=maxttwoanl(i,j)
        else
           finalmax(i,j)=maxttwobg(i,j)
        endif
     elseif (bnum.eq.2) then
        finalmax(i,j)=maxttwoanl(i,j)
     elseif (bnum.eq.1) then
        finalmax(i,j)=maxttwobg(i,j)
     endif
  enddo

  if ((finalmax(1,1).eq.missing).or.(finalmax(1,1).le.0).or.(finalmax(1,1).ge.1000)) then
     print*, "WARNING: MISSING/INVALID VALUES IN FINAL MAX!  GSI MAY FAIL!"
  endif
  
  open(61,file="maxt_"//trim(dname)//"_bg.bin",form="unformatted")
  write(61) finalmax
  close(61)
  
  !deallocate internal arrays
  deallocate(hourlytempbg)
  deallocate(hourlytempanl)
  deallocate(threehour)
  deallocate(sixhour)
  deallocate(tzsec)
  deallocate(tzhour)
  deallocate(temps)
  deallocate(alltz)
  deallocate(allsix)
  deallocate(maxtoneanl)
  deallocate(maxtonebg)
  deallocate(maxttwoanl)
  deallocate(maxttwobg)
  deallocate(finalmax)
  deallocate(allhrlybg)
  deallocate(allhrlyanl)
  
end subroutine readmaxb
