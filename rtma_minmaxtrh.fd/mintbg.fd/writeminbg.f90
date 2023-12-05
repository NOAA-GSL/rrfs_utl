subroutine readminb(nx,ny,hours,threehours,dname,fnum,bnum)

  !abstract: Read in background values from binary files for an URMA
  !calculate min from them and write to an output file
  
  !need explicit interface in main program
  
  implicit none
  
  include 'param.incl'
  
  !input variables
  !real(4),intent(in)::misreal,missing !values will not be changed
  !real(4),intent("in"),parameter::missing=9.999e20
  integer(4),intent(in)::nx,ny,fnum!,nhours,max !,singledim calc singledim on the fly
  !integer(4),intent(in)::tzsec(nx,ny),tzhour(nx,ny)
  integer(4),intent(in)::hours(nhours),threehours(8)
  character(len=7),intent(in)::dname
  
  !internal variables
  integer(4)::i,j,k,l,lun,ptnum,singledim,bnum !counters
  real(4),allocatable::tzsec(:,:),tzhour(:,:)
  real(4),allocatable::temps(:,:)
  real(4),allocatable::hourlytempbg(:,:,:),hourlytempanl(:,:,:)
  real(4),allocatable::threehour(:,:,:)
  real(4),allocatable::sixhour(:,:,:)
  real(4),allocatable::allhrlybg(:,:),allhrlyanl(:,:),allsix(:,:)
  real(4),allocatable::finalmin(:,:)
  real(4),allocatable::mintonebg(:),minttwobg(:,:),minttwoanl(:,:),mintoneanl(:)
  integer(4),allocatable::alltz(:)
  !real(4),allocatable::maxtone(:),maxttwo(:,:)
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
  !allocate(temps(nx,ny))
  singledim=nx*ny
  print*, "maximum singledim=",singledim
  !singledim=1000
  allocate(alltz(singledim))
  allocate(allhrlybg(singledim,nhours))
  allocate(allhrlyanl(singledim,nhours))
  allocate(allsix(singledim,4))
  allocate(mintonebg(singledim))
  allocate(mintoneanl(singledim))
  allocate(minttwobg(nx,ny))
  allocate(minttwoanl(nx,ny))
  allocate(finalmin(nx,ny))
  !allocate(maxtone(singledim))
  !allocate(maxttwo(nx,ny))
  
  !initialize arrays
  
  !temps(:,:)=misreal
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
  mintonebg(:)=misreal
  mintoneanl(:)=misreal
  minttwobg(:,:)=misreal
  minttwoanl(:,:)=misreal
  finalmin(nx,ny)=misreal
  
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
     inquire(file=trim(dname)//"_ges_valid18_prevday.bin",exist=fexist)
     if (fexist) then
        open(12,file=trim(dname)//"_ges_valid18_prevday.bin",form="unformatted")
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
              !hourlytemp(i,j,1)=temps(i,j)
              ptnum=(ny*(i-1))+j
              if (ptnum.gt.singledim) then
                 print*, "FATAL ERROR: ptnum error (line 154 - 1st hour) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,1
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
  if (bnum.eq.2.or.bnum.eq.3) then
     inquire(file=trim(dname)//"_anl_valid18_prevday.bin",exist=fexist)
     if (fexist) then
        open(12,file=trim(dname)//"_anl_valid18_prevday.bin",form="unformatted")
        read(12) temps
        close(12)
        !convert 2d into 1d array, match with hour number (1) to but back into 2d array needed for MDL routine
        do i=1,nx
           do j=1,ny
              if ((temps(i,j).gt.10000).or.(temps(i,j).lt.-10000)) then
                 hourlytempanl(i,j,1)=misreal
              else
                 hourlytempanl(i,j,1) = temps(i,j)
              endif
              !hourlytemp(i,j,1)=temps(i,j)
              ptnum=(ny*(i-1))+j
              if (ptnum.gt.singledim) then
                 print*, "FATAL ERROR: ptnum error (line 182 - 1st hour) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,1
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
     lun=k+12
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
                    print*, "FATAL ERROR: ptnum error (line 209 - hourly) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,k
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
           open(lun,file=trim(dname)//"_anl_valid"//hchar//".bin",form="unformatted")
           read(lun) temps
           close(lun)
           do i=1,nx
              do j=1,ny
                 ptnum=(ny*(i-1))+j
                 if (ptnum.gt.singledim) then
                    print*, "FATAL ERROR: ptnum error (line 236 - hourly) at i,j,ptnum,singledim,dname,hour=",i,j,ptnum,singledim,dname,k
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
           allhrlyanl(:,k)=misreal
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
           print*, "FATAL ERROR: ptnum error (line 265 - timezone) at i,j,ptnum,singledim,dname=",i,j,ptnum,singledim,dname
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
     call maxmin(fnum,allhrlybg,allsix,alltz,names,desc,mintonebg,singledim,4,25,1,2,singledim)
  endif
  if (bnum.eq.2.or.bnum.eq.3) then
     call maxmin(fnum+1,allhrlyanl,allsix,alltz,names,desc,mintoneanl,singledim,4,25,1,2,singledim)
  endif
  do l=1,singledim
     !do l=1,1000
     j=mod(l,ny)
     if (j.eq.0) then
        j=ny
     endif
     i=((l-j)/ny)+1
     if (bnum.eq.1.or.bnum.eq.3) then
        if (mintonebg(l).eq.misreal) then
           minttwobg(i,j)=missing
        elseif (nint(mintonebg(l)).eq.0) then
           minttwobg(i,j)=missing
        else
           minttwobg(i,j)=mintonebg(l)
        endif
     endif
     if (bnum.eq.2.or.bnum.eq.3) then
        if (mintoneanl(l).eq.misreal) then
           minttwoanl(i,j)=missing
        elseif (nint(mintoneanl(l)).eq.0) then
           minttwoanl(i,j)=missing
        else
           minttwoanl(i,j)=mintoneanl(l)
        endif
     endif
     !assign max value based on max between anl and ges at i,j
     if (bnum.eq.3) then
     !   if (minttwoanl(i,j).le.minttwobg(i,j)) then
     !      finalmin(i,j)=minttwoanl(i,j)
     !   else
     !      finalmin(i,j)=minttwobg(i,j)
     !   endif
     !elseif (bnum.eq.2) then
     !   finalmin(i,j)=minttwoanl(i,j)
     !elseif (bnum.eq.1) then
     !   finalmin(i,j)=minttwobg(i,j)
     !endif
        if ((minttwoanl(i,j).eq.missing).and.(minttwobg(i,j).eq.missing)) then
           finalmin(i,j)=missing
        elseif ((minttwoanl(i,j).eq.missing).and.(minttwobg(i,j).ne.missing)) then
           finalmin(i,j)=minttwobg(i,j)
        elseif ((minttwoanl(i,j).ne.missing).and.(minttwobg(i,j).eq.missing)) then
           finalmin(i,j)=minttwoanl(i,j)
        elseif (minttwoanl(i,j).le.minttwobg(i,j)) then
           finalmin(i,j)=minttwoanl(i,j)
        else
           finalmin(i,j)=minttwobg(i,j)
        endif
     elseif (bnum.eq.2) then
        finalmin(i,j)=minttwoanl(i,j)
     elseif (bnum.eq.1) then
        finalmin(i,j)=minttwobg(i,j)
     endif
  enddo

  if ((finalmin(1,1).eq.missing).or.(finalmin(1,1).le.0).or.(finalmin(1,1).ge.1000)) then
     print*, "WARNING: FINALMIN CONTAINS MISSING/INVALID VALUES!  GSI MAY FAIL!"
  endif
  
  open(61,file="mint_"//trim(dname)//"_bg.bin",form="unformatted")
  write(61) finalmin
  close(61)
  
  !deallocate internal arrays
  deallocate(hourlytempbg)
  deallocate(hourlytempanl)
  deallocate(threehour)
  deallocate(sixhour)
  deallocate(tzsec)
  deallocate(tzhour)
  deallocate(temps)
  !singledim=nx*ny
  !print*, "maximum singledim=",singledim
  !singledim=1000
  deallocate(alltz)
  !deallocate(allhrly)
  deallocate(allsix)
  deallocate(mintoneanl)
  deallocate(mintonebg)
  deallocate(minttwoanl)
  deallocate(minttwobg)
  deallocate(finalmin)
  deallocate(allhrlybg)
  deallocate(allhrlyanl)
  
end subroutine readminb
