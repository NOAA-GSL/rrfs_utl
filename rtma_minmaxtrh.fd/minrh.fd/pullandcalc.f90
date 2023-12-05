subroutine pullandcalc(nx,ny,fnuma,bnum,dname)

  use kinds, only: r_kind,i_kind
  use constants, only: r100,r1000,tiny_r_kind,zero

  implicit none

  include 'param.incl'

  !input variables
  integer(i_kind),intent(in)::nx,ny,fnuma,bnum
  character(len=7),intent(in)::dname

  !no output variables, as output is in file

  !internal variables

  integer(i_kind)::i,j,k,singledim,ptnum,lun
  integer(i_kind),parameter::nhrs=13
  real(4),allocatable::tempanl(:,:,:),dewptanl(:,:,:),rhanl(:,:,:),tempbg(:,:,:),dewptbg(:,:,:),rhbg(:,:,:)
  real(4),allocatable::minrh(:,:),tempfld(:,:),dewfld(:,:)
  real(4)::intemp,indew,outrh,anlmin,bgmin
  logical::fexist1,fexist2,fexist3,fexist4
  character(len=2)::hchar

  !allocate and initialize arrays
  !use value tiny_r_kind from constants module for missing value

  allocate(tempanl(nx,ny,nhrs))
  allocate(dewptanl(nx,ny,nhrs))
  allocate(rhanl(nx,ny,nhrs))
  allocate(minrh(nx,ny))
  allocate(tempfld(nx,ny))
  allocate(dewfld(nx,ny))
  if (bnum.eq.1.or.bnum.eq.3) then
     allocate(tempbg(nx,ny,nhrs))
     allocate(dewptbg(nx,ny,nhrs))
     allocate(rhbg(nx,ny,nhrs))
     rhbg(:,:,:)=-r1000!tiny_r_kind
     dewptbg(:,:,:)=-r1000!tiny_r_kind
     tempbg(:,:,:)=-r1000!tiny_r_kind
  endif
  minrh(:,:)=-r1000!tiny_r_kind
  tempanl(:,:,:)=-r1000!tiny_r_kind
  dewptanl(:,:,:)=-r1000!tiny_r_kind
  rhanl(:,:,:)=-r1000!tiny_r_kind
  !tempfld(:,:)=-r1000!tiny_r_kind
  dewfld(:,:)=-r1000!tiny_r_kind
  singledim=nx*ny

  print*, "In pullandcalc:  Initializing nx,ny,singledim=",nx,ny,singledim

  lun=11
  do k=1,nhrs
     write(hchar,"(I2.2)") k
     print*, "Looking at hour:",hchar
     if (bnum.eq.1.or.bnum.eq.3) then
        inquire(file=trim(dname)//"_temp_ges_hour_"//hchar//".bin",exist=fexist1)
        inquire(file=trim(dname)//"_dwpt_ges_hour_"//hchar//".bin",exist=fexist2)
        if (fexist1.and.fexist2) then
           open(lun,file=trim(dname)//"_temp_ges_hour_"//hchar//".bin",form="unformatted",status="old",access="sequential")
           open(lun+1,file=trim(dname)//"_dwpt_ges_hour_"//hchar//".bin",form="unformatted",status="old",access="sequential")
           read(lun) tempfld
           read(lun+1) dewfld
           close(lun)
           close(lun+1)
           do i=1,nx
              do j=1,ny
                 ptnum=(ny*(i-1))+j
                 if (ptnum.gt.singledim) then
                    print*, "FATAL ERROR: ptnum calculation error in guess read in (line 58) at i,j,ptnum,dname,hour=",i,j,ptnum,dname,k
                    stop
                 endif
                 if (tempfld(i,j).gt.zero.and.tempfld(i,j).lt.r1000) then
                    tempbg(i,j,k)=tempfld(i,j)
                 else
                    print*, "WARNING: Invalid temperature guess at dname,hour,i,j,value=",dname,i,j,tempfld(i,j)
                 endif
                 if (dewfld(i,j).gt.zero.and.dewfld(i,k).lt.r1000) then
                    dewptbg(i,j,k)=dewfld(i,j)
                 else
                    print*, "WARNING: Invalid dew point guess at dname,hour,i,j,value=",dname,i,j,dewfld(i,j)
                 endif
              enddo !j=1,ny
           enddo !i=1,nx
        elseif (fexist1) then
           print*, "FATAL ERROR: Dew Point binary guess file missing for domain,hour:",dname,hchar
           stop
        elseif (fexist2) then
           print*, "FATAL ERROR: Temperature binary guess file missing for domain,hour:",dname,hchar
           stop
        else !
           print*, "FATAL ERROR: Temperature and dew point binary guess files missing for domain,hour:",dname,hchar
           stop
        endif
     endif
     if (bnum.eq.2.or.bnum.eq.3) then
        inquire(file=trim(dname)//"_temp_anl_hour_"//hchar//".bin",exist=fexist3)
        inquire(file=trim(dname)//"_dwpt_anl_hour_"//hchar//".bin",exist=fexist4)
        if (fexist3.and.fexist4) then
           open(lun+2,file=trim(dname)//"_temp_anl_hour_"//hchar//".bin",form="unformatted")!,status="old",access="sequential")
           open(lun+3,file=trim(dname)//"_dwpt_anl_hour_"//hchar//".bin",form="unformatted")!,status="old",access="sequential")
           read(lun+2) tempfld
           read(lun+3) dewfld
           close(lun+2)
           close(lun+3)
           do i=1,nx
              do j=1,ny
                 ptnum=(ny*(i-1))+j
                 if (ptnum.gt.singledim) then
                    print*, "FATAL ERROR: ptnum calculation error in analysis read in (line 58) at i,j,ptnum,dname,hour=",i,j,ptnum,dname,k
                    stop
                 endif
                 if (tempfld(i,j).gt.zero.and.tempfld(i,j).lt.r1000) then
                    tempanl(i,j,k)=tempfld(i,j)
                 else
                    print*, "WARNING: Invalid temperature analysis at dname,hour,i,j,value=",dname,i,j,tempfld(i,j)
                 endif
                 if (dewfld(i,j).gt.zero.and.dewfld(i,k).lt.r1000) then
                    dewptanl(i,j,k)=dewfld(i,j)
                 else
                    print*, "WARNING: Invalid dew point analysis at dname,hour,i,j,value=",dname,i,j,dewfld(i,j)
                 endif
              enddo !j=1,ny
           enddo !i=1,nx
        elseif (fexist3) then
           print*, "FATAL ERROR: Dew Point binary analysis file missing for domain,hour:",dname,hchar
           stop
        elseif (fexist4) then
           print*, "FATAL ERROR: Temperature binary analysis file missing for domain,hour:",dname,hchar
           stop
        else !
           print*, "FATAL ERROR: Temperature and dew point binary analysis files missing for domain,hour:",dname,hchar
           stop
        endif
     endif
     lun=lun+4
  enddo
  deallocate(tempfld)
  !now load RH's
  do k=1,nhrs
     do i=1,nx
        do j=1,ny
           if (bnum.eq.2.or.bnum.eq.3) then
              intemp=tempanl(i,j,k)
              indew=dewptanl(i,j,k)
              call calcrh(intemp,indew,outrh)
              rhanl(i,j,k)=outrh
           endif
           if (bnum.eq.1.or.bnum.eq.3) then
              intemp=tempbg(i,j,k)
              indew=dewptbg(i,j,k)
              call calcrh(intemp,indew,outrh)
              rhbg(i,j,k)=outrh
           endif
        enddo
     enddo
  enddo
  deallocate(dewfld)

  if (bnum.eq.1) then
     do i=1,nx
        do j=1,ny
           minrh(i,j)=minval(rhbg(i,j,:))
           if (j.eq.100.and.i.eq.100) then
              print*, "Min RH value at 100,100=",minrh(100,100)
           endif
        enddo
     enddo
     write(fnuma) minrh
     deallocate(tempbg)
     deallocate(dewptbg)
     deallocate(rhbg)
  elseif (bnum.eq.2) then
     do i=1,nx
        do j=1,ny
           minrh(i,j)=minval(rhanl(i,j,:))
           if (j.eq.100.and.i.eq.100) then
              print*, "Min RH value at 100,100=",minrh(100,100)
           endif
        enddo
     enddo
     write(fnuma) minrh
     deallocate(tempanl)
     deallocate(dewptanl)
     deallocate(rhanl)
  !endif
  elseif (bnum.eq.3) then
     do i=1,nx
        do j=1,ny
           bgmin=minval(rhbg(i,j,:))
           anlmin=minval(rhanl(i,j,:))
           if (bgmin.le.anlmin) then
              minrh(i,j)=bgmin
           else
              minrh(i,j)=anlmin
           endif
        enddo
     enddo
     write(fnuma) minrh
     deallocate(tempanl)
     deallocate(dewptanl)
     deallocate(rhanl)
     deallocate(tempbg)
     deallocate(dewptbg)
     deallocate(rhbg)
  endif
  close(fnuma)
  deallocate(minrh)
  !close(fnumb)

end subroutine pullandcalc
