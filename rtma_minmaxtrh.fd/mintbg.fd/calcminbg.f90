program bgmin

!***************************************************************************************  
! abstract: use previous URMA temperature first guess analysis and guess
!           files to compute minimum temperatuer according to NDFD definition
!           (7P-8A).  Account for local time zone using gridded time zone file (fixed file).
!
! paramaters: names of grids via grids_input
!
! inputs: *_ges_valid${hchar}.bin - last 24 guess files for temperature.  hchar = 2 character hour
!         *_anl_valid${hchar}.bin - last 24 analysis files for temperature.  hchar = 2 character hour 
!
! outputs: rtma3d.${PDY}.maxT.bin, rtma3d.${PDY}.minT.grb2
!          akurma.${PDY}.minT.bin, akurma.${PDY}.minT.grb2
!          hiurma.${PDY}.minT.bin, hiurma.${PDY}.minT.grb2
!          prurma.${PDY}.minT.bin, prurma.${PDY}.minT.grb2
!
! program history log:
!    2017-04-24   levine
!
!***************************************************************************************

  implicit none
  
  interface
     subroutine readminb(nx,ny,hours,threehours,dname,fnum,bnum)
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
       
     end subroutine readminb
  end interface
  
  include 'param.incl'
  
  integer(4)::g,i,j,k,l,nx,ny,lun,fnum,bnum!,singledim,ptnum
  real(4)::ds
  integer(4)::hours(nhours),threehours(8)
  integer(4),parameter::grids_total=10  !process up to this number of grids
  character(len=40)::cgrid,mdlpath
  character(len=2)::hchar
  character(len=7)::dname
  character(len=10)::gridnames(grids_total)
  logical::fexist

  integer(4) n,grids

  namelist/blendinput/bnum
  namelist/gridsinfo/grids,gridnames

  open(15,file="blend_input",form="formatted")
  read(15,blendinput)
  close(15)

  if (bnum.ne.1.and.bnum.ne.2.and.bnum.ne.3) then
     print*, "FATAL ERROR: Invalid bnum option.  Must be 1,2, or 3.  Current value is:",bnum
     stop
  endif
  
  inquire(file='gridsinfo_input',exist=fexist)
  if(fexist) then
     open (15,file='gridsinfo_input',form='formatted')
     read (15,gridsinfo)
     close(15)
     if (grids.gt.10) then
        print*, "FATAL ERROR: Maximum of 10 grids (grid names) allowed in dominfo_input file.  Current value is:",grids
        stop
     endif
  else
     print*, "WARNING: no dominfo_input file present.  Defaulting to cohreswexp, akhres, hawaii, and prico grids."
     print*, "Check ex-script to make sure dominfo_input file is being written properly."
     grids=4
     gridnames(1)="cohreswexp"
     gridnames(2)="akhres"
     gridnames(3)="hawaii"
     gridnames(4)="prico"
  endif
  
  print*,'in calcmaxbg: grids=',grids

  do n=1, grids
     print*,'in calcmaxbg: n,gridnames(n)=',n,gridnames(n)
  enddo
  
  do g=1,grids
     cgrid=trim(gridnames(g))
     
     call domain_dims(cgrid,nx,ny,ds)
     print*, "In create_firstguess: nx,ny,ds=",nx,ny,ds
     
     if (cgrid.eq."cohresext") then
        dname="rtma3d"
     elseif (cgrid.eq."cohreswexp") then
        dname="rtma3d"
     elseif (cgrid.eq."cohres") then
        dname="rtma3d"
     elseif (cgrid.eq."hrrr") then
        dname="rtma3d"
     elseif (cgrid.eq."akhres") then
        dname="akurma"
     elseif (cgrid.eq."prico") then
        dname="prurma"
     elseif (cgrid.eq."hawaii") then
        dname="hiurma"
     elseif (cgrid.eq."guam") then
        dname="guurma"
     endif
     
     if (cgrid.ne."guam") then
        
        hours=(/18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
        threehours=(/21,0,3,6,9,12,15,18/)
        
        if (dname.eq."rtma3d") then
           fnum=71
           call readminb(nx,ny,hours,threehours,dname,fnum,bnum)
        elseif (dname.eq."akurma") then
           fnum=73
           call readminb(nx,ny,hours,threehours,dname,fnum,bnum)
        elseif (dname.eq."prurma") then
           fnum=75
           call readminb(nx,ny,hours,threehours,dname,fnum,bnum)
        elseif (dname.eq."hiurma") then
           fnum=77
           call readminb(nx,ny,hours,threehours,dname,fnum,bnum)
        else
           print*, "FATAL ERROR: Invalid domain name:",dname
           stop
        endif
     else
        print*, "Skipping Guam for now..."
     endif
  enddo
  
end program bgmin
