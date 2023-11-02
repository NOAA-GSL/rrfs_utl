program minrhgrid

!***************************************************************************************
! abstract: use previous 12 URMA temperature and dew point analysis and/or guess files
!           files to compute mimimum RH over the previous 12 hours.
!           Note there is no dependency on time zone here.
!
! parameters: bnum - namelist variable to define whether to use guess, analysis or both (via blend_input)
!             grids - number of grids (via grid_input)
!             gridnames - names of grids to calculate for (via grid_input)
!             r100,r1000,tiny_r_kind,zero - constants from constants.f90
!             various grid names via param.incl
!
! inputs:   *_temp_ges_hour_"//hchar//".bin" - last 12 binary files of dew point and
!           temperature for URMA analysis,guess, or both.
!
! outputs:  mimRH field (binary IEEE file) for analysis and/or guess file.
!           Naming style is: hiurma.${PDYm1}.mimrh_anl.dat
!           (replace hiurma with other domain name and/or ges with anl)
!
! program history log:
!    2017-04-24   levine
!
!***************************************************************************************

  use constants, only: r100,r1000,tiny_r_kind,zero
  use kinds, only: r_kind,i_kind

  implicit none

  include 'param.incl'

  integer(i_kind),parameter::grids_total=10 !process up to this many grids
  integer(i_kind)::n,grids,bnum,g,nx,ny,fnuma,fnumb,ds
  character(len=40)::cgrid
  character(len=7)::dname
  character(len=10)::gridnames(grids_total)
  logical::fexist1,fexist2,fexist3,fexist4
  character(len=2)::hchar

  namelist/blendinput/bnum
  namelist/gridsinfo/grids,gridnames

  inquire(file="blend_input",exist=fexist1)
  if (fexist1) then
     open(11,file="blend_input",form="formatted")
     read(11,blendinput)
     close(11)
     if (bnum.ne.1.and.bnum.ne.2.and.bnum.ne.3) then
        print*, "WARNING: Invalid blend number input!  Must be 1,2, or 3.  Value given is:",bnum
        print*, "Defaulting to max analysis from analyses only (option 2)."
     endif
  else
     print*, "WARNING: No blended input value given as input.  Defaulting to max analysis from analyses only (option 2)."
     bnum=2
  endif

  inquire(file="gridsinfo_input",exist=fexist2)
  if (fexist2) then
     open(12,file="gridsinfo_input",form="formatted")
     read(12,gridsinfo)
     close(12)
     if (grids.gt.10) then
        print*, "FATAL ERROR: Manimum of 10 grids (grid names) allowed in gridsinfo_input file.  Number of grids=",grids
        stop
     endif
  else
     print*, "WARNING: no girdsinfo_input file present.  Defaulting to cohreswexp, akhres, hawaii, and prico grids."
     print*, "Check ex-script to makesure gridsinfo_input is being writted properly."
     grids=4
     gridnames(1)="cohreswexp"
     gridnames(2)="akhres"
     gridnames(3)="hawaii"
     gridnames(4)="prico"
  endif

  !print*, "in main: grids=",grids
  do n=1,grids
     print*, "In main: n,gridnames(n)=",n,gridnames(n)
  enddo
  
  do g=1,grids
     cgrid=trim(gridnames(g))

     call domain_dims(cgrid,nx,ny,ds)
     print*, "In main: nx,ny,ds=",nx,ny,ds

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
     else
        print*, "WARNING: Invalid domain name.  Grid number,cgrid=",g,cgrid
        print*, "Skipping this domain."
        cycle
     endif

  !no need for hour number to be included, as files are order by hour number (1-13), rather than hour (0-23)

     if (dname.eq."rtma3d") then
        fnuma=61
        fnumb=62
        !call pullandcalc(nx,ny,fnuma,bnum)
     elseif (dname.eq."akurma") then
        fnuma=63
        fnumb=64
        !call pullandcalc(nx,ny,fnuma,bnum)
     elseif (dname.eq."prurma") then
        fnuma=65
        fnumb=66
        !call pullandcalc(nx,ny,fnuma,bnum)
     elseif (dname.eq."hiurma") then
        fnuma=67
        fnumb=68
        !call pullandcalc(nx,ny,fnuma,bnum)
     elseif (dname.eq."guurma") then
        fnuma=69
        fnumb=70
        !call pullandcalc(nx,ny,fnuma,bnum)
        !call pullandcalc(nx,ny,fnuma,bnum)
     endif
     call pullandcalc(nx,ny,fnuma,bnum,dname)
  enddo

end program
