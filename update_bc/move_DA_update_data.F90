!-----------------------------------------------------------------------
      program move_DA_update_data
!-----------------------------------------------------------------------
!***  Full fields including the integration domain's boundary rows
!***  have been updated by the DA.  Now the interior points of the
!***  core and tracer restart files must be put back into the
!***  standard sized files and the boundary rows must be taken 
!***  from the full fields and put back into the BC files.
!-----------------------------------------------------------------------
      use netcdf
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
      integer,parameter :: double=selected_real_kind(p=13,r=200)           !<-- Define KIND for double precision real
!
      integer,parameter :: halo_integrate=3                                !<-- Halo depth for FV3 integration
!
      integer,parameter :: ndims_bc=8                                   &  !<-- # of dimensions in the BC files
                          ,ndims_res=6                                  &  !<-- # of dimensions in the core restart file
                          ,num_fields_update_core=6                     &  !<-- # of updated fields in the core restart file
                          ,num_fields_core_bc=6                            !<-- # of core fields copied into the BC files
!
      real,parameter :: ptop=200.                                          !<-- The domain's top pressure (Pa)
      real,parameter :: grav=9.81                                          !<-- g from fv_diagnostics
      real,parameter :: rd=287.04                                       &  !<-- rd from fv_diagnostics
                       ,rv=461.51                                          !<-- rv from fv_diagnostics
      real,parameter :: f608=rv/rd-1.
      real,parameter :: rgrav=1./grav
!
      real(kind=double),parameter :: pi_r_180=acos(-1._double)/180._double
!
      integer :: dimid_halo,halo,i,istat,j,k,kend &
                ,lat,latp,len_x,len_y,length,loc1,lon,lonp &
                ,n,na,natts,ncid_bc,ncid_bc_new,ncid_core_combined      &
                ,ncid_core_res,ncid_tracer_combined,ncid_tracer_res     &
                ,nctype,ndims,ngatts,nrows_blend                        &
                ,nrows_bndry,num_vars_bc,num_tracers_bc                 &
                ,nvars_res,unlimdimid                                   &
                ,var_id_bc,var_id_delp,var_id_delz,var_id_phis,var_id_ps            &
                ,var_id_res,var_id_sphum,var_id_t,var_id_tracer         &
                ,width_halo_total
!
      integer :: istart_bc,iend_bc,jstart_bc,jend_bc
      integer :: istart_res,iend_res,jstart_res,jend_res
      integer :: ie_combined,je_combined
!
      integer,dimension(ndims_bc) :: dimid_bc                               !<-- Dimension IDs of the BC file.
      integer,dimension(ndims_res) :: dimid_core                            !<-- Dimension IDs of the core restart dimensions.
!
      integer,dimension(ndims_bc) :: dimsize_bc                             !<-- Dimensions in the BC file.
      integer,dimension(ndims_res) :: dimsize_res                       &   !<-- Dimensions of the core restart domain.
                                     ,dimsize_combined                      !<-- Restart dimensions with bndry rows added.
!
      real,dimension(:),allocatable :: row4_east,row4_north             &
                                      ,row4_south,row4_west
!
      real,dimension(:,:),allocatable :: field_combined                 &
                                        ,phis_combined                  &
                                        ,sphum_combined
!
      real,dimension(:,:),allocatable,target :: ps_east,ps_north        &
                                               ,ps_south,ps_west
!
      character(len=3) :: hour
      character(len=50) :: name_att,name_var
      character(len=128) :: filename_bc_new='gfs_bndy.tile7.hhh_gsi.nc'           &   !<-- The new BC file on forecast model layers.
                           ,filename_bc='gfs_bndy.tile7.hhh.nc'                   &   !<-- The original BC file on input model layers.
                           ,filename_core_combined='fv_core.res.tile1_new.nc'     &
                           ,filename_core_restart='fv_core.res.tile1.nc'          &
                           ,filename_tracer_combined='fv_tracer.res.tile1_new.nc' &
                           ,filename_tracer_restart='fv_tracer.res.tile1.nc' 
!
      character(len=25) :: varname_update_bc
!
      character(len=7),dimension(ndims_res) :: dimname_core=(/           &
                                                              'xaxis_1'  &
                                                             ,'xaxis_2'  &
                                                             ,'yaxis_1'  &
                                                             ,'yaxis_2'  &
                                                             ,'zaxis_1'  &
                                                             ,'Time'     &
                                                                 /)
!
      character(len=5),dimension(ndims_bc) :: dimname_bc=(/             &
                                                           'lon'        &
                                                          ,'lat'        &
                                                          ,'lonp'       &
                                                          ,'latm'       &
                                                          ,'halo'       &
                                                          ,'halop'      &
                                                          ,'lev'        &
                                                          ,'levp'       &
                                                          /)
!
      character(len=4),dimension(num_fields_update_core) :: varname_update_core=(/       &
                                                                                  'u'    &
                                                                                 ,'v'    &
                                                                                 ,'T'    &
                                                                                 ,'delp' &
                                                                                 ,'DZ'   &
                                                                                 ,'W'    &  !<-- Not updated but needed in BC file.
                                                                                 /)
!
      character(len=5),dimension(num_fields_core_bc) :: varname_update_core_bc=(/       &
                                                                                 'u'    &
                                                                                ,'v'    &
                                                                                ,'t'    & 
                                                                                ,'delp' &
                                                                                ,'delz' &
                                                                                ,'w'    &
                                                                                /)
!
      character(len=50),dimension(:),allocatable :: varname_tracers_bc
!
      logical :: scalar_vbl
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Open the core restart file with combined interior and BC data 
!***  and get its ID.
!-----------------------------------------------------------------------
!
      call check(nf90_open(filename_core_combined,nf90_nowrite,ncid_core_combined))  !<-- The combined restart file's ID
!
!-----------------------------------------------------------------------
!***  Open the original core restart file.
!-----------------------------------------------------------------------
!
      call check(nf90_open(filename_core_restart,nf90_write,ncid_core_res))
!
!-----------------------------------------------------------------------
!***  Obtain the IDs of the combined core restart file's dimensions
!***  and read their values.
!-----------------------------------------------------------------------
!
      do n=1,ndims_res                                                     !<-- The # of dimensions in the core restart files.
        call check(nf90_inq_dimid(ncid_core_combined                    &  !<-- The combined core_restart file ID
                                 ,dimname_core(n)                       &  !<-- The dimension's name
                                 ,dimid_core(n)))                          !<-- The dimensions's ID
!
        call check(nf90_inquire_dimension(ncid_core_combined            &  !<-- The combined core restart file ID
                                         ,dimid_core(n)                 &  !<-- The dimensions's ID
                                         ,len=dimsize_combined(n)))        !<-- The dimension's value in the combined core file
      enddo
!
!-----------------------------------------------------------------------
!***  Subtract the halo from the combined restart file's dimensions.
!-----------------------------------------------------------------------
!
      dimsize_res(1)=dimsize_combined(1)-2*halo_integrate                  !<-- x dimension
      dimsize_res(2)=dimsize_res(1)+1                                      !<-- x+1 dimension
      dimsize_res(3)=dimsize_combined(3)-2*halo_integrate                  !<-- y+1 dimension
      dimsize_res(4)=dimsize_combined(4)-2*halo_integrate                  !<-- y dimension
      dimsize_res(5)=dimsize_combined(5)                                   !<-- z dimension
!
!-----------------------------------------------------------------------
!***  Open the normal BC file at the restart time.  The command 
!***  line argument is the forecast hour in the file name and
!***  cannot exceed 3 characters.
!-----------------------------------------------------------------------
!
      call get_command_argument(1,hour,length,istat)
      if(len_trim(hour)==0)then
        write(0,*)' User did not specify the BC hour in the command line!'
        write(0,*)' Aborting.'
        stop
      endif
      if(length>3)then
        write(0,*)' BC hour in command line must be <= 3 characters.'
        write(0,*)' Aborting.'
        stop
      endif
      if(istat/=0)then
        write(0,*)' Failed to read BC hour in command line, istat=',istat
        write(0,*)' Aborting.'
        stop
      endif
!
      hour=adjustr(hour)
      if(hour(1:1)==' ')hour(1:1)='0'
      if(hour(2:2)==' ')hour(2:2)='0'
!
      loc1=index(filename_bc,'hhh')
      filename_bc(loc1:loc1+2)=hour                                        !<-- Insert hour into the original BC filename
      loc1=index(filename_bc_new,'hhh')
      filename_bc_new(loc1:loc1+2)=hour                                    !<-- Insert hour into the new BC filename
!
      call check(nf90_open(filename_bc,nf90_write,ncid_bc))                !<-- Open the original BC file.
!
!-----------------------------------------------------------------------
!***  Find the number of blending rows that are in the current
!***  BC file.  The halo dimension is the sum of the number of
!***  boundary rows and blending rows.  Use the number of 
!***  blending rows when transferring data from the updated 
!***  enlarged restart files to the new BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_inq_dimid(ncid_bc                                 &  !<-- The original BC file's file ID
                               ,'halo'                                  &  !<-- The dimension's name
                               ,dimid_halo))                               !<-- The dimensions's ID
!
      call check(nf90_inquire_dimension(ncid_bc                         &  !<-- The original BC file's file ID
                                       ,dimid_halo                      &  !<-- The halo dimensions's ID
                                       ,len=width_halo_total))             !<-- # of boundary plus blending rows
!
      nrows_bndry=halo_integrate+1
      nrows_blend=width_halo_total-nrows_bndry                             !<-- # of blending rows in the original BC file
!
      write(0,101)nrows_blend
  101 format(' There are ',i3,' blending rows in the BC file.')
!
!-----------------------------------------------------------------------
!***  Create the new post-GSI BC file.
!-----------------------------------------------------------------------
!
      call create_new_bc_file(nrows_blend)                                 !<-- The new BC file with GSI update.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Begin reading fields from the updated core restart file.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Only one task is doing this work.  Do the reads layer-by-layer
!***  in order to avoid any memory issues.  Loop through the variables
!***  and read them from the restart file written by the DA.  We need
!***  to allocate arrays for the full field including boundary rows.
!-----------------------------------------------------------------------
!
      vbls_core: do n=1,num_fields_core_bc                                 !<-- The # of core restart fields written to the BC file.
!
!-----------------------------------------------------------------------
!
        call check(nf90_inq_varid(ncid_core_combined,varname_update_core(n),var_id_res))  !<-- Get this variable's ID.
!
        select case (varname_update_core(n))
!
          case ('u')
            ie_combined=dimsize_combined(1)
            je_combined=dimsize_combined(3)
!
          case ('v')
            ie_combined=dimsize_combined(2)
            je_combined=dimsize_combined(4)
!
          case ('T')
            ie_combined=dimsize_combined(1)
            je_combined=dimsize_combined(4)
!
          case ('delp')
            ie_combined=dimsize_combined(1)
            je_combined=dimsize_combined(4)

          case ('DZ')
            ie_combined=dimsize_combined(1)
            je_combined=dimsize_combined(4)
!
          case ('W')
            ie_combined=dimsize_combined(1)
            je_combined=dimsize_combined(4)
!
        end select
!
        kend=dimsize_res(5)
!
   write(0,*)' p1'
        if(varname_update_core(n)=='u'.or.varname_update_core(n)=='v')then
          scalar_vbl=.false.
          allocate(field_combined(1:ie_combined,1:je_combined))
        else
          scalar_vbl=.true.
          allocate(field_combined(0:ie_combined+1,0:je_combined+1))
        endif
!       
!       allocate(field_combined(1:ie_combined,1:je_combined))
        field_combined(:,:)=9.e9
!
   write(0,*)' p2'
!-----------------------------------------------------------------------
!***  Loop through the core restart field's layers. 
!-----------------------------------------------------------------------
!
        layers_core: do k=1,kend
!
          istart_res=1+halo_integrate                                      !<--
          iend_res  =ie_combined-halo_integrate                            !  These indices are 
          jstart_res=1+halo_integrate                                      !  recomputed for each k
          jend_res  =je_combined-halo_integrate                            !  because their values 
!                                                                          !  change in the boundary
          len_x=iend_res-istart_res+1                                      !  section below.
          len_y=jend_res-jstart_res+1                                      !<--
!
   write(0,*)' p2.1.1 vbl=',trim(varname_update_core(n)),' k=',k
   write(0,*)' istart_res=',istart_res,' iend_res=',iend_res,' jstart_res=',jstart_res,' jend_res=',jend_res
   write(0,*)' len_x=',len_x,' len_y=',len_y
   write(0,*)' ie_combined=',ie_combined,' je_combined=',je_combined
   write(0,*)' lbnd(field_combined)=',lbound(field_combined),' ubnd=',ubound(field_combined)
!-----------------------------------------------------------------------
!***  Read the field from the updated combined restart file.
!-----------------------------------------------------------------------
!
          call check(nf90_get_var(ncid_core_combined,var_id_res         &
                    ,field_combined(1:ie_combined,1:je_combined)        &  !<-- The full field with BC rows for layer k
                    ,start=(/1,1,k/)                                    &
                    ,count=(/ie_combined,je_combined,1/)))
!
   write(0,*)' p2.1 vbl=',trim(varname_update_core(n)),' k=',k
!-----------------------------------------------------------------------
!***  Write the interior data to the standard size core restart file.
!***  The GSI does not update W so skip writing it back into the
!***  original restart file although it will be written to the BC
!***  file.
!-----------------------------------------------------------------------
!
          if(trim(varname_update_core(n))/='W')then
            call check(nf90_put_var(ncid_core_res,var_id_res                                &
                                   ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,start=(/1,1,k/)                                         &
                                   ,count=(/len_x,len_y,1/)))
          endif
!
!-----------------------------------------------------------------------
!***  Insert the boundary rows back into the BC file.  There are
!***  three very important things to keep in mind. 
!
!***  (1) The wind components will be inserted into the BC files 
!***      oriented to the model grid and on model layers but the 
!***      winds in the BC files are oriented to geographic lat/lon
!***      on levels from the external forecast.  To handle this
!***      the model will be told to NOT reorient the BC winds or
!***      do any remapping when using the new BC file.
!
!***  (2) The regional FV3 has delp as a boundary array and that
!***      quantity is updated by the GSI so it will be copied 
!***      into the new BC file.  Although the interface height
!***      zh is in the original BC files the layer depth delz
!***      is a model boundary array so compute delz from the
!***      updated delp and also write it into the new BC file.
!***      We will use the hydrostatic approximation to derive 
!***      delz from delp.  We need the surface pressure to 
!***      integrate upward to get z on each side of the domain
!***      so simply sum delp as we move downward in k from the
!***      top. 
!
!***  (3) Blending data is simply an extension of boundary data.
!***      When writing data from the enlarged restart files into
!***      the new BC file then include the integration rows that
!***      correspond to blending rows.
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
   write(0,*)' p2.2 vbl=',trim(varname_update_core_bc(n)),' k=',k
          call get_bc_limits(varname_update_core_bc(n),'north',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res       &
                            ,istart_bc,jstart_bc,len_x,len_y               &
                            ,var_id_bc )
!
          if(trim(varname_update_core(n))=='delp')then
            if(.not.allocated(ps_north))then
              allocate(ps_north(istart_res-1:iend_res+1,jstart_res-1:jend_res))
              ps_north(:,:)=ptop
            endif
!
            do j=jstart_res,jend_res
            do i=istart_res,iend_res
              ps_north(i,j)=ps_north(i,j)+field_combined(i,j)              !<-- Pressure at bottom of layer k
            enddo
            enddo
!
            if(k==kend)then
!
              do i=istart_res,iend_res
                ps_north(i,jstart_res-1)=2.*ps_north(i,jstart_res)-ps_north(i,jstart_res+1)
              enddo
              do j=jstart_res-1,jend_res
                ps_north(istart_res-1,j)=2.*ps_north(istart_res,j)-ps_north(istart_res+1,j)
                ps_north(iend_res+1,j)=2.*ps_north(iend_res,j)-ps_north(iend_res-1,j)
              enddo
!
              call check(nf90_inq_varid(ncid_bc,'ps_bottom',var_id_ps))
              call check(nf90_put_var(ncid_bc_new,var_id_ps                             &  !<-- Write psfc into the new BC file
                                    ,ps_north(istart_res-1:iend_res+1,jstart_res-1:jend_res) &
                                    ,start=(/istart_bc-1,jstart_bc-1,k/)                   &
                                    ,count=(/len_x+2,len_y+1,1/)))
            endif
          endif
!
          if(scalar_vbl)then
            do i=istart_res,iend_res
              field_combined(i,jstart_res-1)=2.*field_combined(i,jstart_res)    &   !<-- The northern 4th BC row
                                               -field_combined(i,jstart_res+1)
            enddo
            do j=jstart_res-1,jend_res
               field_combined(istart_res-1,j)=2.*field_combined(istart_res,j)    &  !<-- The eastern 4th BC row
                                                -field_combined(istart_res+1,j)
               field_combined(iend_res+1,j)=2.*field_combined(iend_res,j)       &   !<-- The western 4th BC row
                                              -field_combined(iend_res-1,j)
            enddo
!
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &  !<-- Write field n into the new BC file
                                   ,field_combined(istart_res-1:iend_res+1,jstart_res-1:jend_res) &
                                   ,start=(/istart_bc-1,jstart_bc-1,k/)                         &
                                   ,count=(/len_x+2,len_y+1,1/)))
!
          else                                                                  !<-- Not adding outer (4th) BC row for winds
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &  !<-- Write field n into the new BC file
                                   ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,start=(/istart_bc,jstart_bc,k/)                         &
                                   ,count=(/len_x,len_y,1/)))
!
          endif
!
!-----------
!***  South
!-----------
!
          call get_bc_limits(varname_update_core_bc(n),'south',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res       &
                            ,istart_bc,jstart_bc,len_x,len_y               &
                            ,var_id_bc )
!
          if(trim(varname_update_core(n))=='delp')then
            if(.not.allocated(ps_south))then
              allocate(ps_south(istart_res-1:iend_res+1,jstart_res:jend_res+1))
              ps_south(:,:)=ptop
            endif
!
            do j=jstart_res,jend_res
            do i=istart_res,iend_res
              ps_south(i,j)=ps_south(i,j)+field_combined(i,j)              !<-- Pressure at bottom of layer k
            enddo
            enddo
!
            if(k==kend)then
!
              do i=istart_res,iend_res
                 ps_south(i,jend_res+1)=2.*ps_south(i,jend_res)-ps_south(i,jend_res-1)
              enddo
              do j=jstart_res,jend_res+1
                ps_south(istart_res-1,j)=2.*ps_south(istart_res,j)        &
                                           -ps_south(istart_res+1,j)
                ps_south(iend_res+1,j)=2.*ps_south(iend_res,j)            &
                                         -ps_south(iend_res-1,j)
              enddo
!
              call check(nf90_inq_varid(ncid_bc,'ps_top',var_id_ps))
              call check(nf90_put_var(ncid_bc_new,var_id_ps                             &
                                    ,ps_south(istart_res-1:iend_res+1,jstart_res:jend_res+1) &
                                    ,start=(/istart_bc-1,jstart_bc,k/)                   &
                                    ,count=(/len_x+2,len_y+1,1/)))
            endif
          endif
!
          if(scalar_vbl)then
            do i=istart_res,iend_res
              field_combined(i,jend_res+1)=2.*field_combined(i,jend_res)    &   !<-- The southern 4th BC row
                                             -field_combined(i,jend_res-1)
            enddo
            do j=jstart_res,jend_res+1
              field_combined(istart_res-1,j)=2.*field_combined(istart_res,j) &  !<-- The eastern 4th BC row
                                               -field_combined(istart_res+1,j)
              field_combined(iend_res+1,j)=2.*field_combined(iend_res,j)     &  !<-- The western 4th BC row
                                             -field_combined(iend_res-1,j)
            enddo
!
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res-1:iend_res+1,jstart_res:jend_res+1) &
                                   ,start=(/istart_bc-1,jstart_bc,k/)                         &
                                   ,count=(/len_x+2,len_y+1,1/)))
!
          else                                                                  !<-- Not adding outer (4th) BC row for winds
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,start=(/istart_bc,jstart_bc,k/)                         &
                                   ,count=(/len_x,len_y,1/)))
!
          endif
!
!----------
!***  East
!----------
!
          call get_bc_limits(varname_update_core_bc(n),'east ',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res       &
                            ,istart_bc,jstart_bc,len_x,len_y               &
                            ,var_id_bc )
!
          if(trim(varname_update_core(n))=='delp')then
            if(.not.allocated(ps_east))then
              allocate(ps_east(istart_res-1:iend_res,jstart_res:jend_res))
              ps_east(:,:)=ptop
            endif
!
            do j=jstart_res,jend_res
            do i=istart_res,iend_res
              ps_east(i,j)=ps_east(i,j)+field_combined(i,j)                !<-- Pressure at bottom of layer k
            enddo
            enddo
!
            if(k==kend)then
!
              do j=jstart_res,jend_res
                ps_east(istart_res-1,j)=2.*ps_east(istart_res,j)   &       !<-- The eastern 4th BC row
                                          -ps_east(istart_res+1,j)
              enddo
!
             call check(nf90_inq_varid(ncid_bc,'ps_left',var_id_ps))
             call check(nf90_put_var(ncid_bc_new,var_id_ps                            &
                                    ,ps_east(istart_res-1:iend_res,jstart_res:jend_res) &
                                    ,start=(/istart_bc-1,jstart_bc,k/)                  &
                                    ,count=(/len_x+1,len_y,1/)))
            endif
          endif
!
          if(scalar_vbl)then
            do j=jstart_res,jend_res
              field_combined(istart_res-1,j)=2.*field_combined(istart_res,j)   &   !<-- The eastern 4th BC row
                                               -field_combined(istart_res+1,j)
            enddo
!
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res-1:iend_res,jstart_res:jend_res) &
                                   ,start=(/istart_bc-1,jstart_bc,k/)                         &
                                   ,count=(/len_x+1,len_y,1/)))
!
          else                                                                  !<-- Not adding outer (4th) BC row for winds
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,start=(/istart_bc,jstart_bc,k/)                         &
                                   ,count=(/len_x,len_y,1/)))
         endif
!
!----------
!***  West
!----------
!
          call get_bc_limits(varname_update_core_bc(n),'west ',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res       &
                            ,istart_bc,jstart_bc,len_x,len_y               &
                            ,var_id_bc )
!
          if(trim(varname_update_core(n))=='delp')then
            if(.not.allocated(ps_west))then
              allocate(ps_west(istart_res:iend_res+1,jstart_res:jend_res))
              ps_west(:,:)=ptop
            endif
!
            do j=jstart_res,jend_res
            do i=istart_res,iend_res
              ps_west(i,j)=ps_west(i,j)+field_combined(i,j)                !<-- Pressure at bottom of layer k
            enddo
            enddo
!
            if(k==kend)then
!
              do j=jstart_res,jend_res
                ps_west(iend_res+1,j)=2.*ps_west(iend_res,j)    &
                                        -ps_west(iend_res-1,j)
              enddo
!
             call check(nf90_inq_varid(ncid_bc,'ps_right',var_id_ps))
             call check(nf90_put_var(ncid_bc_new,var_id_ps                                   &
                                    ,ps_west(istart_res:iend_res+1,jstart_res:jend_res) &
                                    ,start=(/istart_bc,jstart_bc,k/)                         &
                                    ,count=(/len_x+1,len_y,1/)))
            endif
          endif
!
          if(scalar_vbl)then
            do j=jstart_res,jend_res
              field_combined(iend_res+1,j)=2.*field_combined(iend_res,j)    &   !<-- The western 4th BC row
                                             -field_combined(iend_res-1,j)
            enddo
!
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res:iend_res+1,jstart_res:jend_res) &
                                   ,start=(/istart_bc,jstart_bc,k/)                         &
                                   ,count=(/len_x+1,len_y,1/)))
!
          else                                                                  !<-- Not adding outer (4th) BC row for winds
            call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                   ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,start=(/istart_bc,jstart_bc,k/)                         &
                                   ,count=(/len_x,len_y,1/)))
!
          endif
!
!-----------------------------------------------------------------------
!***  The updated D-grid wind has been moved into the BC file from 
!***  the enlarged core restart file but the C-grid winds are needed
!***  there too.  Interpolate from the D-grid wind and then write the
!***  C-grid wind into the BC file for this model layer.
!-----------------------------------------------------------------------
!
          if(trim(varname_update_core_bc(n))=='u')then
            call dgrid_to_cgrid('u')
          elseif(trim(varname_update_core_bc(n))=='v')then
            call dgrid_to_cgrid('v')
          endif
!
!-----------------------------------------------------------------------
!
        enddo layers_core
!
        deallocate(field_combined)
!
!-----------------------------------------------------------------------
!
      enddo vbls_core
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  The updated D-grid u, D-grid v, and T have been written back
!***  into the standard sized core restart file containing only the
!***  domain's interior points and they have been written into the
!***  BC file.  Now repeat the process with the tracer restart file.
!***  The specific humidity is the only variable updated by the GSI.
!***  However the updated BC file will now be holding variables on
!***  the model (not input) layers so the remaining contents of
!***  the BC file will also need to be within model layers.  Take
!***  all boundary tracers from the updated tracer restart file
!***  and write them into the BC file.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      call check(nf90_open(filename_tracer_combined,nf90_nowrite,ncid_tracer_combined))  !<-- The combined tracer file's ID
      call check(nf90_open(filename_tracer_restart,nf90_write,ncid_tracer_res))          !<-- The original tracer file's ID
!
!-----------------------------------------------------------------------
!***  Find the number of tracers and their names in the tracer
!***  restart files.
!-----------------------------------------------------------------------
!
      call tracer_info(ncid_tracer_res,ncid_bc                          &
                      ,num_tracers_bc,varname_tracers_bc)
!
!-----------------------------------------------------------------------
!***  Only one task is doing this work.  Do the reads layer-by-layer
!***  in order to avoid any memory issues.  Loop through the variables
!***  and read them from the tracer restart file written by the DA.  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      vbls_tracer: do n=1,num_tracers_bc
!
!-----------------------------------------------------------------------
!
        call check(nf90_inq_varid(ncid_tracer_combined,varname_tracers_bc(n),var_id_tracer))  !<-- Get this tracer's ID.
!
        ie_combined=dimsize_combined(1)
        je_combined=dimsize_combined(4)
!
        kend=dimsize_res(5)
!
        if(.not.allocated(field_combined))then
!xxx      allocate(field_combined(1:ie_combined,1:je_combined))
          allocate(field_combined(1:ie_combined+2,1:je_combined+2))        !<-- Add the outermost boundary row
        endif
        field_combined(:,:)=9.e9
!
!-----------------------------------------------------------------------
!***  Loop through the layers in the tracer restart file.  
!-----------------------------------------------------------------------
!
        layers_tracer: do k=1,kend
!
          istart_res=1+halo_integrate
          iend_res  =ie_combined-halo_integrate
          jstart_res=1+halo_integrate
          jend_res  =je_combined-halo_integrate
!
          len_x=iend_res-istart_res+1
          len_y=jend_res-jstart_res+1
!
!-----------------------------------------------------------------------
!***  Read the tracer from the updated combined tracer file.
!-----------------------------------------------------------------------
!
          call check(nf90_get_var(ncid_tracer_combined,var_id_tracer    &
!                   ,field_combined(:,:)                                &  !<-- The full field with BC rows for layer k
                    ,field_combined(2:ie_combined+1,2:je_combined+1)    &  !<-- The full field with BC rows for layer k
                    ,start=(/1,1,k/)                                    &
                    ,count=(/ie_combined,je_combined,1/)))
!
!-----------------------------------------------------------------------
!***  Write the interior data to the standard size tracer restart file
!***  for the only updated quantity (sphum).
!-----------------------------------------------------------------------
!
          if(trim(varname_tracers_bc(n))=='sphum')then
            call check(nf90_put_var(ncid_tracer_res,var_id_tracer                           &
!                                  ,field_combined(istart_res:iend_res,jstart_res:jend_res) &
                                   ,field_combined(istart_res+1:iend_res+1,jstart_res+1:jend_res+1) &
                                   ,start=(/1,1,k/)                                         &
                                   ,count=(/len_x,len_y,1/)))
          endif
!
!
!-----------------------------------------------------------------------
!***  Insert the boundary rows of each tracer back into the BC file
!***  on each side of the domain.  These data are on the forecast 
!***  model layers whereas in the original BC file the data is on
!***  the input model layers.
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
          call get_bc_limits(trim(varname_tracers_bc(n)),'north',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res         &
                            ,istart_bc,jstart_bc,len_x,len_y                 &
                            ,var_id_bc )
!
          do i=2,ie_combined+1
            field_combined(i,1)=2.*field_combined(i,2)-field_combined(i,3)
          enddo
          do j=jstart_bc-1,jstart_bc+len_y-1
            field_combined(1,j)=2.*field_combined(2,j)-field_combined(3,j)
            field_combined(ie_combined+2,j)=2.*field_combined(ie_combined+1,j) &
                                              -field_combined(ie_combined,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                 ,field_combined(1:ie_combined+2,1:len_y+1)               &
                                 ,start=(/1,1,k/)                                         &
                                 ,count=(/len_x+2,len_y+1,1/)))
!
!-----------
!***  South
!-----------
!
          call get_bc_limits(trim(varname_tracers_bc(n)),'south',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res         &
                            ,istart_bc,jstart_bc,len_x,len_y                 &
                            ,var_id_bc )
!

          do i=2,ie_combined+1
            field_combined(i,je_combined+2)=2.*field_combined(i,je_combined+1) &
                                              -field_combined(i,je_combined)
          enddo
          do j=jstart_res+1,jend_res+2
            field_combined(1,j)=2.*field_combined(2,j)-field_combined(3,j)
            field_combined(ie_combined+2,j)=2.*field_combined(ie_combined+1,j) &
                                              -field_combined(ie_combined,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_bc                                      &
                                 ,field_combined(1:ie_combined+2,jstart_res+1:je_combined+2) &
                                 ,start=(/1,1,k/)                                            &
                                 ,count=(/len_x+2,len_y+1,1/)))
!
!----------
!***  East
!----------
!
          call get_bc_limits(trim(varname_tracers_bc(n)),'east ',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res         &
                            ,istart_bc,jstart_bc,len_x,len_y                 &
                            ,var_id_bc )
!
          do j=jstart_res+1,jend_res+1
            field_combined(1,j)=2.*field_combined(2,j)-field_combined(3,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_bc                                &
                                 ,field_combined(1:len_x+1,jstart_res+1:jstart_res+len_y) &
                                 ,start=(/1,1,k/)                                      &
                                 ,count=(/len_x+1,len_y,1/)))
!
!----------
!***  West
!----------
!
          call get_bc_limits(trim(varname_tracers_bc(n)),'west ',nrows_blend &
                            ,istart_res,iend_res,jstart_res,jend_res         &
                            ,istart_bc,jstart_bc,len_x,len_y                 &
                            ,var_id_bc )
!
          do j=jstart_res+1,jend_res+1
            field_combined(ie_combined+2,j)=2.*field_combined(ie_combined+1,j) &
                                              -field_combined(ie_combined,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_bc                                   &
                                 ,field_combined(istart_res+1:iend_res+2,jstart_res+1:jstart_res+len_y) &
                                 ,start=(/1,1,k/)                                         &
                                 ,count=(/len_x+1,len_y,1/)))
!
!-----------------------------------------------------------------------
!
        enddo layers_tracer
!
!-----------------------------------------------------------------------
!
      enddo vbls_tracer
!
      deallocate(field_combined)
!
      if(allocated(row4_north))then
        deallocate(row4_north)
      endif
      if(allocated(row4_south))then
        deallocate(row4_south)
      endif
      if(allocated(row4_east))then
        deallocate(row4_east)
      endif
      if(allocated(row4_west))then
        deallocate(row4_west)
      endif
!
!-----------------------------------------------------------------------
!***  To get the layer thickness (m) needed in the BC file we need
!***  the sensible temperature, delta p, specific humidity, and
!***  surface elevation to integrate upward hydrostatically.
!-----------------------------------------------------------------------
!
      call check(nf90_inq_varid(ncid_tracer_combined,'sphum',var_id_sphum))
      call check(nf90_inq_varid(ncid_core_combined,'T',var_id_T))
      call check(nf90_inq_varid(ncid_core_combined,'delp',var_id_delp))
!
!-----------------------------------------------------------------------
!***  Allocate phis and then read it from the combined restart file
!***  so we can use it to start integrating z upward.
!-----------------------------------------------------------------------
!
      allocate(phis_combined(1:dimsize_combined(1),1:dimsize_combined(4)))
!
      call check(nf90_inq_varid(ncid_core_combined,'phis',var_id_phis))
!
      call check(nf90_get_var(ncid_core_combined,var_id_phis        &
                ,phis_combined(:,:)                                 &
                ,start=(/1,1/)                                      &
                ,count=(/dimsize_combined(1),dimsize_combined(4)/)))
!
!-----------------------------------------------------------------------
!***  We need to insert layer thicknesses in the boundary rows
!***  based on the DA updates of delp and T in the core restart
!***  file and sphum in the core tracer file.  Extract the boundary
!***  points of those variables and integrate z upward from the
!***  surface elevation.  Write delz into the BC file.
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
!x      call write_delz_to_bc_file('north',nrows_blend)
!
!-----------
!***  South
!-----------
!
!x      call write_delz_to_bc_file('south',nrows_blend)
!
!----------
!***  East
!----------
!
!x      call write_delz_to_bc_file('east ',nrows_blend)
!
!----------
!***  West
!----------
!
!x      call write_delz_to_bc_file('west ',nrows_blend)
!
!-----------------------------------------------------------------------
!
      deallocate(ps_north)
      deallocate(ps_south)
      deallocate(ps_east)
      deallocate(ps_west)
      deallocate(phis_combined)
!
!-----------------------------------------------------------------------
!
      call check(nf90_close(ncid_core_res))
      call check(nf90_close(ncid_core_combined))
      call check(nf90_close(ncid_tracer_res))
      call check(nf90_close(ncid_tracer_combined))
      call check(nf90_close(ncid_bc))
      call check(nf90_close(ncid_bc_new))
!
!-----------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------
!
      subroutine check(status)
!
      integer,intent(in) :: status
!
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
!
      end subroutine check
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine create_new_bc_file(nrows_blend)
!
!-----------------------------------------------------------------------
!***  Create a new BC file and prepare its dimensions and variables.
!***  The number of layers is one less than in the original BC file
!***  since the top dummy layer is removed.  All fields will be on
!***  the forecast model layers and not input model layers.
!-----------------------------------------------------------------------
!
      integer,intent(in) :: nrows_blend                                    !<-- # of blending rows
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: halo_bc=halo_integrate+1
      integer :: var_id,var_id_new
!
      integer,dimension(1:3) :: dimids=(/0,0,0/)                        &
                               ,dimids_north=(/0,0,0/)                  &
                               ,dimids_south=(/0,0,0/)                  &
                               ,dimids_east =(/0,0,0/)                  &
                               ,dimids_west =(/0,0,0/) 
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create the new BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_create(filename_bc_new,nf90_clobber,ncid_bc_new))
!
!-----------------------------------------------------------------------
!***  The variables' dimensions.
!-----------------------------------------------------------------------
!
      dimsize_bc(1)=dimsize_res(1)+2*halo_bc
      dimsize_bc(2)=dimsize_res(4)
      dimsize_bc(3)=dimsize_bc(1)+1
      dimsize_bc(4)=dimsize_bc(2)-1
      dimsize_bc(5)=halo_bc+nrows_blend
      dimsize_bc(6)=dimsize_bc(5)+1
      dimsize_bc(7)=dimsize_res(5)
      dimsize_bc(8)=dimsize_bc(7)+1
!
!-----------------------------------------------------------------------
!***  Define the dimensions in the new BC file.
!-----------------------------------------------------------------------
!
      do n=1,ndims_bc
        call check(nf90_def_dim(ncid_bc_new                             &
                               ,dimname_bc(n)                           &
                               ,dimsize_bc(n)                           &
                               ,dimid_bc(n)))
      enddo
!
!-----------------------------------------------------------------------
!***  How many variables are in a normal BC file?  Get that number
!***  then loop through the variables in a normal BC file and define
!***  those in the new post-GSI BC file.  The exception is that zh
!***  will be skipped since it will be replaced by delz.
!-----------------------------------------------------------------------
!
      call check(nf90_inquire(ncid_bc                                   &
                             ,nvariables=num_vars_bc))                     !<-- Total # of vbls in a normal BC file.
!
!-----------------------------------------------------------------------
!***  The new file's variables must be defined while that file is
!***  in define mode.  Define each variable in the new file using
!***  those in the original file but exclude zh since instead we
!***  we will want delz.
!-----------------------------------------------------------------------
!
      var_id_new=0
!
      do n=1,num_vars_bc
        var_id=n
        call check(nf90_inquire_variable(ncid_bc,var_id,name_var,nctype &  !<-- Name and type of this variable
                  ,ndims,dimids,natts))                                    !<-- # of dimensions and attributes in this variable
!
        if(name_var(1:2)/='zh')then                                        !<-- We do not need zh in the new BC file
          var_id_new=var_id_new+1
          call check(nf90_def_var(ncid_bc_new                           &
                                 ,name_var                              &  !<-- The variable's name
                                 ,nctype                                &  !<-- The variable's type
                                 ,dimids(1:ndims)                       &  !<-- The IDs of the variable's dimensions
                                 ,var_id_new))                             !<-- The variable's ID
!
!-----------------------------------------------------------------------
!***  Copy each variable's attributes to the new file's variable.
!-----------------------------------------------------------------------
!
          if(natts>0)then
            do na=1,natts
              call check(nf90_inq_attname(ncid_bc,var_id,na,name_att))     !<-- Get the attribute's name and ID from restart.
              call check(nf90_copy_att(ncid_bc                          &  
                                      ,var_id                           & 
                                      ,name_att                         &
                                      ,ncid_bc_new                      &
                                      ,var_id_new))  
            enddo
          endif
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!***  We want to add delp and delz as variables in the new BC file. 
!***  In the original regional FV3 the boundary values for delp are 
!***  computed during the vertical remapping of the input data.  In 
!***  the DA process the remapping and wind rotation are bypassed so 
!***  the delp and delz values must be present in the BC file to fill
!***  the BC arrays in the model.  Use the same dimensions as t.
!-----------------------------------------------------------------------
!
!-----------------------
!***  Delp on all sides
!-----------------------
!
      call check(nf90_inq_varid(ncid_bc,'t_bottom',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_north))    !<-- Use T's dimension IDs
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_bottom'                             &  !<-- Define delp on the domain bottom (north)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_north(1:ndims)                     &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp bottom bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_top',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_south))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_top'                                &  !<-- Define delp on the domain top (south)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_south(1:ndims)                     &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp top bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_left',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_east))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_left'                               &  !<-- Define delp on the domain's left (east)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_east(1:ndims)                      &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp left bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
      call check(nf90_inq_varid(ncid_bc,'t_right',var_id_t))
      call check(nf90_inquire_variable(ncid_bc,var_id_t,dimids=dimids_west))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delp_right'                              &  !<-- Define delp on the domain's right (west)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_west(1:ndims)                      &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delp right bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","pascals"))
!
!-----------------------
!***  Delz on all sides
!-----------------------
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_bottom'                             &  !<-- Define delz on the domain bottom (north)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_north(1:ndims)                     &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz bottom bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_top'                                &  !<-- Define delz on the domain top (south)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_south(1:ndims)                     &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz top bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_left'                               &  !<-- Define delz on the domain's left (east)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_east(1:ndims)                      &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz left bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
      call check(nf90_def_var(ncid_bc_new                               &
                             ,'delz_right'                              &  !<-- Define delz on the domain's right (west)
                             ,NF90_FLOAT                                &  !<-- The variable's type
                             ,dimids_west(1:ndims)                      &  !<-- The IDs of the variable's dimensions
                             ,var_id))                                     !<-- The variable's ID
      call check(nf90_put_att(ncid_bc_new,var_id,"long_name","delz right bndy"))
      call check(nf90_put_att(ncid_bc_new,var_id,"units","m"))
!
!-----------------------------------------------------------------------
!***  Find the number of global attributes in the original BC file
!***  and copy them to the new BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_inquire(ncid_bc,nattributes=ngatts))
!
      do n=1,ngatts
        call check(nf90_inq_attname(ncid_bc,NF90_GLOBAL,n,name_att))
        call check(nf90_copy_att(ncid_bc                                &
                                ,NF90_GLOBAL                            &
                                ,name_att                               &
                                ,ncid_bc_new                            &
                                ,NF90_GLOBAL))
      enddo
!
!-----------------------------------------------------------------------
!***  We are finished with definitions so put the new BC file
!***  into data mode.
!-----------------------------------------------------------------------
!
      call check(nf90_enddef(ncid_bc_new))
!
!-----------------------------------------------------------------------
!
      end subroutine create_new_bc_file
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine get_bc_limits(field,side,nrows_blend                   &
                              ,istart_res,iend_res,jstart_res,jend_res  &
                              ,istart_bc,jstart_bc,len_x,len_y          &
                              ,var_id_bc )
!
!-----------------------------------------------------------------------
!***  Determine the index limits of the boundary variables.  These
!***  values will include blending rows which are extensions of the
!***  actual boundary rows.
!-----------------------------------------------------------------------
!
      integer,intent(in) :: nrows_blend                                    !<-- # of blending rows
!
      character(len=5),intent(in) :: side
      character(len=*),intent(in) :: field
!
      integer,intent(out) :: istart_res,iend_res,jstart_res,jend_res    &
                            ,istart_bc,jstart_bc,len_x,len_y            &
                            ,var_id_bc
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if(trim(field)=='u')then
!
        if(side=='north')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=1                                     !   including the BC/integration line.
          jend_res=halo_integrate+nrows_blend+1            !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(1)+1                    !   Data limits in the BC file array
          jstart_bc=2                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend     !<--
!
          varname_update_bc=trim(field)//'_s_bottom'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='south')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=dimsize_combined(3)-halo_integrate-nrows_blend    !   including the BC/integration line.
          jend_res=dimsize_combined(3)                     !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(1)+1                    !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend     !<--
!
          varname_update_bc=trim(field)//'_s_top'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='east')then
          istart_res=1                                     !<--
          iend_res=halo_integrate+nrows_blend              !   Data limits in the combined array
          jstart_res=halo_integrate+2                      !   including the BC/integration line.
          jend_res=dimsize_combined(3)-4                   !<--
!
          istart_bc=2                                      !<--
          iend_bc=halo_integrate+nrows_blend+1             !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(3)-2*halo_integrate-2   !<--
!
          varname_update_bc=trim(field)//'_s_left'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='west')then
          istart_res=dimsize_combined(1)-halo_integrate-nrows_blend+1  !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=halo_integrate+2                      !   including the BC/integration line.
          jend_res=dimsize_combined(3)-4                   !<--
!
          istart_bc=1                                      !<--
          iend_bc=istart_bc+halo_integrate+nrows_blend-1   !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(3)-2*halo_integrate-2   !<--
!
          varname_update_bc=trim(field)//'_s_right'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        endif
!
      elseif(trim(field)=='v')then
!
        if(side=='north')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(2)                     !   Data limits in the combined array
          jstart_res=1                                     !   including the BC/integration line.
          jend_res=halo_integrate+nrows_blend              !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(2)+1                    !   Data limits in the BC file array
          jstart_bc=2                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend-1   !<--
!
          varname_update_bc=trim(field)//'_w_bottom'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='south')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(2)                     !   Data limits in the combined array
          jstart_res=dimsize_combined(4)-nrows_blend-2     !   including the BC/integration line.
          jend_res=dimsize_combined(4)                     !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(2)+1                    !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend-1   !<--
!
          varname_update_bc=trim(field)//'_w_top'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='east')then
          istart_res=1                                     !<--
          iend_res=halo_integrate+nrows_blend+1            !   Data limits in the combined array
          jstart_res=halo_integrate+1                      !   including the BC/integration line.
          jend_res=dimsize_combined(4)-halo_integrate      !<--
!
          istart_bc=2                                      !<--
          iend_bc=halo_integrate+nrows_blend+2             !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(4)-2*halo_integrate     !<--
!
          varname_update_bc=trim(field)//'_w_left'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='west')then
          istart_res=dimsize_combined(2)-halo_integrate-nrows_blend    !<--
          iend_res=dimsize_combined(2)                     !   Data limits in the combined array
          jstart_res=halo_integrate+1                      !   including the BC/integration line.
          jend_res=dimsize_combined(4)-halo_integrate      !<--
!
          istart_bc=1                                      !<--
          iend_bc=istart_bc+halo_integrate+nrows_blend     !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(4)-2*halo_integrate     !<--
!
          varname_update_bc=trim(field)//'_w_right'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        endif
!
      else                                                 !<-- All mass-point variables
!
        if(side=='north')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=1                                     !   including the BC/integration line.
          jend_res=halo_integrate+nrows_blend              !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(1)+1                    !   Data limits in the BC file array
          jstart_bc=2                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend-1   !<--
!
          varname_update_bc=trim(field)//'_bottom'
          call check(nf90_inq_varid(ncid_bc_new,trim(varname_update_bc),var_id_bc))
!
        elseif(side=='south')then
          istart_res=1                                     !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=dimsize_combined(4)-halo_integrate-nrows_blend+1  !   including the BC/integration line.
          jend_res=dimsize_combined(4)                     !<--
!
          istart_bc=2                                      !<--
          iend_bc=dimsize_combined(1)+1                    !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=jstart_bc+halo_integrate+nrows_blend-1   !<--
!
          varname_update_bc=trim(field)//'_top'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='east')then
          istart_res=1                                     !<--
          iend_res=halo_integrate+nrows_blend              !   Data limits in the combined array
          jstart_res=halo_integrate+1                      !   including the BC/integration line.
          jend_res=dimsize_combined(4)-halo_integrate      !<--
!
          istart_bc=2                                      !<--
          iend_bc=halo_integrate+nrows_blend+1             !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(4)-2*halo_integrate     !<--
!
          varname_update_bc=trim(field)//'_left'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        elseif(side=='west')then
          istart_res=dimsize_combined(1)-halo_integrate-nrows_blend+1  !<--
          iend_res=dimsize_combined(1)                     !   Data limits in the combined array
          jstart_res=halo_integrate+1                      !   including the BC/integration line.
          jend_res=dimsize_combined(4)-halo_integrate      !<--
!
          istart_bc=1                                      !<--
          iend_bc=istart_bc+halo_integrate+nrows_blend-1   !   Data limits in the BC file array
          jstart_bc=1                                      !   including the BC/integration line.
          jend_bc=dimsize_combined(4)-2*halo_integrate     !<--
!
          varname_update_bc=trim(field)//'_right'
          call check(nf90_inq_varid(ncid_bc_new,varname_update_bc,var_id_bc))
!
        endif
!
      endif
!
      len_x=iend_bc-istart_bc+1
      len_y=jend_bc-jstart_bc+1
!
!-----------------------------------------------------------------------
!
      end subroutine get_bc_limits
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine tracer_info(ncid_tracer_res,ncid_bc                    &
                            ,num_tracers_bc,varname_tracers_bc)
!
!-----------------------------------------------------------------------
!***  Determine the number of tracers in the BC files and save
!***  their names.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: ncid_tracer_res                             &  !<-- File ID of the normal tracer file.
                           ,ncid_bc                                        !<-- File ID of the current BC file.
!
      character(len=50),dimension(:),allocatable,intent(inout) ::       &
                                                    varname_tracers_bc     !<-- Names of the tracers in the BC files.
!
      integer,intent(out) :: num_tracers_bc                                !<-- The # of tracers in the BC files.
!
!---------------------
!***  Local variables
!---------------------
!
      integer :: n,nn,num_bc_vbls,num_tracer_vbls,num_tracers_res
!
      character(len=50) :: tracer_name
!
      character(len=50),dimension(:),allocatable :: bc_vbl_names
!
      type linked_list
        character(len=50) :: save_name
        type(linked_list),pointer :: next_link
      end type linked_list
!
      type(linked_list),pointer :: hold_bc_tracer_name                  &
                                  ,head,tail
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create a list of all variables in the BC file.
!-----------------------------------------------------------------------
!
      call check(nf90_inquire(ncid_bc                                   &  
                             ,nvariables=num_bc_vbls))                     !<-- Total # of tracers in BC file
!
      allocate(bc_vbl_names(1:num_bc_vbls)) 
!
      do n=1,num_bc_vbls 
        call check(nf90_inquire_variable(ncid_bc,n                      &
                                        ,name=bc_vbl_names(n)))            !<-- Save the BC file vbl names
      enddo
!
!-----------------------------------------------------------------------
!***  Extract the tracer names from the tracer restart file and
!***  compare them against the variable names in the BC file.
!***  The ones that match are the tracers in the BC file.
!***  We do not know how many BC tracers there will be so we
!***  must use a linked list to store the names as we find them.
!-----------------------------------------------------------------------
!
      call check(nf90_inquire(ncid_tracer_res                           &  
                             ,nvariables=num_tracer_vbls))                 !<-- Total # of tracers in tracer restart file
      num_tracers_res=num_tracer_vbls-4                                    !<-- The 1st 4 variables are dimensions.
!
      num_tracers_bc=0
      tail=>null()
!
      outer: do n=5,num_tracer_vbls 
        call check(nf90_inquire_variable(ncid_tracer_res,n              &
                                        ,name=tracer_name))                !<-- Tracer name in the tracer restart file
!        
        do nn=1,num_bc_vbls
          if(index(trim(bc_vbl_names(nn)),trim(tracer_name))/=0)then       !<-- If true we found a tracer in the BC file
            num_tracers_bc=num_tracers_bc+1
            if(num_tracers_bc==1)then
              allocate(hold_bc_tracer_name)
              tail=>hold_bc_tracer_name                                    !<-- Point at the top of the linked list
              nullify(hold_bc_tracer_name%next_link)
            else
              allocate(tail%next_link)
              tail=>tail%next_link                                         !<-- Point at the next link
              nullify(tail%next_link)
            endif
            tail%save_name=tracer_name                                     !<-- Save the BC tracer name in this link.
            cycle outer
          endif
        enddo
!
      enddo outer
!
!-----------------------------------------------------------------------
!***  Now dereference the linked list and store the BC tracer names
!***  in an array.
!-----------------------------------------------------------------------
!
      allocate(varname_tracers_bc(1:num_tracers_bc))
!
      n=0
      head=>hold_bc_tracer_name
!
      dealloc: do 
        n=n+1
        if(n>num_tracers_bc)then
          exit dealloc
        endif
        tail=>null()
        varname_tracers_bc(n)=head%save_name
        if(associated(head%next_link))then
          tail=>head%next_link
        endif
!
        if(n>1)then
          deallocate(head,stat=istat)
          if(istat/=0)then
            write(0,10001)n
10001       format(' Failed to deallocate link ',i3) 
          endif
        endif
!
        if(associated(tail))then
          head=>tail
        else
          exit dealloc
        endif
!
      enddo dealloc
!
      end subroutine tracer_info
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine write_delz_to_bc_file(side,nrows_blend)
!
!-----------------------------------------------------------------------
!***  Given delp and T from the updated combined core restart file
!***  and sphum from the tracer file compute the layer thicknesses.
!***  and write them to the new BC file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      integer,intent(in) :: nrows_blend                                     !<-- # of blending rows
!
      character(len=5),intent(in) :: side                                   !<-- Current side of the domain
!
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,j,k
      integer :: istart_res,iend_res,jstart_res,jend_res     
      integer :: istart_bc,jstart_bc,len_x,len_y
      integer :: ie_combined,je_combined
      integer :: var_id_sphum_bc,var_id_delz_bc
!
      real :: pmid
!
      real,dimension(:,:),allocatable :: delp_side                      &
                                        ,delz_side                      &
                                        ,sphum_side                     &
                                        ,T_side                         &
                                        ,zh_bot                         &
                                        ,zh_top
!
      real,dimension(:),allocatable :: row4_east,row4_north             &
                                      ,row4_south,row4_west
!     real,dimension(:,:),pointer :: pint
      real,dimension(:,:),allocatable :: pint
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Find the index limits and the ID for zh on the current side
!***  of the domain.  The '_res' indices are on the enlarged restart
!***  file's grid that are valid for the relevant variable's boundary
!***  rows.
!-----------------------------------------------------------------------
!
      call get_bc_limits('delz',side,nrows_blend                        &
                        ,istart_res,iend_res,jstart_res,jend_res        &
                        ,istart_bc,jstart_bc,len_x,len_y                &
                        ,var_id_delz_bc)
!
      ie_combined=iend_res-istart_res+1
      je_combined=jend_res-jstart_res+1
!
      allocate(sphum_side(istart_res:iend_res,jstart_res:jend_res))
      allocate(T_side    (istart_res:iend_res,jstart_res:jend_res))
      allocate(delp_side (istart_res:iend_res,jstart_res:jend_res))
!     allocate(delz_side (istart_res:iend_res,jstart_res:jend_res))
      allocate(zh_bot    (istart_res:iend_res,jstart_res:jend_res))
      allocate(zh_top    (istart_res:iend_res,jstart_res:jend_res))
      allocate(pint      (istart_res:iend_res,jstart_res:jend_res))
!
      if(side=='north')then
        allocate(delz_side(istart_res-1:iend_res+1,jstart_res-1:jend_res))
      elseif(side=='south')then
        allocate(delz_side(istart_res-1:iend_res+1,jstart_res:jend_res+1))
      elseif(side=='east ')then
        allocate(delz_side(istart_res-1:iend_res,jstart_res:jend_res))
      elseif(side=='west ')then
        allocate(delz_side(istart_res:iend_res+1,jstart_res:jend_res))
      endif
!-----------------------------------------------------------------------
!***  The elevation (m) of the ground surface.
!-----------------------------------------------------------------------
!
      do j=jstart_res,jend_res
      do i=istart_res,iend_res
        zh_bot(i,j)=rgrav*phis_combined(i,j)
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***  Obtain the mid-layer and interface pressures by beginning with
!***  the pressure at the ground surface.
!-----------------------------------------------------------------------
!
      if(trim(side)=='north')then
!       pint=>ps_north
!       pint(:,:)=ps_north
        do j=jstart_res,jend_res
        do i=istart_res,iend_res
          pint(i,j)=ps_north(i,j)
        enddo
        enddo
      elseif(trim(side)=='south')then
!       pint=>ps_south
!       pint(:,:)=ps_south
        do j=jstart_res,jend_res
        do i=istart_res,iend_res
          pint(i,j)=ps_south(i,j)
        enddo
        enddo
      elseif(trim(side)=='east')then
!       pint=>ps_east
!       pint(:,:)=ps_east
        do j=jstart_res,jend_res
        do i=istart_res,iend_res
          pint(i,j)=ps_east(i,j)
        enddo
        enddo
      elseif(trim(side)=='west')then
!       pint=>ps_west
!       pint(:,:)=ps_west
        do j=jstart_res,jend_res
        do i=istart_res,iend_res
          pint(i,j)=ps_west(i,j)
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
      do k=kend,1,-1
!-----------------------------------------------------------------------
!
        call check(nf90_get_var(ncid_tracer_combined,var_id_sphum       &
                  ,sphum_side(:,:)                                      &
                  ,start=(/istart_res,jstart_res,k/)                    &
                  ,count=(/ie_combined,je_combined,1/)))
!
        call check(nf90_get_var(ncid_core_combined,var_id_T             &
                  ,T_side(:,:)                                          &
                  ,start=(/istart_res,jstart_res,k/)                    &
                  ,count=(/ie_combined,je_combined,1/)))
!
        call check(nf90_get_var(ncid_core_combined,var_id_delp          &
                  ,delp_side(:,:)                                       &
                  ,start=(/istart_res,jstart_res,k/)                    &
                  ,count=(/ie_combined,je_combined,1/)))
!
        do j=jstart_res,jend_res
        do i=istart_res,iend_res
          pmid=pint(i,j)-0.5*delp_side(i,j)
          zh_top(i,j)=zh_bot(i,j)-rd*(1.+f608*sphum_side(i,j))*T_side(i,j)        &
                                  *(log(pint(i,j)-delp_side(i,j))-log(pint(i,j))) &
                                  *rgrav
          delz_side(i,j)=zh_bot(i,j)-zh_top(i,j)
          if(delz_side(i,j)>0.)then
            write(0,10101)i,j,k,delz_side(i,j),zh_top(i,j),zh_bot(i,j),trim(side)
            write(0,*)' ABORTING'
10101       format(' BAD DELZ(',i3,',',i3,',',i3,')=',e12.5,' zh_top=',e12.5,' zh_bot=',e12.5,' side=',a)
            stop
          endif
          pint(i,j)=pint(i,j)-delp_side(i,j)
          zh_bot(i,j)=zh_top(i,j)
        enddo
        enddo
!
!-----------------------------------------------------------------------
!***  Add the outermost boundary row then put all boundary points
!***  of delz into the new BC file.
!-----------------------------------------------------------------------
!
        if(side=='north')then
          do i=istart_res,iend_res
            delz_side(i,jstart_res-1)=2.*delz_side(i,jstart_res)        &
                                        -delz_side(i,jstart_res+1)
          enddo
          do j=jstart_res-1,jend_res
            delz_side(istart_res-1,j)=2.*delz_side(istart_res,j)        &
                                        -delz_side(istart_res+1,j)
            delz_side(iend_res+1,j)=2.*delz_side(iend_res,j)            &
                                      -delz_side(iend_res-1,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_delz_bc                         &
                                 ,delz_side(istart_res-1:iend_res+1,jstart_res-1:jend_res) &
                                 ,start=(/istart_bc-1,jstart_bc-1,k/)                    &
                                 ,count=(/len_x+2,len_y+1,1/)))
!
        elseif(side=='south')then
          do i=istart_res,iend_res
            delz_side(i,jend_res+1)=2.*delz_side(i,jend_res)            &
                                        -delz_side(i,jend_res-1)
          enddo
          do j=jstart_res,jend_res+1
            delz_side(istart_res-1,j)=2.*delz_side(istart_res,j)        &
                                        -delz_side(istart_res+1,j)
            delz_side(iend_res+1,j)=2.*delz_side(iend_res,j)            &
                                      -delz_side(iend_res-1,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_delz_bc                         &
                                 ,delz_side(istart_res-1:iend_res+1,jstart_res:jend_res+1) &
                                 ,start=(/istart_bc-1,jstart_bc,k/)                    &
                                 ,count=(/len_x+2,len_y+1,1/)))
!
        elseif(side=='east ')then
          do j=jstart_res,jend_res
            delz_side(istart_res-1,j)=2.*delz_side(istart_res,j)        &
                                        -delz_side(istart_res+1,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_delz_bc                         &
                                 ,delz_side(istart_res-1:iend_res,jstart_res:jend_res) &
                                 ,start=(/istart_bc-1,jstart_bc,k/)                    &
                                 ,count=(/len_x+1,len_y,1/)))
!
        elseif(side=='west ')then
          do j=jstart_res,jend_res
            delz_side(iend_res+1,j)=2.*delz_side(iend_res,j)            &
                                        -delz_side(iend_res-1,j)
          enddo
!
          call check(nf90_put_var(ncid_bc_new,var_id_delz_bc                         &
                                 ,delz_side(istart_res:iend_res+1,jstart_res:jend_res) &
                                 ,start=(/istart_bc,jstart_bc,k/)                    &
                                 ,count=(/len_x+1,len_y,1/)))
!
        endif
!
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      deallocate(sphum_side)
      deallocate(T_side)
      deallocate(delp_side)
      deallocate(delz_side)
      deallocate(zh_bot)
      deallocate(zh_top)
      deallocate(pint)
!
!-----------------------------------------------------------------------
!
      end subroutine write_delz_to_bc_file
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      subroutine dgrid_to_cgrid(var)
!
!-----------------------------------------------------------------------
!***  The GSI updates the D-grid winds but the BC file also needs the
!***  C-grid winds.   Compute the C-grid winds in the boundary rows
!***  by interpolating from the D-grid winds and write them into the 
!***  BC file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      character(len=1),intent(in) :: var
!
!--------------------
!*** Local variables
!--------------------
!
      integer :: i,ix,j,jx,var_id
      integer :: istart_bc,jstart_bc,len_x,len_y
!
      real,dimension(:,:),allocatable :: uc_east,uc_north               &
                                        ,uc_south,uc_west               &
                                        ,vc_east,vc_north               &
                                        ,vc_south,vc_west
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The following 'start' and 'end' values are locations in the
!***  C-grid boundary rows relative to 1,1 since the values will
!***  use put_var to write them into the new BC netcdf file.  The
!***  loop indices start at 2,2 because the BC file arrays have 
!***  4 rows while the field_combined data from the enlarged restart 
!***  file has 3 boundary rows thus the new BC file rows are only 
!***  updated beginning in their 2nd row.  
!
!***  For C-grid components that lie on the outer edge of the 3rd
!***  boundary row simply extrapolate from the two adjacent C-grid
!***  components that have already been computed.
!-----------------------------------------------------------------------
!
      u_or_v: if(var=='u')then
!
!-----------------------------------------------------------------------
!
!--------------
!***  North uc
!--------------
!
        istart_bc=2
        iend_bc=dimsize_combined(2)+1
        jstart_bc=2
!xxx    jend_bc=halo_integrate+1
        jend_bc=halo_integrate+nrows_blend+1
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(uc_north))then
          allocate(uc_north(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
          do i=istart_bc+1,iend_bc-1
            uc_north(i,j)=0.25*(field_combined(i-2,j-1)                 &
                               +field_combined(i-2,j)                   &
                               +field_combined(i-1,j-1)                 &
                               +field_combined(i-1,j))
          enddo
          uc_north(istart_bc,j)=2.*uc_north(istart_bc+1,j)-uc_north(istart_bc+2,j)
          uc_north(iend_bc,j)=2.*uc_north(iend_bc-1,j)-uc_north(iend_bc-2,j)
        enddo
!
        call check(nf90_inq_varid(ncid_bc_new,'u_w_bottom',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,uc_north(:,:)                           &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!--------------
!***  South uc
!--------------
!
        istart_bc=2
        iend_bc=dimsize_combined(2)+1
        jstart_bc=1
!xxx    jend_bc=halo_integrate
        jend_bc=halo_integrate+nrows_blend
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(uc_south))then
          allocate(uc_south(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
!xxx      jx=j+dimsize_res(4)+halo_integrate
          jx=j+dimsize_res(4)+halo_integrate-nrows_blend
          do i=istart_bc+1,iend_bc-1
            uc_south(i,j)=0.25*(field_combined(i-2,jx)                  &
                               +field_combined(i-2,jx+1)                &
                               +field_combined(i-1,jx)                  &
                               +field_combined(i-1,jx+1))
          enddo
          uc_south(istart_bc,j)=2.*uc_south(istart_bc+1,j)-uc_south(istart_bc+2,j)
          uc_south(iend_bc,j)=2.*uc_south(iend_bc-1,j)-uc_south(iend_bc-2,j)
        enddo
!
!xxx    call check(nf90_inq_varid(ncid_bc,'u_w_top',var_id))
        call check(nf90_inq_varid(ncid_bc_new,'u_w_top',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,uc_south(:,:)                           &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!-------------
!***  East uc
!-------------
!
        istart_bc=2
!xxx    iend_bc=halo_integrate+2
        iend_bc=halo_integrate+nrows_blend+2
        jstart_bc=1
        jend_bc=dimsize_res(4)
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(uc_east))then
          allocate(uc_east(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
          jx=j+halo_integrate
          do i=istart_bc+1,iend_bc
            uc_east(i,j)=0.25*(field_combined(i-2,jx)                   &
                              +field_combined(i-2,jx+1)                 &
                              +field_combined(i-1,jx)                   &
                              +field_combined(i-1,jx+1))
          enddo
          uc_east(istart_bc,j)=2.*uc_east(istart_bc+1,j)-uc_east(istart_bc+2,j)
        enddo
!
!xxx    call check(nf90_inq_varid(ncid_bc,'u_w_left',var_id))
        call check(nf90_inq_varid(ncid_bc_new,'u_w_left',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,uc_east(:,:)                            &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!-------------
!***  West uc
!-------------
!
        istart_bc=1
!xxx    iend_bc=halo_integrate+1
        iend_bc=halo_integrate+nrows_blend+1
        jstart_bc=1
        jend_bc=dimsize_res(4)
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(uc_west))then
          allocate(uc_west(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
          jx=j+halo_integrate
          do i=istart_bc,iend_bc-1
!xxx        ix=i+dimsize_res(1)+halo_integrate
            ix=i+dimsize_res(1)+halo_integrate-nrows_blend
            uc_west(i,j)=0.25*(field_combined(ix-1,jx)                  &
                              +field_combined(ix-1,jx+1)                &
                              +field_combined(ix,jx)                    &
                              +field_combined(ix,jx+1))
          enddo
          uc_west(iend_bc,j)=2.*uc_west(iend_bc-1,j)-uc_west(iend_bc-2,j)
        enddo
!
!xxx    call check(nf90_inq_varid(ncid_bc,'u_w_right',var_id))
        call check(nf90_inq_varid(ncid_bc_new,'u_w_right',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,uc_west(:,:)                            &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!-----------------------------------------------------------------------
!
      elseif(var=='v')then
!
!-----------------------------------------------------------------------
!
!--------------
!***  North vc
!--------------
!
        istart_bc=2
        iend_bc=dimsize_combined(1)+1
        jstart_bc=2
!xxx    jend_bc=halo_integrate+2
        jend_bc=halo_integrate+nrows_blend+2
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(vc_north))then
          allocate(vc_north(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc+1,jend_bc
          do i=istart_bc,iend_bc
            vc_north(i,j)=0.25*(field_combined(i-1,j-2)                 &
                               +field_combined(i-1,j-1)                 &
                               +field_combined(i,j-2)                   &
                               +field_combined(i,j-1))
          enddo
        enddo
        do i=istart_bc,iend_bc
          vc_north(i,jstart_bc)=2.*vc_north(i,jstart_bc+1)-vc_north(i,jstart_bc+2)
        enddo
!
        call check(nf90_inq_varid(ncid_bc_new,'v_s_bottom',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,vc_north(:,:)                           &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!--------------
!***  South vc
!--------------
!
        istart_bc=2
        iend_bc=dimsize_combined(1)+1
        jstart_bc=1
!xxx    jend_bc=halo_integrate+1
        jend_bc=halo_integrate+nrows_blend+1
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(vc_south))then
          allocate(vc_south(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc-1
!xxx      jx=j+dimsize_res(4)+halo_integrate
          jx=j+dimsize_res(4)+halo_integrate-nrows_blend
          do i=istart_bc,iend_bc
            vc_south(i,j)=0.25*(field_combined(i-1,jx-1)                &
                               +field_combined(i-1,jx)                  &
                               +field_combined(i,jx-1)                  &
                               +field_combined(i,jx))
          enddo
        enddo
        do i=istart_bc,iend_bc
          vc_south(i,jend_bc)=2.*vc_south(i,jend_bc-1)-vc_south(i,jend_bc-2)
        enddo
!
        call check(nf90_inq_varid(ncid_bc_new,'v_s_top',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,vc_south(:,:)                           &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!-------------
!***  East vc
!-------------
!
        istart_bc=2
!xxx    iend_bc=halo_integrate+1
        iend_bc=halo_integrate+nrows_blend+1
        jstart_bc=1
        jend_bc=dimsize_res(3)-2
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(vc_east))then
          allocate(vc_east(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
          jx=j+halo_integrate
          do i=istart_bc,iend_bc
            vc_east(i,j)=0.25*(field_combined(i-1,jx)                   &
                              +field_combined(i-1,jx+1)                 &
                              +field_combined(i,jx)                     &
                              +field_combined(i,jx+1))
          enddo
        enddo
!
        call check(nf90_inq_varid(ncid_bc_new,'v_s_left',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,vc_east(:,:)                            &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
!-------------
!***  West vc
!-------------
!
        istart_bc=1
!xxx    iend_bc=halo_integrate
        iend_bc=halo_integrate+nrows_blend
        jstart_bc=1
        jend_bc=dimsize_res(3)-2
        len_x=iend_bc-istart_bc+1
        len_y=jend_bc-jstart_bc+1
!
        if(.not.allocated(vc_west))then
          allocate(vc_west(istart_bc:iend_bc,jstart_bc:jend_bc))
        endif
!
        do j=jstart_bc,jend_bc
          jx=j+halo_integrate
          do i=istart_bc,iend_bc
!xxx        ix=i+dimsize_res(1)+halo_integrate
            ix=i+dimsize_res(1)+halo_integrate-nrows_blend
            vc_west(i,j)=0.25*(field_combined(ix,jx)                    &
                              +field_combined(ix,jx+1)                  &
                              +field_combined(ix+1,jx)                  &
                              +field_combined(ix+1,jx+1))
          enddo
        enddo
!
        call check(nf90_inq_varid(ncid_bc_new,'v_s_right',var_id))
!
        call check(nf90_put_var(ncid_bc_new,var_id                      &
                               ,vc_west(:,:)                            &
                               ,start=(/istart_bc,jstart_bc,k/)         &
                               ,count=(/len_x,len_y,1/)))
!
      endif u_or_v
!
!-----------------------------------------------------------------------
!
      end subroutine dgrid_to_cgrid
!
!-----------------------------------------------------------------------
!
      end program move_DA_update_data
!
!-----------------------------------------------------------------------
