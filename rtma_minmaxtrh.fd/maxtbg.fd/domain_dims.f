       subroutine domain_dims(cgrid,nx,ny,ds)

       implicit none
       
       character(30) cgrid
       include 'param.incl'

       integer(4) nx,ny
       real(4) ds

       if (trim(cgrid) == 'conus') then
           nx=nx_conus   
           ny=ny_conus
           ds=da_conus

        elseif (trim(cgrid) == 'alaska') then 
           nx=nx_alaska  
           ny=ny_alaska
           ds=da_alaska

        elseif (trim(cgrid) == 'hawaii') then 
           nx=nx_hawaii  
           ny=ny_hawaii
           ds=da_hawaii

        elseif (trim(cgrid) == 'prico') then 
           nx=nx_prico
           ny=ny_prico
           ds=da_prico

        elseif (trim(cgrid) == 'guam') then 
           nx=nx_guam    
           ny=ny_guam
           ds=da_guam

        elseif (trim(cgrid) == 'cohres') then
           nx=nx_cohres
           ny=ny_cohres
           ds=da_cohres

        elseif (trim(cgrid) == 'akhres') then
           nx=nx_akhres
           ny=ny_akhres
           ds=da_akhres

        elseif (trim(cgrid) == 'hrrr') then
           nx=nx_hrrr
           ny=ny_hrrr
           ds=da_hrrr

        elseif (trim(cgrid) == 'juneau') then
           nx=nx_juneau
           ny=ny_juneau
           ds=da_juneau

        elseif (trim(cgrid) == 'cohresext') then
           nx=nx_cohresext
           ny=ny_cohresext
           ds=da_cohresext

        elseif (trim(cgrid) == 'cohreswexp') then
           nx=nx_cohreswexp
           ny=ny_cohreswexp
           ds=da_cohreswexp

        else
           print*,'in domain_dims: unknown grid ',cgrid,'...aborting'
           call abort
        endif

        return
        end
