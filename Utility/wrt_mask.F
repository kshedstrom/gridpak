#include "griddefs.h"

      subroutine wrt_mask
!
!=======================================================================
!                                                                    ===
!  This routine writes the masks into grid NetCDF file.              ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global data.
!-----------------------------------------------------------------------
!
      use netcdf

#include "bathy.h"
#include "ncgrid.h"
!
!-----------------------------------------------------------------------
!  Define local data.
!-----------------------------------------------------------------------
!
      BIGREAL bvar
!
!=======================================================================
!  Begin Executable code.
!=======================================================================
!
!  Open restart file for read/write.
!
        rcode=nf90_open(TRIM(gridfile),nf90_write,ncgridid)
        if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
          write(stdout,10) TRIM(gridfile)
          call crash ('WRT_MASK',1)
        endif
!
!  Write out bathymetry clipping depths.
!
      bvar = depthmin
      rcode = nf90_inq_varid(ncgridid,'depthmin',varid)
      rcode = nf90_put_var(ncgridid,varid,bvar)
      bvar = depthmax
      rcode = nf90_inq_varid(ncgridid,'depthmax',varid)
      rcode = nf90_put_var(ncgridid,varid,bvar)
!
!  Write out mask at RHO-points.
!
      rcode = nf90_inq_varid(ncgridid,'mask_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_rho', TRIM(gridfile)
        call crash ('WRT_MASK',1)
      endif
!
!  Write out mask at U-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_u',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_u)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_u', TRIM(gridfile)
        call crash ('WRT_MASK',1)
      endif
!
!  Write out mask at V-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_v',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_v)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_v', TRIM(gridfile)
        call crash ('WRT_MASK',1)
      endif
!
!  Write out mask at PSI-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_psi', TRIM(gridfile)
        call crash ('WRT_MASK',1)
      endif
!
!  Synchronize restart NetCDF file to disk to allow other processes to
!  access data immediately after it is written.
!
      rcode = nf90_sync(ncgridid)
      if (rcode.ne.0) then
        write(stdout,20)
        call crash ('WRT_MASK',1)
      endif
!
  10  format(/,' WRT_MASK - error while writing variable: ',a,/,11x, &
     &       'into grid NetCDF file for time record: ',i4)
  20  format(/,' WRT_MASK - unable to synchronize grid NetCDF to ', &
     &       'disk.')
      return
      end
