#include "griddefs.h"

      subroutine wrt_fhmn(f0,beta)
!
!=======================================================================
!                                                                    ===
!  This routine writes fields into the grid NetCDF file.             ===
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
#include "proj.h"
!
!-----------------------------------------------------------------------
!  Define local data.
!-----------------------------------------------------------------------
!
      BIGREAL f0, beta
!
!=======================================================================
!  Begin Executable code.
!=======================================================================
!
!  Open restart file for read/write.
!
      rcode = nf90_open(TRIM(gridfile),nf90_write,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,10) TRIM(gridfile)
        call crash ('WRT_FHMN',1)
      endif
!
!  Write out Coriolis parameter.
!
      if (spherical == 0) then
        rcode = nf90_inq_varid(ncgridid,'spherical',varid)
        rcode = nf90_put_var(ncgridid,varid,0)
        if (rcode.ne.0) then
          write(stdout,10) 'spherical', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'f0',varid)
        rcode = nf90_put_var(ncgridid,varid,f0)
        if (rcode.ne.0) then
          write(stdout,10) 'f0', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'dfdy',varid)
        rcode = nf90_put_var(ncgridid,varid,beta)
        if (rcode.ne.0) then
          write(stdout,10) 'dfdy', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
      endif
      rcode = nf90_inq_varid(ncgridid,'f',varid)
      rcode = nf90_put_var(ncgridid,varid,f)
      if (rcode.ne.0) then
        write(stdout,10) 'f', TRIM(gridfile)
        call crash ('WRT_FHMN',1)
      endif
!
!  Write out bathymetry.
!
      rcode = nf90_inq_varid(ncgridid,'h',varid)
      rcode = nf90_put_var(ncgridid,varid,h)
      if (rcode.ne.0) then
        write(stdout,10) 'h', TRIM(gridfile)
        call crash ('WRT_FHMN',1)
      endif
!
!  Write out grid metrics.
!
      rcode = nf90_inq_varid(ncgridid,'spherical',varid)
      if (spherical == 1) then
        rcode = nf90_put_var(ncgridid,varid,1)
        if (rcode.ne.0) then
          write(stdout,20) 'spherical', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'pm',varid)
        rcode = nf90_put_var(ncgridid,varid,pm)
        if (rcode.ne.0) then
          write(stdout,10) 'pm', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'pn',varid)
        rcode = nf90_put_var(ncgridid,varid,pn)
        if (rcode.ne.0) then
          write(stdout,10) 'pn', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'dndx',varid)
        rcode = nf90_put_var(ncgridid,varid,dndx)
        if (rcode.ne.0) then
          write(stdout,10) 'dndx', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
        rcode = nf90_inq_varid(ncgridid,'dmde',varid)
        rcode = nf90_put_var(ncgridid,varid,dmde)
        if (rcode.ne.0) then
          write(stdout,10) 'dmde', TRIM(gridfile)
          call crash ('WRT_FHMN',1)
        endif
      endif
!
!  Write out grid angle.
!
      rcode = nf90_inq_varid(ncgridid,'angle',varid)
      rcode = nf90_put_var(ncgridid,varid,angle)
      if (rcode.ne.0) then
        write(stdout,10) 'angle', TRIM(gridfile)
        call crash ('WRT_FHMN',1)
      endif
!
!  Synchronize restart NetCDF file to disk to allow other processes to
!  access data immediately after it is written.
!
      rcode = nf90_sync(ncgridid)
      if (rcode.ne.0) then
        write(stdout,20)
        call crash ('WRT_FHMN',1)
      endif
!
  10  format(/,' WRT_FHMN - error while writing variable: ',a,/,11x, &
     &       'into grid NetCDF file for time record: ',i4)
  20  format(/,' WRT_FHMN - unable to synchronize grid NetCDF to ', &
     &       'disk.')
      return
      end
