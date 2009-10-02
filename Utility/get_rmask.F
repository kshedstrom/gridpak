#include "griddefs.h"
      subroutine get_rmask
!
!=======================================================================
!                                                                    ===
!  This subroutine reads mask information from grid NetCDF file.     ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global variables.
!-----------------------------------------------------------------------
!
      use netcdf

# include "bathy.h"
# include "ncgrid.h"
!
!-----------------------------------------------------------------------
!  Define local variables.
!-----------------------------------------------------------------------
!
      logical gotrmask
      integer grdrmaskid, i
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      call opencdf (TRIM(gridfile))
!
!  Scan variable list from input NetCDF and check for rho mask.
!
      do i=1,nvars
        if (TRIM(varnam(i)).eq.'mask_rho') then
          grdrmaskid=i
          gotrmask=.true.
        endif
      enddo
!
!  Terminate execution if essential grid variables are not found.
!
      if (.not.gotrmask) then
        write(stdout,10) 'mask_rho', TRIM(gridfile)
        call crash ('GET_RMASK',1)
      endif
!
!  Open grid NetCDF file for reading.
!
      rcode=nf90_open(TRIM(gridfile),nf90_nowrite,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,20) TRIM(gridfile)
        call crash ('GET_RMASK',1)
      endif
!
!  Read in rho mask.
!
      if (gotrmask) then
        rcode = nf90_get_var(ncgridid,grdrmaskid,mask_rho)
        if (rcode.ne.0) then
          write(stdout,30) 'mask_rho', TRIM(gridfile)
          call crash ('GET_RMASK',1)
        endif
      endif
!
!  Close file
!
      rcode = nf90_close(ncgridid)
!
  10  format(/' GET_RMASK - unable to find grid variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  20  format(/' GET_RMASK - unable to open grid NetCDF file: ',a)
  30  format(/' GET_RMASK - error while reading variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
      return
      end
