#include "griddefs.h"
      subroutine get_mn
!
!=======================================================================
!                                                                    ===
!  This subroutine reads grid metric information from grid NetCDF    ===
!  file.                                                             ===
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
      logical gotpm, gotpn
      integer grdpmid, grdpnid, i
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
      call opencdf(TRIM(gridfile))
!
!  Scan variable list from input NetCDF and check for pm and pn
!  variables.
!
      do i=1,nvars
        if (TRIM(varnam(i)).eq.'pm') then
          grdpmid=i
          gotpm=.true.
        endif
        if (TRIM(varnam(i)).eq.'pn') then
          grdpnid=i
          gotpn=.true.
        endif
      enddo
!
!  Terminate execution if essential grid variables are not found.
!
      if (.not.gotpm) then
        write(stdout,10) 'pm', TRIM(gridfile)
        call crash ('GET_MN',1)
      endif
      if (.not.gotpn) then
        write(stdout,10) 'pn', TRIM(gridfile)
        call crash ('GET_MN',1)
      endif
!
!  Open grid NetCDF file for reading.
!
      rcode=nf90_open(TRIM(gridfile),nf90_nowrite,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,20) TRIM(gridfile)
        call crash ('GET_MN',1)
      endif
!
!  Read in grid metrics.
!
      if (gotpm) then
        rcode = nf90_get_var(ncgridid,grdpmid,pm)
        if (rcode.ne.0) then
          write(stdout,30) 'pm', TRIM(gridfile)
          call crash ('GET_MN',1)
        endif
      endif
      if (gotpn) then
        rcode = nf90_get_var(ncgridid,grdpnid,pn)
        if (rcode.ne.0) then
          write(stdout,30) 'pn', TRIM(gridfile)
          call crash ('GET_MN',1)
        endif
      endif
!
!  Close file
!
      rcode = nf90_close(ncgridid)
!
  10  format(/' GET_MN - unable to find grid variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  20  format(/' GET_MN - unable to open grid NetCDF file: ',a)
  30  format(/' GET_MN - error while reading variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
      return
      end
