#include "griddefs.h"
      subroutine get_xy
!
!=======================================================================
!                                                                    ===
!  This subroutine reads x,y grid information from grid NetCDF file. ===
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
      logical gotxp, gotyp, gotxu, gotyu, gotxv, gotyv, &
     &        gotxr, gotyr
      integer grdxpid, grdypid, grdxuid, grdyuid, grdxvid, grdyvid, &
     &        grdxrid, grdyrid, i
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Inquire about the contents of SCRUM grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      call opencdf (TRIM(gridfile))
!
!  Scan variable list from input NetCDF and check for x,y variables.
!
      do i=1,nvars
        if (TRIM(varnam(i)).eq.'x_rho') then
          grdxrid=i
          gotxr=.true.
        elseif (TRIM(varnam(i)).eq.'y_rho') then
          grdyrid=i
          gotyr=.true.
        elseif (TRIM(varnam(i)).eq.'x_psi') then
          grdxpid=i
          gotxp=.true.
        elseif (TRIM(varnam(i)).eq.'y_psi') then
          grdypid=i
          gotyp=.true.
        elseif (TRIM(varnam(i)).eq.'x_u') then
          grdxuid=i
          gotxu=.true.
        elseif (TRIM(varnam(i)).eq.'y_u') then
          grdyuid=i
          gotyu=.true.
        elseif (TRIM(varnam(i)).eq.'x_v') then
          grdxvid=i
          gotxv=.true.
        elseif (TRIM(varnam(i)).eq.'y_v') then
          grdyvid=i
          gotyv=.true.
        endif
      enddo
!
!  Terminate execution if essential grid variables are not found.
!
      if (.not.gotxr) then
        write(stdout,10) 'xr', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotyr) then
        write(stdout,10) 'yr', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotxp) then
        write(stdout,10) 'xp', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotyp) then
        write(stdout,10) 'yp', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotxu) then
        write(stdout,10) 'xu', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotyu) then
        write(stdout,10) 'yu', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotxv) then
        write(stdout,10) 'xv', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
      if (.not.gotyv) then
        write(stdout,10) 'yv', TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
!
!  Open grid NetCDF file for reading.
!
      rcode = nf90_open(TRIM(gridfile),nf90_nowrite,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,20) TRIM(gridfile)
        call crash ('GET_XY',1)
      endif
!
!  Read in (x,y) coordinates at RHO-points.
!
      if (gotxr) then
        rcode = nf90_get_var(ncgridid,grdxrid,xr)
        if (rcode.ne.0) then
          write(stdout,30) 'x_rho', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
      if (gotyr) then
        rcode = nf90_get_var(ncgridid,grdyrid,yr)
        if (rcode.ne.0) then
          write(stdout,30) 'y_rho', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
!
!  Read in (x,y) coordinates at PSI-points.
!
      if (gotxp) then
        rcode = nf90_get_var(ncgridid,grdxpid,xp)
        if (rcode.ne.0) then
          write(stdout,30) 'x_psi', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
      if (gotyp) then
        rcode = nf90_get_var(ncgridid,grdypid,yp)
        if (rcode.ne.0) then
          write(stdout,30) 'y_psi', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
!
!  Read in (x,y) coordinates at U-points.
!
      if (gotxu) then
        rcode = nf90_get_var(ncgridid,grdxuid,xu)
        if (rcode.ne.0) then
          write(stdout,30) 'x_u', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
      if (gotyu) then
        rcode = nf90_get_var(ncgridid,grdyuid,yu)
        if (rcode.ne.0) then
          write(stdout,30) 'y_u', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
!
!  Read in (x,y) coordinates at V-points.
!
      if (gotxv) then
        rcode = nf90_get_var(ncgridid,grdxvid,xv)
        if (rcode.ne.0) then
          write(stdout,30) 'x_v', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
      if (gotyv) then
        rcode = nf90_get_var(ncgridid,grdyvid,yv)
        if (rcode.ne.0) then
          write(stdout,30) 'y_v', TRIM(gridfile)
          call crash ('GET_XY',1)
        endif
      endif
!
!  Close file
!
      rcode = nf90_close(ncgridid)
!
  10  format(/' GET_XY - unable to find grid variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  20  format(/' GET_XY - unable to open grid NetCDF file: ',a)
  30  format(/' GET_XY - error while reading variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
      return
      end
