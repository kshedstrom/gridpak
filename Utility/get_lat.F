#include "griddefs.h"
      subroutine get_lat
!
!=======================================================================
!                                                                    ===
!  This subroutine reads lat/lon information from grid NetCDF file.  ===
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
      logical gotlatp, gotlonp, gotlatr, gotlonr, gotlatu, gotlonu, &
     &        gotlatv, gotlonv
      integer grdlapid, grdlopid, grdlauid, grdlouid, grdlavid, &
     &        grdlovid, grdlarid, grdlorid, i
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
!  Scan variable list from input NetCDF and check for latitude and
!  longitude variables.
!
      do i=1,nvars
        if (TRIM(varnam(i)).eq.'lat_rho') then
          grdlarid=i
          gotlatr=.true.
        elseif (TRIM(varnam(i)).eq.'lon_rho') then
          grdlorid=i
          gotlonr=.true.
        elseif (TRIM(varnam(i)).eq.'lat_psi') then
          grdlapid=i
          gotlatp=.true.
        elseif (TRIM(varnam(i)).eq.'lon_psi') then
          grdlopid=i
          gotlonp=.true.
        elseif (TRIM(varnam(i)).eq.'lat_u') then
          grdlauid=i
          gotlatu=.true.
        elseif (TRIM(varnam(i)).eq.'lon_u') then
          grdlouid=i
          gotlonu=.true.
        elseif (TRIM(varnam(i)).eq.'lat_v') then
          grdlavid=i
          gotlatv=.true.
        elseif (TRIM(varnam(i)).eq.'lon_v') then
          grdlovid=i
          gotlonv=.true.
        endif
      enddo
!
!  Terminate execution if essential grid variables are not found.
!
      if (.not.gotlatr) then
        write(stdout,10) 'lat_rho', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlonr) then
        write(stdout,10) 'lon_rho', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlatp) then
        write(stdout,10) 'lat_psi', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlonp) then
        write(stdout,10) 'lon_psi', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlatu) then
        write(stdout,10) 'lat_u', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlonu) then
        write(stdout,10) 'lon_u', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlatv) then
        write(stdout,10) 'lat_v', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
      if (.not.gotlonv) then
        write(stdout,10) 'lon_v', TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
!
!  Open grid NetCDF file for reading.
!
      rcode = nf90_open(TRIM(gridfile),nf90_nowrite,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,20) TRIM(gridfile)
        call crash ('GET_LAT',1)
      endif
!
!  Read in lat,lon coordinates at RHO-points.
!
      if (gotlatr) then
        rcode = nf90_get_var(ncgridid,grdlarid,lat_rho)
        if (rcode.ne.0) then
          write(stdout,30) 'lat_rho', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
      if (gotlonr) then
        rcode = nf90_get_var(ncgridid,grdlorid,lon_rho)
        if (rcode.ne.0) then
          write(stdout,30) 'lon_rho', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
!
!  Read in lat,lon coordinates at PSI-points.
!
      if (gotlatp) then
        rcode = nf90_get_var(ncgridid,grdlapid,lat_psi)
        if (rcode.ne.0) then
          write(stdout,30) 'lat_psi', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
      if (gotlonp) then
        rcode = nf90_get_var(ncgridid,grdlopid,lon_psi)
        if (rcode.ne.0) then
          write(stdout,30) 'lon_psi', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
!
!  Read in lat,lon coordinates at U-points.
!
      if (gotlatu) then
        rcode = nf90_get_var(ncgridid,grdlauid,lat_u)
        if (rcode.ne.0) then
          write(stdout,30) 'lat_u', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
      if (gotlonu) then
        rcode = nf90_get_var(ncgridid,grdlouid,lon_u)
        if (rcode.ne.0) then
          write(stdout,30) 'lon_u', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
!
!  Read in lat,lon coordinates at V-points.
!
      if (gotlatv) then
        rcode = nf90_get_var(ncgridid,grdlavid,lat_v)
        if (rcode.ne.0) then
          write(stdout,30) 'lat_v', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
      if (gotlonv) then
        rcode = nf90_get_var(ncgridid,grdlovid,lon_v)
        if (rcode.ne.0) then
          write(stdout,30) 'lon_v', TRIM(gridfile)
          call crash ('GET_LAT',1)
        endif
      endif
!
!  Close file
!
      rcode = nf90_close(ncgridid)
!
  10  format(/' GET_LAT - unable to find grid variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  20  format(/' GET_LAT - unable to open grid NetCDF file: ',a)
  30  format(/' GET_LAT - error while reading variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
      return
      end
