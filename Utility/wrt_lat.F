#include "griddefs.h"

      subroutine wrt_lat
!
!=======================================================================
!                                                                    ===
!  This routine writes latitude/longitude fields into grid           ===
!  NetCDF file.                                                      ===
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
      real    lat(2)
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
          call crash ('WRT_LAT',1)
        endif
!
!  Write out map projection parameters.
!
      rcode=nf90_inq_varid(ncgridid,'JPRJ',varid)
      rcode = nf90_put_var(ncgridid,varid,JPRJ)
      lat(1) = PLAT
      if (JPRJ .eq. 'LC') then
        lat(2) = ROTA
      else
        lat(2) = 0.
      endif
      rcode=nf90_inq_varid(ncgridid,'PLAT',varid)
      rcode = nf90_put_var(ncgridid,varid,lat)
      rcode=nf90_inq_varid(ncgridid,'PLONG',varid)
      rcode = nf90_put_var(ncgridid,varid,PLONG)
      if (JPRJ .eq. 'LC') then
        lat(1) = 0.
      else
        lat(1) = ROTA
      endif
      rcode=nf90_inq_varid(ncgridid,'ROTA',varid)
      rcode = nf90_put_var(ncgridid,varid,lat)

      rcode=nf90_inq_varid(ncgridid,'JLTS',varid)
      rcode = nf90_put_var(ncgridid,varid,JLTS)
      rcode=nf90_inq_varid(ncgridid,'P1',varid)
      rcode = nf90_put_var(ncgridid,varid,P1)
      rcode=nf90_inq_varid(ncgridid,'P2',varid)
      rcode = nf90_put_var(ncgridid,varid,P2)
      rcode=nf90_inq_varid(ncgridid,'P3',varid)
      rcode = nf90_put_var(ncgridid,varid,P3)
      rcode=nf90_inq_varid(ncgridid,'P4',varid)
      rcode = nf90_put_var(ncgridid,varid,P4)
      rcode=nf90_inq_varid(ncgridid,'XOFF',varid)
      rcode = nf90_put_var(ncgridid,varid,XOFF)
      rcode=nf90_inq_varid(ncgridid,'YOFF',varid)
      rcode = nf90_put_var(ncgridid,varid,YOFF)
!
!  Write out lat/lon at RHO-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_rho', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_rho', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
!
!  Write out lat/lon at PSI-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_psi', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_psi', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
!
!  Write out lat/lon at U-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_u',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_u)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_u', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_u',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_u)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_u', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
!
!  Write out lat/lon at V-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_v',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_v)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_v', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_v',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_v)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_v', TRIM(gridfile)
        call crash ('WRT_LAT',1)
      endif
!
!  Synchronize restart NetCDF file to disk to allow other processes to
!  access data immediately after it is written.
!
      rcode = nf90_sync(ncgridid)
      if (rcode.ne.0) then
        write(stdout,20)
        call crash ('WRT_LAT',1)
      endif
!
  10  format(/,' WRT_LAT - error while writing variable: ',a,/,11x, &
     &       'into grid NetCDF file for time record: ',i4)
  20  format(/,' WRT_LAT - unable to synchronize grid NetCDF to ', &
     &       'disk.')
      return
      end
