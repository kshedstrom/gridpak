#include "griddefs.h"

      subroutine wrt_all
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
      integer level
      integer start(3), count(3)
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
        call crash ('WRT_ALL',1)
      endif
!
!  Write out Coriolis parameter.
!
      rcode = nf90_inq_varid(ncgridid,'f',varid)
      rcode = nf90_put_var(ncgridid,varid,f)
      if (rcode.ne.0) then
        write(stdout,10) 'f', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out bathymetry.
!
      rcode = nf90_inq_varid(ncgridid,'h',varid)
      rcode = nf90_put_var(ncgridid,varid,h)
      if (rcode.ne.0) then
        write(stdout,10) 'h', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Find out what slice to write into.
!
      print *, 'There are ', bathsize, ' bathymetries.'
      print *, 'Which level would you like to write to?'
      read (5,*) level
      if (level .lt. 1) then
        call crash('Illegal level', level)
      endif
      if (level .gt. bathsize) then
        level = bathsize+1
        bathsize = level
      endif
!
!  Write out next bathymetry into hraw.
!
      start(1)=1
      count(1)=Lp
      start(2)=1
      count(2)=Mp
      start(3)=level
      count(3)=1
      rcode=nf90_inq_varid(ncgridid,'hraw',varid)
      rcode = nf90_put_var(ncgridid,varid,h,start,count)
      if (rcode.ne.0) then
        write(stdout,10) 'hraw', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out grid metrics.
!
      rcode=nf90_inq_varid(ncgridid,'spherical',varid)
      if (spherical == 1) then
        rcode = nf90_put_var(ncgridid,varid,1)
      else
        rcode = nf90_put_var(ncgridid,varid,0)
      endif
      if (rcode.ne.0) then
        write(stdout,20) 'spherical', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'pm',varid)
      rcode = nf90_put_var(ncgridid,varid,pm)
      if (rcode.ne.0) then
        write(stdout,20) 'pm', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'pn',varid)
      rcode = nf90_put_var(ncgridid,varid,pn)
      if (rcode.ne.0) then
        write(stdout,20) 'pn', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'dndx',varid)
      rcode = nf90_put_var(ncgridid,varid,dndx)
      if (rcode.ne.0) then
        write(stdout,20) 'dndx', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'dmde',varid)
      rcode = nf90_put_var(ncgridid,varid,dmde)
      if (rcode.ne.0) then
        write(stdout,20) 'dmde', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  x,y grids.
!
      rcode=nf90_inq_varid(ncgridid,'x_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,xr)
      if (rcode.ne.0) then
        write(stdout,20) 'x_rho', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,yr)
      if (rcode.ne.0) then
        write(stdout,20) 'y_rho', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,xp)
      if (rcode.ne.0) then
        write(stdout,20) 'x_psi', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,yp)
      if (rcode.ne.0) then
        write(stdout,20) 'y_psi', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_u',varid)
      rcode = nf90_put_var(ncgridid,varid,xu)
      if (rcode.ne.0) then
        write(stdout,20) 'x_u', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_u',varid)
      rcode = nf90_put_var(ncgridid,varid,yu)
      if (rcode.ne.0) then
        write(stdout,20) 'y_u', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_v',varid)
      rcode = nf90_put_var(ncgridid,varid,xv)
      if (rcode.ne.0) then
        write(stdout,20) 'x_v', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_v',varid)
      rcode = nf90_put_var(ncgridid,varid,yv)
      if (rcode.ne.0) then
        write(stdout,20) 'y_v', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out lat/lon at RHO-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_rho', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_rho', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out lat/lon at PSI-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_psi', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_psi', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out lat/lon at U-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_u',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_u)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_u', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_u',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_u)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_u', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out lat/lon at V-points.
!
      rcode=nf90_inq_varid(ncgridid,'lat_v',varid)
      rcode = nf90_put_var(ncgridid,varid,lat_v)
      if (rcode.ne.0) then
        write(stdout,10) 'lat_v', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'lon_v',varid)
      rcode = nf90_put_var(ncgridid,varid,lon_v)
      if (rcode.ne.0) then
        write(stdout,10) 'lon_v', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out mask at RHO-points.
!
      rcode = nf90_inq_varid(ncgridid,'mask_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_rho)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_rho', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out mask at U-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_u',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_u)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_u', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out mask at V-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_v',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_v)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_v', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out mask at PSI-points.
!
      rcode=nf90_inq_varid(ncgridid,'mask_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,mask_psi)
      if (rcode.ne.0) then
        write(stdout,10) 'mask_psi', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Write out grid angle.
!
      rcode = nf90_inq_varid(ncgridid,'angle',varid)
      rcode = nf90_put_var(ncgridid,varid,angle)
      if (rcode.ne.0) then
        write(stdout,10) 'angle', TRIM(gridfile)
        call crash ('WRT_ALL',1)
      endif
!
!  Synchronize restart NetCDF file to disk to allow other processes to
!  access data immediately after it is written.
!
      rcode = nf90_sync(ncgridid)
      if (rcode.ne.0) then
        write(stdout,20)
        call crash ('WRT_ALL',1)
      endif
!
  10  format(/,' WRT_ALL - error while writing variable: ',a,/,11x, &
     &       'into grid NetCDF file for time record: ',i4)
  20  format(/,' WRT_ALL - unable to synchronize grid NetCDF to ', &
     &       'disk.')
      return
      end
