#include "griddefs.h"

      subroutine wrt_h
!
!=======================================================================
!                                                                    ===
!  This routine writes bathymetry into grid NetCDF file.             ===
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
      integer count(3), start(3)
      integer level
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
          call crash ('WRT_H',1)
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
        call crash ('WRT_H',1)
      endif
!
!  Synchronize restart NetCDF file to disk to allow other processes to
!  access data immediately after it is written.
!
      rcode = nf90_sync(ncgridid)
      if (rcode.ne.0) then
        write(stdout,20)
        call crash ('WRT_H',1)
      endif
!
  10  format(/,' WRT_H - error while writing variable: ',a,/,11x, &
     &       'into grid NetCDF file for time record: ',i4)
  20  format(/,' WRT_H - unable to synchronize grid NetCDF to ', &
     &       'disk.')
      return
      end
