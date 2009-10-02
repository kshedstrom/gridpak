#include "griddefs.h"
      subroutine get_h
!
!=======================================================================
!                                                                    ===
!  This subroutine reads bathymetry information from grid NetCDF     ===
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
      logical gothraw
      integer grdhrid, i, level
      integer count(3), start(3)
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
      if (bathsize .lt. 1) then
        write(stdout,40) TRIM(gridfile)
        call crash ('GET_H',1)
      endif
!
!  Find out which level to read.
!
      print *, 'There are ', bathsize, ' bathymetries.'
      print *, 'Which level would you like to read?'
      read (5,*) level
      if (level .lt. 1 .or. level .gt. bathsize) then
        call crash('Illegal level', level)
      endif
!
!  Scan variable list from input NetCDF and check for raw bathymetry.
!
      do i=1,nvars
        if (TRIM(varnam(i)).eq.'hraw') then
          grdhrid=i
          gothraw=.true.
        endif
      enddo
!
!  Terminate execution if essential grid variables are not found.
!
      if (.not.gothraw) then
        write(stdout,10) 'hraw', TRIM(gridfile)
        call crash ('GET_H',1)
      endif
!
!  Open grid NetCDF file for reading.
!
      rcode=nf90_open(TRIM(gridfile),nf90_nowrite,ncgridid)
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,20) TRIM(gridfile)
        call crash ('GET_H',1)
      endif
!
!  Read in last bathymetry in hraw.
!
      start(1)=1
      count(1)=Lp
      start(2)=1
      count(2)=Mp
      start(3)=level
      count(3)=1
      if (gothraw) then
        rcode = nf90_get_var(ncgridid,grdhrid,h,start,count)
        if (rcode.ne.0) then
          write(stdout,30) 'hraw', TRIM(gridfile)
          call crash ('GET_H',1)
        endif
      endif
!
!  Close file
!
      rcode = nf90_close(ncgridid)
!
  10  format(/' GET_H - unable to find grid variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  20  format(/' GET_H - unable to open grid NetCDF file: ',a)
  30  format(/' GET_H - error while reading variable: ',a, &
     &       /12x,'in grid NetCDF file: ',a)
  40  format(/' GET_H - no bathymetries have been written yet: ',a)
      return
      end
