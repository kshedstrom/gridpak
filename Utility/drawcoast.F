! *********************************************************************
!   Copyright (c) 1991, 1993 Rutgers University
!
! *********************************************************************

      subroutine drawcoast
#include "griddefs.h"
      real           lat, lon

!  MAX is maximum number of segments in data file
      integer        MAX
      parameter    ( MAX=100000    )
      character*80   datafile
      integer        icst, nn

      icst = 40
      call getenv('XCOASTDATA', datafile)
      open(icst,file=datafile,err=130)

!  I want a "while (scanf() != EOF)", but I have to fake it.
      do 100 nn=1,MAX
        read (icst,*,err=120,end=110) lat, lon
        call mapit(lat, lon, 0)
        do while (abs(lat) .lt. 89.99)
          call mapit(lat, lon, 2)
          read (icst,*,err=120,end=110) lat, lon
        enddo
  100 continue

      print *,'Did not reach end of file in drawcoast'
      print *,'Increase MAX to draw the rest of your file'
  110 continue
      call mapiq
      close(icst)
      return

  120 print *,'read error in drawcoast'
#if NO_EXIT
      stop
#else
      call exit(1)
#endif  /* NO_EXIT */
  130 print *, 'error opening file in drawcoast ', datafile
#if NO_EXIT
      stop
#else
      call exit(1)
#endif  /* NO_EXIT */
      end
