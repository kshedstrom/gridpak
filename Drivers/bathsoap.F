      program bathsoap

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!  This program reads and writes to the hraw variable in
!  the grid netCDF file.
!
!  It uses a Shapiro filter to smooth the bathymetry
! *******************************************************************

#include "griddefs.h"
#include "bathy.h"
#include "proj.h"
#include "ncgrid.h"
      BIGREAL         wrk1(L2d), wrk2(L2d)

      logical         colour, grover
      integer         i, j, k, ii, jj, imap, npass
      BIGREAL         vmin, vmax, rmax, rv
      real            U1, U2, V1, V2
      BIGREAL         udeg, uscale, DTOR, RTOD, &
     &                REarth, enlarge
      parameter (     REarth=6.3708e6  )
      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

#include "gridid.h"

      call get_xy
      call get_mn
      call get_h
      call get_rmask

!  Plot after each pass of Shapiro filter

      xmin = vmin(xp,L*M)
      ymin = vmin(yp,L*M)
      xmax = vmax(xp,L*M)
      ymax = vmax(yp,L*M)

      el = ymax-ymin
      xl = xmax-xmin

!  true for colour fill plots (false gives contours)
      colour = .true.

!  true to draw grid over colour plots
      grover = .false.

#if PLOTS
      call start_plot
#if DRAW_COASTS
!  put things on EZMAP space:
      call mapsti('GR',JGRD)
      call mapstc('OU','PS')
      call maproj(JPRJ,PLAT,PLONG,ROTA)
      call mapset(JLTS,P1,P2,P3,P4)
      call mapint

!  rescale xp, yp, xr, yr  (from Roberta Young)
      call maptrn((PLAT+.5),PLONG,U2,V2)
      call maptrn((PLAT-.5),PLONG,U1,V1)
      udeg = sqrt((U2-U1)*(U2-U1) + (V2-V1)*(V2-V1))
      uscale=DTOR*REarth/udeg
      write(6,*)' udeg   =',udeg
      write(6,*)' uscale =',uscale

      do j=1,M
        do i=1,L
          xp(i,j)=(xp(i,j)-XOFF)/uscale
          yp(i,j)=(yp(i,j)-YOFF)/uscale
        enddo
      enddo

      do j=0,M
        do i=0,L
          xr(i,j)=(xr(i,j)-XOFF)/uscale
          yr(i,j)=(yr(i,j)-YOFF)/uscale
        enddo
      enddo

!  find minimum x and y locations:
      xmin = vmin(xp,L*M)
      ymin = vmin(yp,L*M)
      xmax = vmax(xp,L*M)
      ymax = vmax(yp,L*M)
      el = ymax-ymin
      xl = xmax-xmin

!  make them larger for the plots
      enlarge = 0.03
      xmin = xmin - xl*enlarge
      ymin = ymin - el*enlarge
      xmax = xmax + xl*enlarge
      ymax = ymax + el*enlarge
#endif  /* DRAW_COASTS */

!  set foreground to black, background to white
      if (colour) then
        call gscr(1,0,1.,1.,1.)
        call gscr(1,1,0.,0.,0.)
      end if

      imap = 3
      call cpseti('MAP',imap)
      call getxxyy
      call pcseti('QUALITY',1)
#endif  /* PLOTS */

      print *
      print *, 'How many passes of the shapiro filter?'
      read *, npass

#if PLOTS
!     call ploth(gridid(1:40),colour,grover)
#endif  /* PLOTS */

      do k=1,npass
        call shapiro(h,Lp,Mp,4,wrk1,wrk2,mask_rho)
        rmax = 0.
        ii = 1
        jj = 1
        do j=1,M-1
          do i=1,L-1
            if (rv(i,j) .gt. rmax) then
              rmax = rv(i,j)
              ii = i
              jj = j
            endif
          enddo
        enddo
        print *,k,' rmax = ',rmax,' at (i,j): ',ii,jj
      enddo

!  Write the final filtered bathymetry to netCDF file
      call wrt_h

#if PLOTS
!  Plot the final bathymetry
      call ploth(gridid(1:40),colour,grover)
      call end_plot
#endif  /* PLOTS */

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
  100 call crash('bathsoap: read or write error', 0)
      end

! ********************************************************************

      subroutine crash(icrash,ierr)
      character*(*) icrash
      integer       ierr

      print *,icrash
      if (ierr .gt. 1) print *,'ierr = ',ierr
#if NO_EXIT
      stop
#else
      call exit(1)
#endif  /* NO_EXIT */
      return
      end

! ****************************************************************

      subroutine shapiro(u,Lp,Mp,N,tmp,tmp2,mask)
      integer    Lp, Mp, N
      BIGREAL    u(0:Lp-1,0:Mp-1)
      BIGREAL    mask(0:Lp-1,0:Mp-1)
      BIGREAL    tmp(0:Lp-1,0:Mp-1), tmp2(0:Lp-1,0:Mp-1)

! ****************************************************************

!  Version 1.0 April 27,1988   by Kate Hedstrom

! ****************************************************************

!  The user must supply two work arrays, tmp and tmp2, dimensioned
!  the same as the u array.

!  This subroutine will apply a Shapiro filter of order n (defined
!  as twice the order in Shapiro (1970), with N even) to an
!  array u.  The order of the filter is reduced at the boundaries.
!  The size of the u array is 0-L by 0-M.

!  This filter can be used for any rectangular grid.  If the grid
!  spacing is non-uniform then you have to scale by the grid spacing
!  before filtering and scale back after filtering in order to
!  conserve the quantity.  If conservation is not important then the
!  filter can be used as is even for non-uniform spacing.

      integer    d, L, M, Lm, Mm, i, j, k

      L = Lp-1
      M = Mp-1
      Lm = Lp-2
      Mm = Mp-2

      if (mod(N,2).ne.0) then
        print *,'N must be even in the shapiro filter'
#if NO_EXIT
        stop
#else
        call exit(1)
#endif  /* NO_EXIT */
      end if

!  Do the first y pass to initialize the temporary array

      do 100 i=0,L
      do 100 j=1,M-1
        tmp(i,j) = 0.25 * (u(i,j-1)*mask(i,j-1) + u(i,j+1)*mask(i,j+1)  &
     &            - 2*u(i,j)*mask(i,j)) * &
     &         mask(i,j-1) * mask(i,j+1) * mask(i,j)
  100 continue

!  Other passes in the y direction.

      do 120 k=4,N,2
        d = k/2
        do 110 j=d,M-d
        do 110 i=0,L
          tmp2(i,j) = - 0.25 * (tmp(i,j-1)*mask(i,j-1) &
     &            + tmp(i,j+1)*mask(i,j+1) &
     &            - 2*tmp(i,j)*mask(i,j)) * &
     &         mask(i,j-1) * mask(i,j+1) * mask(i,j)
  110   continue
      do 120 j=d,M-d
      do 120 i=0,L
        tmp(i,j) = tmp2(i,j)
  120 continue

!  Add the changes to u

      do 130 j=1,M-1
      do 130 i=0,L
        u(i,j) = u(i,j) + tmp(i,j)
  130 continue

!  Initialize tmp to filter in the x direction.

      do 140 j=0,M
      do 140 i=1,L-1
        tmp(i,j) = 0.25 * (u(i-1,j)*mask(i-1,j) + u(i+1,j)*mask(i+1,j) &
     &            - 2*u(i,j)*mask(i,j)) * &
     &         mask(i-1,j) * mask(i+1,j) * mask(i,j)
  140 continue

!  Other x passes

      do 160 k=4,N,2
        d = k/2
      do 160 j=0,M
        do 150 i=d,L-d
          tmp2(i,j) = - 0.25 * (tmp(i-1,j)*mask(i-1,j) &
     &                       + tmp(i+1,j)*mask(i+1,j) &
     &            - 2*tmp(i,j)*mask(i,j)) * &
     &         mask(i-1,j) * mask(i+1,j) * mask(i,j)
  150   continue
      do 160 i=d,L-d
        tmp(i,j) = tmp2(i,j)
  160 continue

!  Add changes to u

      do 170 j=0,M
      do 170 i=1,L-1
        u(i,j) = u(i,j) + tmp(i,j)
  170 continue

!  Make edge values like interior
      do 180 j=1,Mm
        u(0,j) = u(1,j)
        u(L,j) = u(Lm,j)
  180 continue
      do 190 i=0,L
        u(i,0) = u(i,1)
        u(i,M) = u(i,Mm)
  190 continue

      return
      end

      BIGREAL function rv(i,j)
      integer   i, j

!     calculate the SPEM r-value
#include "bathy.h"
      BIGREAL   dhdxx, dhdyy

      dhdxx  = abs((h(i,j) - h(i-1,j)) / (h(i,j) + h(i-1,j))) &
     &         * mask_rho(i,j) * mask_rho(i-1,j)
      dhdyy  = abs((h(i,j) - h(i,j-1)) / (h(i,j) + h(i,j-1))) &
     &         * mask_rho(i,j) * mask_rho(i,j-1)
      rv     = max(dhdxx,dhdyy)
      return
      end
