#include "griddefs.h"
      program bathsuds

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!  bathsuds (Selective User-Defined Smoothing) is a variation on
!  bathsoap (Smooth Over All Points) which allows the user to fiddle
!  endlessly with a set of weights applied to the Shapiro filter.
!  functions rvwgt{x,y} weight the correction term over the Shapiro
!  filter (whicyh are usually 1).  For example, if repeated sweeps of
!  the filter are being used to eliminate isolated steep regions of the
!  bathymetry, the rvwgt{x,y} weights can be set to zero over the
!  already smooth region to avoid needlessly reducing the domain to a
!  flat abyssal plain.
!
!  The present set of weights does this by computing Aike's r-value
!  and setting the weights to rvwmin where r-value is less than rv0,
!  and rvwmax where r-value is greater than rv0 (with a smooth tanh
!  transition between the two).  rvwmax > 1.0 amounts to a successive
!  over-relaxation of the fliter weights.  Alternatively, making rvwgt
!  simple functions of i,j will apply the filtering to predetermined
!  regions of the domain.
!
!  There is now a KEEP_SHALLOW option to uplift the shallow areas after
!  each pass of the filter.  There are tunable parameters to specify
!  the shallow depths which are being attempted to be kept - the
!  shallower the values you pick, the harder it is to get to a useful
!  r-value.  With KEEP_SHALLOW defined, you are likely to need many
!  more passes of the filter.  For instance, a domain in which you
!  achieved a reasonable r-value after 20 passes may now require
!  several hundred or a thousand passes.  Note that you can do
!  something similar to keep deep areas deep.  This can be useful
!  in regions surrounded by a coasts that would otherwise get too
!  shallow, such as in the Caribbean.
!
!  This program reads and writes to the hraw variable in the
!  grid netCDF file.
!
!  *******************************************************************

#include "bathy.h"
#include "ncgrid.h"
#include "proj.h"
      BIGREAL         wrk1(L2d), wrk2(L2d)
      logical         colour, grover
      integer         imap, npasses, i, j, k, ii, jj
      BIGREAL         vmin, vmax, rv, rmax
      real            U1, U2, V1, V2
      BIGREAL         udeg, uscale, DTOR, RTOD,                         &
     &                REarth, enlarge
      parameter (     REarth=6.3708e6  )
      integer         np, Lm2, Mm2
      common /bidul/  np
#if KEEP_SHALLOW
      logical         mask1(0:L,0:M), mask2(0:L,0:M), mask3(0:L,0:M)
      BIGREAL         depth1, depth2, depth3
      parameter    (  depth1=200, depth2=400, depth3=1000 )
#endif  /* KEEP_SHALLOW */
      parameter (     Lm2=L-2, Mm2=M-2  )
      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

#include "gridid.h"

      call get_xy
      call get_mn
      call get_h
      call get_rmask

      xmin = vmin(xr,Lp*Mp)
      ymin = vmin(yr,Lp*Mp)
      xmax = vmax(xr,Lp*Mp)
      ymax = vmax(yr,Lp*Mp)

      el = ymax-ymin
      xl = xmax-xmin

!     true for colour fill plots (false gives contours)
      colour = .true.

!     true to draw grid over colour plots
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

!     call ploth(gridid(1:40),colour,grover)
#endif  /* PLOTS */

#if KEEP_SHALLOW
      do j=0,M
        do i=0,L
          mask1(i,j) = .false.
          mask2(i,j) = .false.
          mask3(i,j) = .false.
        enddo
      enddo
      do j=0,M
        do i=0,L
          if (h(i,j) .lt. depth1) mask1(i,j) = .true.
          if (h(i,j) .lt. depth2) mask2(i,j) = .true.
          if (h(i,j) .lt. depth3) mask3(i,j) = .true.
        enddo
      enddo
#endif  /* KEEP_SHALLOW */

      print *,' Enter number of passes of filter:'
      read(5,*) npasses

      do 120 k=1,npasses

#if KEEP_SHALLOW
!  Keep shallow areas shallow
        if (k .le. npasses-2) then
          do j=0,M
            do i=0,L
              if (h(i,j) .gt. depth1 .and. mask1(i,j)) h(i,j) = depth1
              if (h(i,j) .gt. depth2 .and. mask2(i,j)) h(i,j) = depth2
              if (h(i,j) .gt. depth3 .and. mask3(i,j)) h(i,j) = depth3
            enddo
          enddo
        endif
#endif  /* KEEP_SHALLOW */

        call shapiro(h,Lp,Mp,2,wrk1,wrk2,mask_rho)
!       call ploth(gridid(1:40),colour,grover)
!  zero slope on h at boundaries
        do i=1,Lm
          h(i,0) = h(i,1)
          h(i,M) = h(i,Mm)
        enddo
        do j=0,M
          h(0,j) = h(1,j)
          h(L,j) = h(Lm,j)
        enddo

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
        print *,'depths: ', h(ii,jj), h(ii-1,jj), h(ii,jj-1)

 120  continue

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
 130  call crash('bathsuds: read or write error  ',0)
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
!  hacked about by John Wilkin 28 June 29, 1991

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

      integer    d, L, M, i, j, k
      BIGREAL    rvwgtx, rvwgty

      L = Lp-1
      M = Mp-1

      if (mod(N,2).ne.0) then
        print *,'N must be even in the shapiro filter'
#if NO_EXIT
        stop
#else
        call exit(1)
#endif  /* NO_EXIT */
      end if

!  Do the first y pass to initialize the temporary array

      do i=0,L
        do j=1,M-1
          if (mask(i,j) == 0 .or.                                       &
     &                (mask(i,j-1)==0 .and. mask(i,j+1)==0)) then
            tmp(i,j) = 0
          else if (mask(i,j-1) == 0 .and. mask(i,j+1) == 1) then
            tmp(i,j) = 0.25 * (u(i,j+1) - u(i,j))
          else if (mask(i,j-1) == 1 .and. mask(i,j+1) == 0) then
            tmp(i,j) = 0.25 * (u(i,j-1) - u(i,j))
          else
            tmp(i,j) = 0.25 * (u(i,j-1) + u(i,j+1) - 2*u(i,j))
          end if
        end do
      end do

!  Other passes in the y direction.

      do 120 k=4,N,2
        d = k/2
        do 110 j=d,M-d
        do 110 i=0,L
          tmp2(i,j) = - 0.25 * (tmp(i,j-1) + tmp(i,j+1)                 &
     &        - 2*tmp(i,j))
  110   continue
      do 120 j=d,M-d
      do 120 i=0,L
        tmp(i,j) = tmp2(i,j)
  120 continue

!  Add the changes to u

      do 130 j=1,M-1
      do 130 i=0,L
        u(i,j) = u(i,j) + rvwgty(i,j)*tmp(i,j)
  130 continue

!  Initialize tmp to filter in the x direction.

      do j=0,M
        do i=1,L-1
          if (mask(i,j) == 0 .or.                                       &
     &                (mask(i-1,j)==0 .and. mask(i+1,j)==0)) then
            tmp(i,j) = 0
          else if (mask(i-1,j) == 0 .and. mask(i+1,j) == 1) then
            tmp(i,j) = 0.25 * (u(i+1,j) - u(i,j))
          else if (mask(i-1,j) == 1 .and. mask(i+1,j) == 0) then
            tmp(i,j) = 0.25 * (u(i-1,j) - u(i,j))
          else
            tmp(i,j) = 0.25 * (u(i-1,j) + u(i+1,j) - 2*u(i,j))
          end if
        end do
      end do

!  Other x passes

      do 160 k=4,N,2
        d = k/2
      do 160 j=0,M
        do 150 i=d,L-d
          tmp2(i,j) = - 0.25 * (tmp(i-1,j) + tmp(i+1,j)                 &
     &        - 2*tmp(i,j))
  150   continue
      do 160 i=d,L-d
        tmp(i,j) = tmp2(i,j)
  160 continue

!  Add changes to u

      do 170 j=0,M
      do 170 i=1,L-1
        u(i,j) = u(i,j) + rvwgtx(i,j)*tmp(i,j)
  170 continue

      return
      end

! ********************************************************************

      BIGREAL function rvwgtx(i,j)
      integer     i, j
      BIGREAL     rv0, rvwscl,rvwmin,  rvwmax, dhdxx1, dhdxx2, dhdxx,   &
     &            rvwavg, rvwrng, h1, h2

#include "bathy.h"

      rv0    = 0.15
      rvwscl = 0.04
      rvwmin = 0.0
      rvwmax = 1.2
      rvwavg = 0.5*(rvwmin+rvwmax)
      rvwrng = 0.5*(rvwmax-rvwmin)

      h1 = max(5.0,h(i,j))
      h2 = max(5.0,h(i-1,j))
      dhdxx1 = abs((h1 - h2) / (h1 + h2))           &
     &         * mask_rho(i,j) * mask_rho(i-1,j)
      h2 = max(5.0,h(i+1,j))
      dhdxx2 = abs((h1 - h2) / (h1 + h2))           &
     &         * mask_rho(i,j) * mask_rho(i+1,j)
      dhdxx  = max(dhdxx1,dhdxx2)

      rvwgtx = rvwavg + rvwrng*(tanh((dhdxx-rv0)/rvwscl))
      if (h(i,j) .le. 5.0) rvwgtx = 0.0
      return
      end

! ********************************************************************

      BIGREAL function rvwgty(i,j)
      integer     i, j
      BIGREAL     rv0, rvwscl,rvwmin,  rvwmax, h1, h2,                  &
     &            rvwavg, rvwrng, dhdyy1, dhdyy2, dhdyy

#include "bathy.h"

      rv0    = 0.15
      rvwscl = 0.04
      rvwmin = 0.0
      rvwmax = 1.2
      rvwavg = 0.5*(rvwmin+rvwmax)
      rvwrng = 0.5*(rvwmax-rvwmin)
      h1 = max(5.0,h(i,j))
      h2 = max(5.0,h(i,j-1))
      dhdyy1 = abs((h1 - h2) / (h1 + h2))           &
     &         * mask_rho(i,j) * mask_rho(i,j-1)
      h2 = max(5.0,h(i,j+1))
      dhdyy2 = abs((h1 - h2) / (h1 + h2))           &
     &         * mask_rho(i,j) * mask_rho(i,j+1)
      dhdyy  = max(dhdyy1,dhdyy2)

      rvwgty = rvwavg + rvwrng*(tanh((dhdyy-rv0)/rvwscl))
      if (h(i,j) .le. 5.0) rvwgty = 0.0
      return
      end

! ********************************************************************

      BIGREAL function rv(i,j)
      integer   i, j

!     calculate the SPEM r-value
#include "bathy.h"
      BIGREAL   dhdxx, dhdyy, h1, h2, h3

      h1 = max(5.0,h(i,j))
      h2 = max(5.0,h(i-1,j))
      h3 = max(5.0,h(i,j-1))
      dhdxx  = abs((h1 - h2) / (h1 + h2))           &
     &         * mask_rho(i,j) * mask_rho(i-1,j)
      dhdyy  = abs((h1 - h3) / (h1 + h3))           &
     &         * mask_rho(i,j) * mask_rho(i,j-1)
      rv     = max(dhdxx,dhdyy)
      if (h(i,j) .le. 5.0) rv = 0.0
      return
      end
