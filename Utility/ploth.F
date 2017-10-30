#include "griddefs.h"
! ********************************************************************

      subroutine cpmpxy(imap,x,y,fx,fy)
      integer    imap
      real       x, y, fx, fy
#include "bathy.h"
#include "griddefs.h"

      if ( imap .eq. 3) then
        fx = x_v(int(x),int(y)) &
     &        + (x_v(int(x)+1,int(y))-x_v(int(x),int(y)))*(x-aint(x)) &
     &        + (x_v(int(x),int(y)+1)-x_v(int(x),int(y)))*(y-aint(y)) &
     &        + (x_v(int(x)+1,int(y)+1)-x_v(int(x),int(y)+1) &
     &        - x_v(int(x)+1,int(y))+x_v(int(x),int(y))) &
     &        *(x-aint(x))*(y-aint(y))
        fy = y_v(int(x),int(y)) &
     &        + (y_v(int(x)+1,int(y))-y_v(int(x),int(y)))*(x-aint(x)) &
     &        + (y_v(int(x),int(y)+1)-y_v(int(x),int(y)))*(y-aint(y)) &
     &        + (y_v(int(x)+1,int(y)+1)-y_v(int(x),int(y)+1) &
     &        - y_v(int(x)+1,int(y))+y_v(int(x),int(y))) &
     &        *(x-aint(x))*(y-aint(y))
      end if
      return
      end

!*******************************************************************c

      subroutine getxxyy
#include "bathy.h"
      integer    i, j

      do 100 i = 1,L
      do 100 j = 1,M
        x_v(i,j) = xp(i,j)
        y_v(i,j) = yp(i,j)
  100 continue
      do 110 i=1,L
        x_v(i,Mp) = x_v(i,M)
        y_v(i,Mp) = y_v(i,M)
  110 continue
      do 120 j=1,Mp
        x_v(Lp,j) = x_v(L,j)
        y_v(Lp,j) = y_v(L,j)
  120 continue
      return
      end

!*******************************************************************c

      subroutine getxyh
#include "bathy.h"
      integer    i, j

      do 100 i = 1,Lp
      do 100 j = 1,Mp
        x_v(i,j) = xr(i-1,j-1)
        y_v(i,j) = yr(i-1,j-1)
  100 continue
      do 110 i=1,Lp
        x_v(i,0) = x_v(i,1)
        y_v(i,0) = y_v(i,1)
        x_v(i,M+2) = x_v(i,Mp)
        y_v(i,M+2) = y_v(i,Mp)
  110 continue
      do 120 j=1,M+2
        x_v(0,j) = x_v(1,j)
        y_v(0,j) = y_v(1,j)
        x_v(L+2,j) = x_v(Lp,j)
        y_v(L+2,j) = y_v(Lp,j)
  120 continue
      return
      end

! ********************************************************************

      subroutine grdplt(x1,x2,y1,y2,gridid)

!  plots the psi points grid
!  if called without advancing frame this can be used to overlay the grid
!  on a colour filled contour plot

#include "bathy.h"
      character*40 gridid
      real         x1, x2, y1, y2
      integer      i,  j

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq (0.5,0.96,gridid,.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
#if DRAW_COASTS
      call mapdrw
#endif  /* DRAW_COASTS */
#if DBLEPREC
      do j = 1,Mm
        do i = 1,Lm
          if (mask_rho(i,j) .eq. 1.) then
            call frstpt(sngl(xp(i,j)),sngl(yp(i,j)))
            call vector(sngl(xp(i+1,j)),sngl(yp(i+1,j)))
            call vector(sngl(xp(i+1,j+1)),sngl(yp(i+1,j+1)))
            call vector(sngl(xp(i,j+1)),sngl(yp(i,j+1)))
            call vector(sngl(xp(i,j)),sngl(yp(i,j)))
          endif
        enddo
      enddo
#else
      do j = 1,Mm
        do i = 1,Lm
          if (mask_rho(i,j) .eq. 1.) then
            call frstpt(xp(i,j),yp(i,j))
            call vector(xp(i+1,j),yp(i+1,j))
            call vector(xp(i+1,j+1),yp(i+1,j+1))
            call vector(xp(i,j+1),yp(i,j+1))
            call vector(xp(i,j),yp(i,j))
          endif
        enddo
      enddo
#endif  /* DBLEPREC */

      return
      end

! ********************************************************************

      subroutine ploth(gridid,colour,grover)

#include "bathy.h"
#include "griddefs.h"
      real            dhdxx(L,0:M), dhdyy(0:L,M), tmp(L,M), &
     &                htmp(0:L,0:M)
      character*24    ltit, lnote, lnote2
      character*40    gridid
      logical         colour, grover, tallflg
      real            x1, x2, y1, y2, depmin, depmax, ratio, slpmax
      integer         i, j
      BIGREAL         a1, a2, av2, vmin, vmax
      real            vsmax

      av2(a1,a2) = 0.5*(a1+a2)

!  find shape of domain to plot
      if (xl .ge. el) then
        x1 = 0.05
        x2 = 0.95
        y1 = -.45*el/xl + .5
        y2 = y1 + el/xl*.9
        tallflg = .false.
      else
        if (colour) then
          y1 = 0.02
          y2 = 0.92
        else
          y1 = 0.04
          y2 = 0.94
        end if
        x1 = -.45*xl/el + .45
        x2 = x1 + xl/el*.9
        tallflg = .true.
      end if

!  draw the grid first
      call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
      call drawcoast
#endif  /* DRAW_COASTS */
      call frame

      write (ltit,100)
  100 format ('Bottom Topography')

      call getxyh
      depmin = vmin(h,Lp*Mp)
      write (lnote,120) depmin
  120 format ('MIN DEPTH =',f9.3)
      depmax = vmax(h,Lp*Mp)
      write (lnote2,131) depmax
  131 format ('MAX DEPTH =',f9.1)

      do j=0,M
        do i=0,L
          htmp(i,j) = h(i,j)
        enddo
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.50,.98,trim(ltit),.012,0.,0.)
      call plchhq(.85,.98,trim(lnote),.012,0.,0.)
      call plchhq(.85,.95,trim(lnote2),.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(htmp,Lp,Lp,Mp,16,lcflag,tallflg,.true.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(htmp,Lp,Lp,Mp,0.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      call getxxyy
      write (ltit,130)
  130 format ('Bottom Slope')

      do j=0,M
        do i=1,L
          dhdxx(i,j) = (h(i,j) - h(i-1,j))*av2(pm(i,j),pm(i-1,j))       &
     &         * mask_rho(i,j) * mask_rho(i-1,j)
        end do
      end do
      do j=1,M
        do i=0,L
          dhdyy(i,j) = (h(i,j) - h(i,j-1))*av2(pn(i,j),pn(i,j-1))       &
     &         * mask_rho(i,j) * mask_rho(i,j-1)
        end do
      end do
      do j=1,M
        do i=1,L
          tmp(i,j) = sqrt((0.5*(dhdxx(i,j)+dhdxx(i,j-1)))**2 +          &
     &        (0.5*(dhdyy(i,j)+dhdyy(i-1,j)))**2)
        end do
      end do

      slpmax = vsmax(tmp,L*M)
      write (lnote,170) slpmax
  170 format ('MAX SLOPE =',f9.3)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq (0.5,0.98,trim(ltit),.012,0.,0.)
      call plchhq(.85,.98,trim(lnote),.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp,L,L,M,0.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

#ifndef NEGATIVE
      write (ltit,180)
  180 format ('r-Value')

      do 190  j=0,M
      do 190  i=1,L
        dhdxx(i,j) = abs((h(i,j) - h(i-1,j)) &
     &        /     (h(i,j) + h(i-1,j))) &
     &        * mask_rho(i,j) * mask_rho(i-1,j)
  190 continue
      do 200  j=1,M
      do 200  i=0,L
        dhdyy(i,j) = abs((h(i,j) - h(i,j-1)) &
     &        /     (h(i,j) + h(i,j-1))) &
     &        * mask_rho(i,j) * mask_rho(i,j-1)
  200 continue
      do 210  j=1,M
      do 210  i=1,L
        tmp(i,j) = max(max(dhdxx(i,j),dhdxx(i,j-1)), &
     &        max(dhdyy(i,j),dhdyy(i-1,j)))
  210 continue

      ratio = vsmax(tmp,L*M)
      write (lnote,220) ratio
  220 format ('MAX RATIO =',f9.3)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq (0.5,0.98,trim(ltit),.012,0.,0.)
      call plchhq(.85,.98,trim(lnote),.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp,L,L,M,0.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame
#endif

      return
      end

! ***********************************************************

      real function vsmax(vect,N)
      integer   N
      real      vect(N)
      real      tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = max(tmp,vect(i))
      enddo
      vsmax = tmp
      return
      end

      real function vsmin(vect,N)
      integer   N
      real      vect(N)
      real      tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = min(tmp,vect(i))
      enddo
      vsmin = tmp
      return
      end

