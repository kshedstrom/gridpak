      program spheres

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!  This program reads in the x,y grid, the lat,long grid,
!  and the topography from the netCDF grid file.  It uses the
!  latitude information to calculate f = 2*Omega*sin(lat). The
!  results are placed into the netCDF grid file.
!
!  Logicals:
!    hflat  - true for constant depth of h0 - don't read in topography.
!    sphere - true for spherical geometry in pm's and pn's or false
!       for planar geometry.
!    colour - true for colour fill plots (false gives contours)
!    grover - true to draw grid over colour filled plots
!    readlat - true to read in the lat, long file produced by tolat.
! *******************************************************************

#include "griddefs.h"
#include "bathy.h"
#include "ncgrid.h"
#include "proj.h"
      BIGREAL         h0, f0, beta
      parameter (     h0 = 800.   )
      parameter (     f0=9.e-4, beta=2.e-11  )
#if DRAW_COASTS
!  Earth angular frequency and radius of sphere of equal volume
      real            U1, U2, V1, V2, udeg, uscale
      BIGREAL         enlarge
#endif  /* DRAW_COASTS */
      real            DTOR, RTOD
      BIGREAL         Aomega, REarth, pi
      parameter (     Aomega=7.292115e-5 , REarth=6.3708e6  )
      parameter (     pi = 3.14159265  )
      logical         hflat, colour, grover, readlat
      BIGREAL         el2, sum, dxmin, dist, galpha
      integer         i, j, imap
      BIGREAL         a1, a2, vmin, vmax

      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

#include "sphereflags.h"
#include "gridid.h"

      call get_xy
      call get_mn
      call get_rmask

      if (readlat) then
        call get_lat
!  calculate the angle between the xi-eta grid and the lon-lat grid
!  at rho points.

        do j = 1, Mm
          do i = 1, Lm
            a1 = lat_u(i+1, j) - lat_u(i, j)
            a2 = lon_u(i+1, j) - lon_u(i, j)
            if (abs(a2) .gt. 180.) then
              if (a2 .lt. -180. ) then
                a2 = a2 + 360.
              else
                a2 = a2 - 360.
              endif
            endif
            a2 = a2 * cos(0.5*DTOR*(lat_u(i, j) + lat_u(i+1, j)))
            angle(i, j) = atan2(a1, a2)
          enddo
        enddo
        do j = 1, Mm
          do i = 1, Lm
            a2 = lat_v(i, j) - lat_v(i, j+1)
            a1 = (lon_v(i, j) - lon_v(i, j+1))
            if (abs(a1) .gt. 180.) then
              if (a1 .lt. -180. ) then
                a1 = a1 + 360.
              else
                a1 = a1 - 360.
              endif
            endif
            a1 = a1 * cos(0.5*DTOR*(lat_v(i, j) + lat_v(i, j+1)))
            angle(i, j) = 0.5*(angle(i, j) + atan2(a1, -a2))
          enddo
        enddo
      else
!  calculate the angle between the xi-eta grid and the lon-lat grid
!  at rho points.
        do j = 1, Mm
          do i = 1, Lm
            a1 = yu(i+1, j) - yu(i, j)
            a2 = xu(i+1, j) - xu(i, j)
            angle(i, j) = atan2(a1, a2)
            a2 = yv(i, j) - yv(i, j+1)
            a1 = xv(i, j) - xv(i, j+1)
            angle(i, j) = 0.5*(angle(i, j) + atan2(a1, -a2))
          enddo
        enddo
      endif

      do i = 1, Lm
        angle(i, 0) = angle(i, 1)
        angle(i, M) = angle(i, Mm)
      enddo

      do j = 0, M
        angle(0, j) = angle(1, j)
        angle(L, j) = angle(Lm, j)
      enddo

      if (hflat) then
        do j=0,M
          do i=0,L
            h(i,j) = h0
          enddo
        enddo
      else
        call get_h
      end if

! If on a sphere then compute spherical pm,pn factors
      if (spherical == 1 .and. readlat) then

        do j = 1,Mm
          do i = 1,Lm
            call geodesic_dist(lon_u(i,j),lat_u(i,j),lon_u(i+1,j), &
     &                     lat_u(i+1,j),1,dist,galpha)
            pm(i,j) = 1./dist
            call geodesic_dist(lon_v(i,j),lat_v(i,j),lon_v(i,j+1), &
     &                     lat_v(i,j+1),1,dist,galpha)
            pn(i,j) = 1./dist
          enddo
        enddo

!  M,N factors outside the boundaries
        do j=1,Mm
          pm(0,j) = pm(1,j)
          pn(0,j) = pn(1,j)
          pm(L,j) = pm(Lm,j)
          pn(L,j) = pn(Lm,j)
        enddo
        do i=0,L
          pm(i,0) = pm(i,1)
          pn(i,0) = pn(i,1)
          pm(i,M) = pm(i,Mm)
          pn(i,M) = pn(i,Mm)
        enddo

!  Compute dndx,dmde
        do i=1,Lm
          do j=1,Mm
            dndx(i,j) = (1./pn(i+1,j) - 1./pn(i-1,j))/2.
          enddo
        enddo
        do i=1,Lm
          do j=1,Mm
            dmde(i,j) = (1./pm(i,j+1) - 1./pm(i,j-1))/2.
          enddo
        enddo
        do j=1,Mm
          dndx(0,j) = 0.
          dmde(0,j) = 0.
          dndx(L,j) = 0.
          dmde(L,j) = 0.
        enddo
        do i=0,L
          dndx(i,0) = 0.
          dmde(i,0) = 0.
          dndx(i,M) = 0.
          dmde(i,M) = 0.
        enddo
      end if

!  Coriolis parameter
      if (readlat) then
        do i=0,L
          do j=0,M
            f(i,j) = 2*Aomega*sin(DTOR*lat_rho(i,j))
          enddo
        enddo
      else
        el2 = yp(1,M/2)
        do i=0,L
          do j=0,M
            f(i,j) = f0 + beta*(yr(1,j) - el2)
          enddo
        enddo
      end if

!  Output for the SPEM

      call wrt_fhmn(spherical,f0,beta)

#if PLOTS
!  Plots just to be sure

      xmin = vmin(xr,Lp*Mp)
      ymin = vmin(yr,Lp*Mp)
      xmax = vmax(xr,Lp*Mp)
      ymax = vmax(yr,Lp*Mp)

      el = ymax-ymin
      xl = xmax-xmin

      call start_plot
!  set foreground to black, background to white
      if (colour) then
        call gscr(1,0,1.,1.,1.)
        call gscr(1,1,0.,0.,0.)
      end if

#if DRAW_COASTS
!  put things on EZMAP space:
      call mapsti('GR',JGRD)
      call mapstc('OU','NO')
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

      imap = 3
      call cpseti('MAP',imap)
      call getxxyy
      call pcseti('QUALITY',1)
      call ploth(gridid(1:40),colour,grover)
      call plotmnf(gridid(1:40),colour,grover)

      call end_plot
#endif  /* PLOTS */

      sum = 0
      do j=1,Mm
        do i=1,Lm
          sum = sum + 1/pm(i,j)
        enddo
      enddo
      print *,'average dx ',sum/(Lm*Mm)
      dxmin = vmax(pm,Lp*Mp)
      dxmin = 1/dxmin
      print *, 'minimum dx = ',dxmin
      sum = 0
      do j=1,Mm
        do i=1,Lm
          sum = sum + 1/pn(i,j)
        enddo
      enddo
      print *,'average dy ',sum/(Lm*Mm)
      dxmin = vmax(pn,Lp*Mp)
      dxmin = 1/dxmin
      print *, 'minimum dy = ',dxmin

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
  200 call crash('sphere: read or write error', 0)
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
      call exit(ierr)
#endif  /* NO_EXIT */
      return
      end

!****************************************************************

      subroutine plotmnf(gridid,colour,grover)
#include "bathy.h"
      real            tmp(L,M), tmp2(L,M), tmp3(L,M), htmp(0:L,0:M)
      character*20    ltit
      character*24    lnote, lnote2
      character*40    gridid
      logical         colour, grover, tallflg
      real            x1, x2, y1, y2, dxmin, dxmax, dymin, dymax
      integer         i, j
      real            vsmin, vsmax

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

      call getxxyy

      write (ltit,100)
  100 format ('DX')

      do j=1,M
        do i=1,L
          tmp(i,j) = 4./(pm(i,j) + pm(i-1,j) + pm(i,j-1) + &
     &        pm(i-1,j-1))
        enddo
      enddo

      dxmin = vsmin(tmp,L*M)
      write (lnote,120) dxmin
  120 format ('MIN DX =',f9.0)
      dxmax = vsmax(tmp,L*M)
      write (lnote2,130) dxmax
  130 format ('MAX DX =',f9.0)

      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,ltit,.012,0.,0.)
      if (colour) then
        call plchhq(.85,.98,trim(lnote),.012,0.,0.)
        call plchhq(.85,.95,trim(lnote2),.012,0.,0.)
      endif
      call set (x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp,L,L,M,-20.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      write (ltit,140)
  140 format ('DY')

      do j=1,M
        do i=1,L
          tmp2(i,j) = 4./(pn(i,j) + pn(i-1,j) + pn(i,j-1) + &
     &        pn(i-1,j-1))
        enddo
      enddo

      dymin = vsmin(tmp2,L*M)
      write (lnote,150) dymin
  150 format ('MIN DY =',f9.0)
      dymax = vsmax(tmp2,L*M)
      write (lnote2,155) dymax
  155 format ('MAX DY =',f9.0)

      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,ltit,.012,0.,0.)
      if (colour) then
        call plchhq(.85,.98,trim(lnote),.012,0.,0.)
        call plchhq(.85,.95,trim(lnote2),.012,0.,0.)
      endif
      call set (x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp2,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp2,L,L,M,-20.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      write (ltit,160)
  160 format ('RMS gridsize (km)')

      do j=1,M
        do i=1,L
          tmp3(i,j)=sqrt(tmp(i,j)**2 + tmp2(i,j)**2)/sqrt(2.0d0)*0.001d0
        enddo
      enddo

      dymin = vsmin(tmp3,L*M)
      write (lnote,161) dymin
  161 format ('MIN SIZE =',f9.2)
      dymax = vsmax(tmp3,L*M)
      write (lnote2,162) dymax
  162 format ('MAX SIZE =',f9.2)

      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,ltit,.012,0.,0.)
      if (colour) then
        call plchhq(.85,.98,trim(lnote),.012,0.,0.)
        call plchhq(.85,.95,trim(lnote2),.012,0.,0.)
      endif
      call set (x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp3,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp3,L,L,M,-20.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      write (ltit,170)
  170 format ('ASPECT RATIO')

      do j=1,M
        do i=1,L
          tmp3(i,j) = tmp(i,j)/tmp2(i,j)
        enddo
      enddo

      dymin = vsmin(tmp3,L*M)
      write (lnote,171) dymin
  171 format ('MIN RATIO =',f9.2)
      dymax = vsmax(tmp3,L*M)
      write (lnote2,172) dymax
  172 format ('MAX RATIO =',f9.2)

      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,ltit,.012,0.,0.)
      if (colour) then
        call plchhq(.85,.98,trim(lnote),.012,0.,0.)
        call plchhq(.85,.95,trim(lnote2),.012,0.,0.)
      endif
      call set (x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(tmp3,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp3,L,L,M,-20.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      write (ltit,180)
  180 format ('F')

      do j=1,M
        do i=1,L
          tmp(i,j) = 0.25*(f(i,j) + f(i-1,j) + f(i,j-1) + f(i-1,j-1))
        enddo
      enddo

      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq (0.5,0.98,ltit(1:1),.012,0.,0.)
      call set (x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill (tmp,L,L,M,8,lcflag,tallflg,.false.)
        if (grover) call grdplt(x1,x2,y1,y2,gridid)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      else
        call cpshift(tmp,L,L,M,-20.,.false.,1.)
#if DRAW_COASTS
        call drawcoast
#endif  /* DRAW_COASTS */
      endif
      call frame

      write (ltit,190)
  190 format ('RHO MASK')

      call getxyh
      do j=0,M
        do i=0,L
          htmp(i,j) = mask_rho(i,j)
        enddo
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ltit(1:8),.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(htmp,Lp,Lp,Mp,2,lcflag,tallflg,.false.)
        call grdplt(x1,x2,y1,y2,gridid)
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

      write (ltit,200)
  200 format ('ANGLE')

      call getxyh
      do j=0,M
        do i=0,L
          htmp(i,j) = angle(i,j)
        enddo
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ltit(1:5),.012,0.,0.)
      call set(x1,x2,y1,y2,xmin,xmax,ymin,ymax,1)
      if (colour) then
        call cpsfill(htmp,Lp,Lp,Mp,8,lcflag,tallflg,.false.)
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

      return
      end
