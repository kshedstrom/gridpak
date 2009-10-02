      program Wgrid

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
! *******************************************************************
!     This version reads digitized rectangular boundaries, uses
!     splines to interpolate to the grid resolution, and fills in
!     the interior.
!
!     Direct any question/problems (and reports of errors!) to:
!     Kate Hedstrom
!     kate@ahab.rutgers.edu
!
!  IBIG is the largest number of points to be read in for one
!    boundary

#include "griddefs.h"
#include "grid.h"
#include "ncgrid.h"
      integer         i, j
      BIGCOMPLEX      z(M2+L2+M2+L2)
      BIGREAL         stmp, area, tarea
      BIGREAL         xb(M2+L2+M2+L2),yb(M2+L2+M2+L2), &
     &                x(0:L2,0:M2),y(0:L2,0:M2)
      character*40    lbl

! label for grid
#include "gridid.h"
      write(lbl,100) gridid(1:40)
  100 format(40a)

!  original distribution of x,y points is preserved on boundary kb1
!  and kb2:
      if(kb1.ne.1.and.kb1.ne.3) &
     &        call crash('boundary index kb1 must be 1 or 3', kb1)
      if(kb2.ne.2.and.kb2.ne.4) &
     &        call crash('boundary index kb2 must be 2 or 4', kb2)

!  set up boundary spline interpolation arrays
      call readbndy

!  initialize vector z (complex) with contour of physical boundary
#if XPOTS1
      read(ipot1,*) stmp
#endif  /* XPOTS1 */
      do i=1,N1
#if XPOTS1
        read(ipot1,*) stmp
        call z1(stmp,xb(i),yb(i),z(i))
#else
        stmp = FLoaT(i)/FLoaT(N1)
        call z1(stmp,xb(i),yb(i),z(i))
#endif  /* XPOTS1 */
      enddo
#if XPOTS2
      read(ipot2,*) stmp
#endif  /* XPOTS2 */
      do i=N1+1,N2
#if XPOTS2
        read(ipot2,*) stmp
        call z2(stmp,xb(i),yb(i),z(i))
#else
        stmp = FLoaT(i-N1)/FLoaT(N2-N1)
        call z2(stmp,xb(i),yb(i),z(i))
#endif  /* XPOTS2 */
      enddo
#if XPOTS3
      read(ipot3,*) stmp
#endif  /* XPOTS3 */
      do i=N2+1,N3
#if XPOTS3
        read(ipot3,*) stmp
        call z3(stmp,xb(i),yb(i),z(i))
#else
        stmp = FLoaT(i-N2)/FLoaT(N3-N2)
        call z3(stmp,xb(i),yb(i),z(i))
#endif  /* XPOTS3 */
      enddo
#if XPOTS4
      read(ipot4,*) stmp
#endif  /* XPOTS4 */
      do i=N3+1,N4
#if XPOTS4
        read(ipot4,*) stmp
        call z4(stmp,xb(i),yb(i),z(i))
#else
        stmp = FLoaT(i-N3)/FLoaT(N4-N3)
        call z4(stmp,xb(i),yb(i),z(i))
#endif  /* XPOTS4 */
      enddo

!  Set boundary values of the grid
      do i=1,N1
        x(0,N1-i)  = xb(i)
        y(0,N1-i)  = yb(i)
      enddo
      do i=N1+1,N2
        x(i-N1,0)  = xb(i)
        y(i-N1,0)  = yb(i)
      enddo
      do i=N2+1,N3
        x(L2,i-N2) = xb(i)
        y(L2,i-N2) = yb(i)
      enddo
      do i=N3+1,N4
        x(N4-i,M2) = xb(i)
        y(N4-i,M2) = yb(i)
      enddo

! interior of the grid
      if (kb1 .eq. 1 .and. kb2 .eq. 2) then
        do j=1,M2
          do i=1,L2
            x(i,j) = x(i,0) + x(0,j) - x(0,0)
            y(i,j) = y(i,0) + y(0,j) - y(0,0)
          enddo
        enddo
      else if (kb1 .eq. 1 .and. kb2 .eq. 4) then
        do j=0,M2-1
          do i=1,L2
            x(i,j) = x(i,M2) + x(0,j) - x(0,M2)
            y(i,j) = y(i,M2) + y(0,j) - y(0,M2)
          enddo
        enddo
      else if (kb1 .eq. 3 .and. kb2 .eq. 2) then
        do j=1,M2
          do i=0,L2-1
            x(i,j) = x(i,0) + x(L2,j) - x(L2,0)
            y(i,j) = y(i,0) + y(L2,j) - y(L2,0)
          enddo
        enddo
      else if (kb1 .eq. 3 .and. kb2 .eq. 4) then
        do j=0,M2-1
          do i=0,L2-1
            x(i,j) = x(i,M2) + x(L2,j) - x(L2,M2)
            y(i,j) = y(i,M2) + y(L2,j) - y(L2,M2)
          enddo
        enddo
      endif

!  Compute pm,pn factors
      do j = 1,Mm
        do i = 1,Lm
          pm(i,j) = 1.0/sqrt( ( x(2*i,2*j-1)-x(2*i-2,2*j-1) )**2 &
     &        + ( y(2*i,2*j-1)-y(2*i-2,2*j-1) )**2 )
          pn(i,j) = 1.0/sqrt( ( x(2*i-1,2*j)-x(2*i-1,2*j-2) )**2 &
     &        + ( y(2*i-1,2*j)-y(2*i-1,2*j-2) )**2 )
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

!  Split up grid solution into separate arrays for the coordinates
!  of the four locations on the Arakawa C grid corresponding to rho,
!  psi, u and v points.

!  psi points
      do j=1,M
        do i=1,L
          xp(i,j) = x(2*i-2,2*j-2)
          yp(i,j) = y(2*i-2,2*j-2)
        enddo
      enddo

      if (subset) then
         write (iout3,*) Mnorth-Msouth+1
         do j=Mnorth,Msouth,-1
            write (iout3,*) xp(Lwest,j),yp(Lwest,j)
         enddo
         write (iout3,*) Least-Lwest+1
         do i=Lwest,Least
            write (iout3,*) xp(i,Msouth),yp(i,Msouth)
         enddo
         write (iout3,*) Mnorth-Msouth+1
         do j=Msouth,Mnorth
            write (iout3,*) xp(Least,j),yp(Least,j)
         enddo
         write (iout3,*) Least-Lwest+1
         do i=Least,Lwest,-1
            write (iout3,*) xp(i,Mnorth),yp(i,Mnorth)
         enddo
      end if

      xl = xp(L,1) - xp(1,1)
      el = yp(1,M) - yp(1,1)

!  u points
      do j=1,Mm
        do i=1,L
          xu(i,j) = x(2*i-2,2*j-1)
          yu(i,j) = y(2*i-2,2*j-1)
        enddo
      enddo
      do i=1,L
        xu(i,0) = 2*xp(i,1) - xu(i,1)
        xu(i,M) = 2*xp(i,M) - xu(i,Mm)
        yu(i,0) = 2*yp(i,1) - yu(i,1)
        yu(i,M) = 2*yp(i,M) - yu(i,Mm)
      enddo

!  vpoints
      do j=1,M
        do i=1,Lm
          xv(i,j) = x(2*i-1,2*j-2)
          yv(i,j) = y(2*i-1,2*j-2)
        enddo
      enddo
      do j=1,M
        xv(0,j) = 2*xp(1,j) - xv(1,j)
        xv(L,j) = 2*xp(L,j) - xv(Lm,j)
        yv(0,j) = 2*yp(1,j) - yv(1,j)
        yv(L,j) = 2*yp(L,j) - yv(Lm,j)
      enddo

!  rho points
      do j=1,Mm
        do i=1,Lm
          xr(i,j) = x(2*i-1,2*j-1)
          yr(i,j) = y(2*i-1,2*j-1)
        enddo
      enddo
      do j=1,Mm
        xr(L,j) = 2*xu(L,j) - xr(Lm,j)
        xr(0,j) = 2*xu(1,j) - xr(1,j)
        yr(L,j) = 2*yu(L,j) - yr(Lm,j)
        yr(0,j) = 2*yu(1,j) - yr(1,j)
      enddo
      do i=0,L
        xr(i,M) = 2*xv(i,M) - xr(i,Mm)
        yr(i,M) = 2*yv(i,M) - yr(i,Mm)
        xr(i,0) = 2*xv(i,1) - xr(i,1)
        yr(i,0) = 2*yv(i,1) - yr(i,1)
      enddo

!  Diagnostics:
!  Compute area of domain from m,n factors
      area = 0.
      do j = 1,Mm
        do i = 1,Lm
          area = area + 1./(pm(i,j)*pn(i,j))
        enddo
      enddo

      write(6,460)area
  460 format(5x,'area         ',g10.4)

!  Output solution
!  (This output format parallels the input format assumed by
!  subroutine getgrid in pemodel)
      call checkdefs
      call def_grid
#if PLOTS
      call start_plot
      call gridplot(lbl,area,tarea)
      call end_plot
#endif  /* PLOTS */
#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
      end

! *******************************************************************

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

! *******************************************************************

#if PLOTS
      subroutine gridplot(lbl,area,tarea)
#include "grid.h"
      character*20    lbl
      BIGREAL         area, tarea

      character*5     lbls
      integer         i, j
      SMALLREAL       tmp(0:L,0:M), xlow, xhigh, ylow, yhigh
      SMALLREAL       x1, x2, y1, y2
      BIGREAL         vmin, vmax
      BIGREAL         dxmin

      call pcseti('QUALITY',1)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,lbl,.012,0.,0.)
      write(lbl,100)area
  100 format(' Sum 1/MN=',g10.4)
      call plchhq(.5,.03,lbl,.012,0.,0.)

!  Draw the psi points grid (boundaries coincide with channel walls)
      xlow  = vmin(xp,L*M)
      xhigh = vmax(xp,L*M)
      yhigh = vmax(yp,L*M)
      ylow  = vmin(yp,L*M)
      xl = xhigh-xlow
      el = yhigh-ylow
      if (xl .ge. el) then
        x1 = 0.05
        x2 = 0.95
        y1 = -.45*el/xl + .5
        y2 = y1 + el/xl*.9
      else
        y1 = 0.05
        y2 = 0.95
        x1 = -.45*xl/el + .5
        x2 = x1 + xl/el*.9
      end if
      call set(x1,x2,y1,y2,xlow,xhigh,ylow,yhigh,1)
#if DBLEPREC
      do j = 1,M
        call frstpt(sngl(xp(1,j)),sngl(yp(1,j)))
        do i = 2,L
          call vector(sngl(xp(i,j)),sngl(yp(i,j)))
        enddo
      enddo
      do i = 1,L
        call frstpt(sngl(xp(i,1)),sngl(yp(i,1)))
        do j = 2,M
          call vector(sngl(xp(i,j)),sngl(yp(i,j)))
        enddo
      enddo
#else
      do j = 1,M
        call frstpt(xp(1,j),yp(1,j))
        do i = 2,L
          call vector(xp(i,j),yp(i,j))
        enddo
      enddo
      do i = 1,L
        call frstpt(xp(i,1),yp(i,1))
        do j = 2,M
          call vector(xp(i,j),yp(i,j))
        enddo
      enddo
#endif  /* DBLEPREC */
      call frame

!  Plot m and n
      write(lbls,130)
  130 format('DX')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,lbls(1:2),.012,0.,0.)
      do j=0,M
        do i=0,L
          tmp(i,j) = 1/pm(i,j)
        enddo
      enddo
      dxmin = vmax(pm,Lp*Mp)
      dxmin = 1/dxmin
      print *, 'minimum dx = ',dxmin
      call set(x1,x2,y1,y2,1.,float(Lp),1.,float(Mp),1)
      call perim(1,L,1,M)
      call cpseti('SET',0)
      call cpcnrc(tmp,Lp,Lp,Mp,0.0,0.0,0.0,1,-1,-682)
      call frame
      write(lbls,150)
  150 format('DY')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,lbls(1:2),.012,0.,0.)
      do j=0,M
        do i=0,L
          tmp(i,j) = 1/pn(i,j)
        enddo
      enddo
      dxmin = vmax(pn,Lp*Mp)
      dxmin = 1/dxmin
      print *, 'minimum dy = ',dxmin
      call set(x1,x2,y1,y2,1.,float(Lp),1.,float(Mp),1)
      call perim(1,L,1,M)
      call cpseti('SET',0)
      call cpcnrc(tmp,Lp,Lp,Mp,0.0,0.0,0.0,1,-1,-682)
      call frame

      return
      end
#endif  /* PLOTS */

! *******************************************************************

      subroutine readbndy
#include "grid.h"
      BIGREAL         arcs(IBIG), huge
      parameter    (  huge = 1.e35   )
      integer         i

      if (rbnd1) then
!  Read x,y data points describing boundary 1
        read(5,*) nb1pts
        if (nb1pts .gt. IBIG) then
          call crash('readbndy: Need to make IBIG at least', nb1pts)
        end if
        do i=1,nb1pts
          read(5,*) x1spl(i),y1spl(i)
        enddo
        if (.not. even1) then
          do i=1,nb1pts
            s1spl(i) = (i-1.)/(nb1pts-1.)
          enddo
        else
          arcs(1) = 0
          do i = 2,nb1pts
            arcs(i) = arcs(i-1) + sqrt((x1spl(i)-x1spl(i-1))**2 + &
     &        (y1spl(i)-y1spl(i-1))**2)
          enddo
          do i=1,nb1pts
            s1spl(i) = arcs(i)/arcs(nb1pts)
          enddo
        end if

        call spline(s1spl,x1spl,nb1pts,huge,huge,b1spl)
        call spline(s1spl,y1spl,nb1pts,huge,huge,c1spl)
      end if

      if (rbnd2) then
!  Read x,y data points describing boundary 2
        read(5,*) nb2pts
        if (nb2pts .gt. IBIG) then
          call crash('readbndy: Need to make IBIG at least', nb2pts)
        end if
        do i=1,nb2pts
          read(5,*) x2spl(i),y2spl(i)
        enddo
        if (.not. even2) then
          do i=1,nb2pts
            s2spl(i) = (i-1.)/(nb2pts-1.)
          enddo
        else
          arcs(1) = 0
          do i = 2,nb2pts
            arcs(i) = arcs(i-1) + sqrt((x2spl(i)-x2spl(i-1))**2 + &
     &        (y2spl(i)-y2spl(i-1))**2)
          enddo
          do i=1,nb2pts
            s2spl(i) = arcs(i)/arcs(nb2pts)
          enddo
        end if

        call spline(s2spl,x2spl,nb2pts,huge,huge,b2spl)
        call spline(s2spl,y2spl,nb2pts,huge,huge,c2spl)
      end if

      if (rbnd3) then
!  Read x,y data points describing boundary 3
        read(5,*) nb3pts
        if (nb3pts .gt. IBIG) then
          call crash('readbndy: Need to make IBIG at least', nb3pts)
        end if
        do i=1,nb3pts
          read(5,*) x3spl(i),y3spl(i)
        enddo
        if (.not. even3) then
          do i=1,nb3pts
            s3spl(i) = (i-1.)/(nb3pts-1.)
          enddo
        else
          arcs(1) = 0
          do i = 2,nb3pts
            arcs(i) = arcs(i-1) + sqrt((x3spl(i)-x3spl(i-1))**2 + &
     &        (y3spl(i)-y3spl(i-1))**2)
          enddo
          do i=1,nb3pts
            s3spl(i) = arcs(i)/arcs(nb3pts)
          enddo
        end if

        call spline(s3spl,x3spl,nb3pts,huge,huge,b3spl)
        call spline(s3spl,y3spl,nb3pts,huge,huge,c3spl)
      end if

      if (rbnd4) then
!  Read x,y data points describing boundary 4
        read(5,*) nb4pts
        if (nb4pts .gt. IBIG) then
          call crash('readbndy: Need to make IBIG at least', nb4pts)
        end if
        do i=1,nb4pts
          read(5,*) x4spl(i),y4spl(i)
        enddo
        if (.not. even4) then
          do i=1,nb4pts
            s4spl(i) = (i-1.)/(nb4pts-1.)
          enddo
        else
          arcs(1) = 0
          do i = 2,nb4pts
            arcs(i) = arcs(i-1) + sqrt((x4spl(i)-x4spl(i-1))**2 + &
     &        (y4spl(i)-y4spl(i-1))**2)
          enddo
          do i=1,nb4pts
            s4spl(i) = arcs(i)/arcs(nb4pts)
          enddo
        end if

        call spline(s4spl,x4spl,nb4pts,huge,huge,b4spl)
        call spline(s4spl,y4spl,nb4pts,huge,huge,c4spl)
      end if
      return
      end

! *******************************************************************

      subroutine spline(x,y,NN,yp1,ypn,y2)

!     The following two subroutines are used to perform the cubic spline
!     interpolation required to match up the distribution of points on
!     opposite sides of the transformed plane rectangle.  The routines
!     are taken from  Press,W.H., B.P.Flannery, S.A.Teukolsky and
!     W.T.Vetterling: "Numerical Recipes, the Art of Scientific Computing"
!     Cambridge University Press, 1986.
#include "grid.h"
      integer         NN
      BIGREAL         x(NN), y(NN), y2(NN), yp1, ypn

      integer         nmax
      parameter (     nmax=M2+L2             )
      BIGREAL         u(nmax), sig, p, qn, un
      integer         i, k

      if (yp1 .gt. .99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,NN-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn .gt. .99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(NN)-x(NN-1)))*(ypn-(y(NN)-y(NN-1))/(x(NN)-x(NN-1)))
      endif
      y2(NN)=(un-qn*u(NN-1))/(qn*y2(NN-1)+1.)
      do k=NN-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end

! *******************************************************************

      subroutine splint(xa,ya,y2a,NN,x,y)
      integer    NN
      BIGREAL    xa(NN), ya(NN), y2a(NN), x, y

      integer    klo, khi, k
      BIGREAL    h, a, b

      klo=1
      khi=NN
  100 if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 100
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) call crash('splint: bad xa input', 0)
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ &
     &        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

! *******************************************************************

      subroutine z1(s1,x,y,z)
      BIGREAL    s1, s2, s3, s4, x, y
      BIGCOMPLEX z
#include "grid.h"

!     Subroutines which specify the boundaries of the physical
!     domain.  They are defined as functions of the variable s which
!     ranges from 0 to 1 on each boundary, proceeding anti-clockwise.

      call splint(s1spl,x1spl,b1spl,nb1pts,s1,x)
      call splint(s1spl,y1spl,c1spl,nb1pts,s1,y)
#if DBLEPREC
      z = dcmplx(x,y)
#else
      z = cmplx(x,y)
#endif  /* DBLEPREC */
      return

      entry z2(s2,x,y,z)
      call splint(s2spl,x2spl,b2spl,nb2pts,s2,x)
      call splint(s2spl,y2spl,c2spl,nb2pts,s2,y)
#if DBLEPREC
      z = dcmplx(x,y)
#else
      z = cmplx(x,y)
#endif  /* DBLEPREC */
      return

      entry z3(s3,x,y,z)
      call splint(s3spl,x3spl,b3spl,nb3pts,s3,x)
      call splint(s3spl,y3spl,c3spl,nb3pts,s3,y)
#if DBLEPREC
      z = dcmplx(x,y)
#else
      z = cmplx(x,y)
#endif  /* DBLEPREC */
      return

      entry z4(s4,x,y,z)
      call splint(s4spl,x4spl,b4spl,nb4pts,s4,x)
      call splint(s4spl,y4spl,c4spl,nb4pts,s4,y)
#if DBLEPREC
      z = dcmplx(x,y)
#else
      z = cmplx(x,y)
#endif  /* DBLEPREC */
      return
      end

