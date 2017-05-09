      program gridmain

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kshedstrom@alaska.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!     This version reads digitized coastal and offshore boundaries and
!     computes an orthogonal curvilinear coordinate grid for use with
!     the Haidvogel Semi-spectral Primitive Equation Model (SPEM).
!     The algorithm used is from David C.Ives and Robert M.Zacharias:
!     "Conformal Mapping and Orthogonal Grid Generation",
!     paper no. 87-2057, AIAA/SAE/ASME/ASEE 23rd Joint Propulsion
!     Conference, San Diego, California, June 1987.
!
!     It now uses sepeli only.
!
!     This version uses double precision complex arithmetic.  If your
!     compiler does not support double precision complex, undefine
!     DCOMPLEX in griddefs.h. Conversion
!     to single precision will produce a significant degradation in
!     performance.  It will probably cause 'tracking error' problems in
!     the conformal mapping calculation unless the machine you are using
!     (e.g. a Cray) has a greater than 32-bit word length.

!     Direct any question/problems (and reports of errors!) to:
!     Kate Hedstrom
!     kshedstrom@alaska.edu
!     email: wilkin@flood.ml.csiro.au

!  ITMAX is the number of iterations to perform
!  IBIG is the largest number of points to be read in for one
!  boundary

      implicit none
#include "griddefs.h"
#include "grid.h"
#include "ncgrid.h"
#if DCOMPLEX
      complex*16      zwrk(M2+L2+M2+L2)
#endif  /* DCOMPLEX */
      BIGCOMPLEX      z(M2+L2+M2+L2), z_old(M2+L2+M2+L2)
      integer         i, j, ikb, ierr, it
      BIGREAL         xb(M2+L2+M2+L2), yb(M2+L2+M2+L2), wrk(0:M2+L2), &
     &                wxi(0:L2), weta(0:M2), s1(0:L2), s2(0:M2), &
     &                xint(0:M2+L2), yint(0:M2+L2), &
     &                x(0:L2big,0:M2big), y(0:L2big,0:M2big), &
     &                rhs(0:L2big,0:M2big), ewrk(nwrk), huge
      parameter  (    huge = 1.e+35   )
      BIGREAL         pertrb, dum, t1, t2, zero
      parameter  (    zero = 0.d0  )
      character*40    lbl
      BIGREAL         errplt(Lm,Mm), stmp, error, area, tarea, error_old
      equivalence     (rhs,errplt)
      external        coef, bndyc, cofx, cofy

! label for grid
#include "gridid.h"
      write(lbl,100) gridid(1:40)
  100 format(40a)

#if PLOTS
      call start_plot
#endif  /* PLOTS */

!  original distribution of x,y points is preserved on boundary kb1
!  and kb2:
      if(kb1.ne.1.and.kb1.ne.3) &
     &        call crash('boundary index kb1 must be 1 or 3', kb1)
      if(kb2.ne.2.and.kb2.ne.4) &
     &        call crash('boundary index kb2 must be 2 or 4', kb2)

!  set up boundary spline interpolation arrays
      call readbndy

!  initialize vector z (complex) with contour of physical boundary
      call z_init(z, xb, yb)

!  Error code
      ierr = 0

!  Calculate departure of contour from rectangular
      error = 0.
#if DBLEPREC
      do i=1,N1
         error = error + dabs(dreal(z(i))-dreal(z(1)))
      enddo
      do i=N1+1,N2
         error = error + dabs(dimag(z(i))-dimag(z(N1+1)))
      enddo
      do i=N2+1,N3
         error = error + dabs(dreal(z(i))-dreal(z(N2+1)))
      enddo
      do i=N3+1,N4
         error = error + dabs(dimag(z(i))-dimag(z(N3+1)))
      enddo
#else
      do i=1,N1
         error = error + abs(real(z(i))-real(z(1)))
      enddo
      do i=N1+1,N2
         error = error + abs(aimag(z(i))-aimag(z(N1+1)))
      enddo
      do i=N2+1,N3
         error = error + abs(real(z(i))-real(z(N2+1)))
      enddo
      do i=N3+1,N4
         error = error + abs(aimag(z(i))-aimag(z(N3+1)))
      enddo
#endif  /* DBLEPREC */
      error = error/FLoaT(N4)
      error_old = 1.2*error

      it = 0
!  Map physical boundary to a rectangle
      do while (ierr==0 .and. error < error_old)
         it = it+1
         error_old = error
         z_old = z
#if DCOMPLEX
         call rect(z,N,N1,N2,N3,N4,ierr)
#else
         call rects(z,N,N1,N2,N3,N4,ierr)
#endif  /* DCOMPLEX */

!  Calculate departure of contour from rectangular
         error = 0.
#if DBLEPREC
         do i=1,N1
            error = error + dabs(dreal(z(i))-dreal(z(1)))
         enddo
         do i=N1+1,N2
            error = error + dabs(dimag(z(i))-dimag(z(N1+1)))
         enddo
         do i=N2+1,N3
            error = error + dabs(dreal(z(i))-dreal(z(N2+1)))
         enddo
         do i=N3+1,N4
            error = error + dabs(dimag(z(i))-dimag(z(N3+1)))
         enddo
#else
         do i=1,N1
            error = error + abs(real(z(i))-real(z(1)))
         enddo
         do i=N1+1,N2
            error = error + abs(aimag(z(i))-aimag(z(N1+1)))
         enddo
         do i=N2+1,N3
            error = error + abs(real(z(i))-real(z(N2+1)))
         enddo
         do i=N3+1,N4
            error = error + abs(aimag(z(i))-aimag(z(N3+1)))
         enddo
#endif  /* DBLEPREC */
         error = error/FLoaT(N4)
         write(6,110) it,error
  110    format(' rectangularity error in mapped contour at', &
     &        ' iteration ',i2,' is ',1pe11.4)
      enddo
      z = z_old

!  Store distribution of xi,eta points along boundaries (used to
!  compute coefficients of elliptic equation)
      do i=0,L2
#if DBLEPREC
         if (kb2 .eq. 2) sxi(i) = dreal(z(N1+i))
         if (kb2 .eq. 4) sxi(i) = dreal(z(N4-i))
#else
         if (kb2 .eq. 2) sxi(i) = real(z(N1+i))
         if (kb2 .eq. 4) sxi(i) = real(z(N4-i))
#endif  /* DBLEPREC */
      enddo

      do j=0,M2-1
#if DBLEPREC
         if (kb1 .eq. 1) seta(j) = dimag(z(N1-j))
         if (kb1 .eq. 3) seta(j) = dimag(z(N2+j))
#else
         if (kb1 .eq. 1) seta(j) = aimag(z(N1-j))
         if (kb1 .eq. 3) seta(j) = aimag(z(N2+j))
#endif  /* DBLEPREC */
      enddo
#if DBLEPREC
      seta(M2) = dimag(z(N4))
#else
      seta(M2) = aimag(z(N4))
#endif  /* DBLEPREC */

!  Cubic spline interpolation of mapping on boundaries 3 and 4 to
!  match distribution of points with those on boundaries 1 and 2
!  This also does the interpolation when Lsmall and L are different.
!  It actually does too much work when L==Lsmall.

!  boundary 1
      do j=0,M2-1
         ikb = N1-j
#if DBLEPREC
         weta(j) = dimag(z(ikb))
#else
         weta(j) = aimag(z(ikb))
#endif  /* DBLEPREC */
         xint(j) = xb(ikb)
         yint(j) = yb(ikb)
      enddo
      ikb = N4
#if DBLEPREC
      weta(N1) = dimag(z(ikb))
#else
      weta(N1) = aimag(z(ikb))
#endif  /* DBLEPREC */
      xint(N1) = xb(ikb)
      yint(N1) = yb(ikb)
      call spline(weta,xint,M2+1,huge,huge,wrk)
      do j=0,M2big
         call splint(weta,xint,wrk,M2+1,seta(j),x(0,j))
      enddo
      call spline(weta,yint,M2+1,huge,huge,wrk)
      do j=0,M2big
         call splint(weta,yint,wrk,M2+1,seta(j),y(0,j))
      enddo

!  boundary 3
      do j=0,M2
         ikb = N2+j
#if DBLEPREC
         weta(j) = dimag(z(ikb))
#else
         weta(j) = aimag(z(ikb))
#endif  /* DBLEPREC */
         xint(j) = xb(ikb)
         yint(j) = yb(ikb)
      enddo
      call spline(weta,xint,M2+1,huge,huge,wrk)
      do j=0,M2big
         call splint(weta,xint,wrk,M2+1,seta(j),x(L2big,j))
      enddo
      call spline(weta,yint,M2+1,huge,huge,wrk)
      do j=0,M2big
         call splint(weta,yint,wrk,M2+1,seta(j),y(L2big,j))
      enddo

!  boundary 2
      do i=0,L2
         ikb = N1+i
#if DBLEPREC
         wxi(i)  = dreal(z(ikb))
#else
         wxi(i)  = real(z(ikb))
#endif  /* DBLEPREC */
         xint(i) = xb(ikb)
         yint(i) = yb(ikb)
      enddo
      call spline(wxi,xint,L2+1,huge,huge,wrk)
      do i=0,L2big
         call splint(wxi,xint,wrk,L2+1,sxi(i),x(i,0))
      enddo
      call spline(wxi,yint,L2+1,huge,huge,wrk)
      do i=0,L2big
         call splint(wxi,yint,wrk,L2+1,sxi(i),y(i,0))
      enddo

!  boundary 4
      do i=0,L2
         ikb = N4-i
#if DBLEPREC
         wxi(i)  = dreal(z(ikb))
#else
         wxi(i)  = real(z(ikb))
#endif  /* DBLEPREC */
         xint(i) = xb(ikb)
         yint(i) = yb(ikb)
      enddo
      call spline(wxi,xint,L2+1,huge,huge,wrk)
      do i=0,L2big
         call splint(wxi,xint,wrk,L2+1,sxi(i),x(i,M2big))
      enddo
      call spline(wxi,yint,L2+1,huge,huge,wrk)
      do i=0,L2big
         call splint(wxi,yint,wrk,L2+1,sxi(i),y(i,M2big))
      enddo

!  Set right hand side of elliptic equation to zero
      do j=1,M2big-1
         do i=1,L2big-1
            x(i,j) = 0.
            y(i,j) = 0.
            rhs(i,j) = 0.
         enddo
      enddo

      ewrk(1) = nwrk
      t1 = L2big
      t2 = M2big
      call sepeli(0,2,zero,t1,L2big,1,dum,dum,dum,dum, &
     &                   zero,t2,M2big,1,dum,dum,dum,dum, &
     &                   cofx,cofy,rhs,x,L2big+1,ewrk,pertrb,ierr)
      if (ierr .gt. 0) then
         call crash('sepeli failed while computing x solution', ierr)
      end if
      ewrk(1) = nwrk
      call sepeli(0,2,zero,t1,L2big,1,dum,dum,dum,dum, &
     &                   zero,t2,M2big,1,dum,dum,dum,dum, &
     &                   cofx,cofy,rhs,y,L2big+1,ewrk,pertrb,ierr)
      if (ierr .gt. 0) then
         call crash('sepeli failed while computing y solution', ierr)
      end if

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

      if (pleft1) then
         write (iout1,*) nb1pts
         do j=1,nb1pts
            write (iout1,*) x1spl(j),y1spl(j)
         enddo
         write (iout1,*) Lmiddle
         do i=1,Lmiddle
            write (iout1,*) xp(i,1),yp(i,1)
         enddo
         write (iout1,*) M
         do j=1,M
            write (iout1,*) xp(Lmiddle,j),yp(Lmiddle,j)
         enddo
         write (iout1,*) Lmiddle
         do i=Lmiddle,1,-1
            write (iout1,*) xp(i,M),yp(i,M)
         enddo
      end if
      if (pleft2) then
         write (iout2,*) M
         do j=M,1,-1
            write (iout2,*) xp(1,j),yp(1,j)
         enddo
      end if

      if (pbot1) then
         write (iout1,*) Mmiddle
         do j=Mmiddle,1,-1
            write (iout1,*) xp(1,j),yp(1,j)
         enddo
         write (iout1,*) nb2pts
         do i=1,nb2pts
            write (iout1,*) x2spl(i),y2spl(i)
         enddo
         write (iout1,*) Mmiddle
         do j=1,Mmiddle
            write (iout1,*) xp(L,j),yp(L,j)
         enddo
         write (iout1,*) L
         do i=L,1,-1
            write (iout1,*) xp(i,Mmiddle),yp(i,Mmiddle)
         enddo
      end if
      if (pbot2) then
         write (iout2,*) L
         do i=1,L
            write (iout2,*) xp(i,1),yp(i,1)
         enddo
      end if

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
      write(6,130)area
  130 format(5x,'area         ',g11.4)

!  Check orthogonality by evaluating dx/dxi*dx/deta+dy/dxi*dy/deta
!  everywhere.  Normalize with respect to local grid cell area.
!  Store result in errplt and pass to gridplot.
      do j=1,Mm
         do i=1,Lm
            errplt(i,j) = ((x(2*i,2*j-1)-x(2*i-2,2*j-1))* &
     &        (x(2*i-1,2*j)-x(2*i-1,2*j-2))+ &
     &        (y(2*i,2*j-1)-y(2*i-2,2*j-1))* &
     &        (y(2*i-1,2*j)-y(2*i-1,2*j-2)))* &
     &        pm(i,j)*pn(i,j)
         enddo
      enddo

!  Output solution
!  (This output format parallels the input format assumed by
!  subroutine getgrid in pemodel)
      call checkdefs
      call def_grid
#if PLOTS
      call gridplot(lbl,area,tarea,errplt)
      call end_plot
#endif  /* PLOTS */

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
      end

! *******************************************************************

      subroutine bndyc (kbdy,xory,alfa,gbdy)
      integer    kbdy
      BIGREAL    xory, alfa, gbdy
      write (6,100)
  100 format (' you are in bndyc: tut tut.')
      stop
      end

! *******************************************************************

      subroutine coef (x,y,cxx,cyy,cx,cy,c)
      BIGREAL    x, y, cxx, cyy, cx, cy, c
#include "grid.h"
      integer    i, j
      BIGREAL    dxdi, dedj

      i = nint(x)
      j = nint(y)
      dxdi = 0.5*(sxi(i+1)-sxi(i-1))
      dedj = 0.5*(seta(j+1)-seta(j-1))
      cxx = 1./dxdi**2
      cyy = 1./dedj**2
      cx  = (-sxi(i+1)+2.*sxi(i)-sxi(i-1))/dxdi**3
      cy  = (-seta(j+1)+2.*seta(j)-seta(j-1))/dedj**3
      c   = 0.
      return
      end

! *******************************************************************

      subroutine cofx(xx,afun,bfun,cfun)
      BIGREAL    xx, afun, bfun, cfun

!  subroutine to compute the coefficients of the elliptic equation
!  solved to fill in the grid.  it is passed (along with cofy) to
!  the elliptic solver sepeli.

#include "grid.h"
      integer    i
      BIGREAL    dxdi
      i = int(xx)
      dxdi = 0.5 * (sxi(i+1)-sxi(i-1))
      afun = 1. / dxdi ** 2
      bfun = (-sxi(i+1)+2.*sxi(i)-sxi(i-1)) / dxdi ** 3
      cfun = 0.
      return
      end

! *******************************************************************

      subroutine cofy(yy,dfun,efun,ffun)
      BIGREAL    yy, dfun, efun, ffun
#include "grid.h"
      integer    j
      BIGREAL    dedj
      j = int(yy)
      dedj = 0.5 * (seta(j+1)-seta(j-1))
      dfun = 1. / dedj ** 2
      efun = (-seta(j+1)+2.*seta(j)-seta(j-1)) / dedj ** 3
      ffun = 0.
      return
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
      subroutine gridplot(lbl,area,tarea,errplt)
#include "grid.h"
      character*40    lbl
      character*5     lbls
      BIGREAL         area, tarea
      BIGREAL         errplt(Lm,Mm)
      SMALLREAL       tmp(0:L,0:M), xlow, xhigh, ylow, yhigh, &
     &                x1, x2, y1, y2, dxmin, tmp2(Lm,Mm)
      integer         i, j
      BIGREAL         vmin, vmax

      call pcseti('QUALITY',1)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,lbl,.012,0.,0.)
      write(lbl,100)area
  100 format(' Sum 1/MN=',g11.4)
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
      write(lbls,110)
  110 format('DX')
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
      call set(.05,.95,.1,.9,1.,float(Lp),1.,float(Mp),1)
      call perim(1,L,1,M)
      call cpseti('SET',0)
      call cpcnrc(tmp,Lp,Lp,Mp,0.0,0.0,0.0,1,-1,-682)
      call frame
      write(lbls,120)
  120 format('DY')
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
      call set(.05,.95,.1,.9,1.,float(Lp),1.,float(Mp),1)
      call perim(1,L,1,M)
      call cpseti('SET',0)
      call cpcnrc(tmp,Lp,Lp,Mp,0.0,0.0,0.0,1,-1,-682)
      call frame

!  plot error in grid
      do j=1,Mm
        do i=1,Lm
          tmp2(i,j) = errplt(i,j)
        enddo
      enddo
      write(lbls,130)
  130 format('Error')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(.5,.98,lbls,.012,0.,0.)
      call set(.05,.95,.1,.9,1.,float(Lm),1.,float(Mm),1)
      call perim(1,Lm-1,1,Mm-1)
      call cpseti('SET',0)
      call cpcnrc(tmp2,Lm,Lm,Mm,0.0,0.0,0.0,1,-1,-682)
      call frame
      return
      end
#endif  /* PLOTS */

! *******************************************************************

      subroutine readbndy
#include "grid.h"
      BIGREAL    arcs(IBIG), huge
      parameter (huge = 1.e+35   )
      integer    i

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
      return
      end

! *******************************************************************

#if DCOMPLEX
      subroutine rect(z,N,N1,N2,N3,N4,err)
      integer    N, N1, N2, N3, N4, err

!     This subroutine is taken directly from Ives, D.C. and
!     R.M. Zacharias "Conformal Mapping and Orthogonal Grid Generation"
!     AIAA-87-2057.

      complex*16 z(N), z0, zd
      integer    BIG
      parameter (BIG=8000)
      double precision r(BIG), t(BIG), pi, alpha, pwr, pmin, &
     &           pmax, tp
      integer    i, j, k, im, ip

!  Check on work array size
      if (N .gt. BIG) then
         call crash('BIG in rect must be at least ', N)
      endif
      pi = 4.d0*datan(1.d0)
      do 160 i=1,N
         im = N-mod(N-i+1,N)
         ip = 1+mod(i,N)
         zd = (z(im)-z(i))/(z(ip)-z(i))
         alpha = datan2(dimag(zd),dreal(zd))
         if(alpha.lt.0)alpha = alpha+pi+pi
         pwr = pi/alpha
         if(i.ne.N1.and.i.ne.N2.and.i.ne.N3.and.i.ne.N4) goto 110
         zd = (z(im)-z(i))/cdabs(z(im)-z(i))
         do 100 j=1,N
            z(j) = dcmplx(0.d0,1.d0)*z(j)/zd
  100    continue
         pwr = pwr/2.
  110    pmin = 100.
         pmax =-100.
         tp   = 0.
         do 140 j=2,N
            zd = z(mod(j+i-2,N)+1)-z(i)
            r(j) = cdabs(zd)
            t(j) = datan2(dimag(zd),dreal(zd))-pi-pi-pi-pi-pi-pi
            do 120 k=1,7
               if(dabs(t(j)-tp).lt.pi)goto 130
               t(j) = t(j)+pi+pi
  120       continue
!        pause ' warning - tracking error '
            err = 1
            return
  130       tp = t(j)
            pmax = dmax1(pmax,t(j)*pwr)
            pmin = dmin1(pmin,t(j)*pwr)
  140    continue
         pwr = dmin1(pwr,1.98d0*pi*pwr/(pmax-pmin))
         z(i) = dcmplx(0.d0,0.d0)
         do 150 j=2,N
#if __sgi
            z(mod(j+i-2,N)+1) = r(j)**pwr*zexp(dcmplx(0.d0,t(j)*pwr))
#else
            z(mod(j+i-2,N)+1) = r(j)**pwr*cdexp(dcmplx(0.d0,t(j)*pwr))
#endif
  150    continue
         zd = 1./(z(N2)-z(N1))
         z0 = z(N1)
      do 160 j=1,N
         z(j) = (z(j)-z0)*zd
  160 continue
      return
      end

#else
! *******************************************************************

      subroutine rects(z,N,N1,N2,N3,N4,err)
      integer    N, N1, N2, N3, N4, err

!     Single precision version of subroutine rect from Ives,D.C. and
!     R.M.Zacharias "Conformal Mapping and Orthogonal Grid Generation"
!     AIAA-87-2057.
      complex    z(N), z0, zd
      integer    BIG
      parameter (BIG=3000)
      real       r(BIG), t(BIG), pi, alpha, pwr, pmin, &
     &           pmax, tp
      integer    i, j, k, im, ip

!  Check on work array size
      if (N .gt. BIG) then
         call crash('BIG in rect must be at least ', N)
      endif
      pi = 4.*atan(1.)
      do 160 i=1,N
         im = N-mod(N-i+1,N)
         ip = 1+mod(i,N)
         zd = (z(im)-z(i))/(z(ip)-z(i))
         alpha = atan2(aimag(zd),real(zd))
         if(alpha.lt.0)alpha = alpha+pi+pi
         pwr = pi/alpha
         if(i.ne.N1.and.i.ne.N2.and.i.ne.N3.and.i.ne.N4)goto 110
         zd = (z(im)-z(i))/cabs(z(im)-z(i))
         do 100 j=1,N
            z(j) = cmplx(0.,1.)*z(j)/zd
  100    continue
         pwr = pwr/2.
  110    pmin = 100.
         pmax =-100.
         tp   = 0.
         do 140 j=2,N
            zd = z(mod(j+i-2,N)+1)-z(i)
            r(j) = cabs(zd)
            t(j) = atan2(aimag(zd),real(zd))-pi-pi-pi-pi-pi-pi
            do 120 k=1,7
               if(abs(t(j)-tp).lt.pi)goto 130
               t(j) = t(j)+pi+pi
  120       continue
!        pause ' warning - tracking error '
            err = 1
            return
  130       tp = t(j)
            pmax = max(pmax,t(j)*pwr)
            pmin = min(pmin,t(j)*pwr)
  140    continue
         pwr = min(pwr,1.98*pi*pwr/(pmax-pmin))
         z(i) = cmplx(0.,0.)
         do 150 j=2,N
            z(mod(j+i-2,N)+1) = r(j)**pwr*cexp(cmplx(0.,t(j)*pwr))
  150    continue
         zd = 1./(z(N2)-z(N1))
         z0 = z(N1)
      do 160 j=1,N
         z(j) = (z(j)-z0)*zd
  160 continue
      return
      end
#endif  /* DCOMPLEX */

! *******************************************************************

      subroutine spline(x,y,NN,yp1,ypn,y2)
      integer    NN
      BIGREAL    x(NN), y(NN), yp1, ypn, y2(NN)

!     The following two subroutines are used to perform the cubic spline
!     interpolation required to match up the distribution of points on
!     opposite sides of the transformed plane rectangle.  The routines
!     are taken from  Press,W.H., B.P.Flannery, S.A.Teukolsky and
!     W.T.Vetterling: "Numerical Recipes, the Art of Scientific Computing"
!     Cambridge University Press, 1986.
#include "grid.h"
      integer         nmax
      parameter (     nmax=M2+L2             )
      BIGREAL    u(nmax), qn, un, sig, p
      integer    i, k

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

      subroutine splint(xa,ya,y2a,N,x,y)
      integer    N
      BIGREAL    xa(N), ya(N), y2a(N), x, y
      integer    klo, khi, k
      BIGREAL    h, a, b

      klo=1
      khi=N
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

! *******************************************************************

      subroutine z_init(z, xb, yb)
#include "grid.h"
      BIGCOMPLEX      z(M2+L2+M2+L2)
      BIGREAL         xb(M2+L2+M2+L2), yb(M2+L2+M2+L2)
      BIGREAL         stmp
      integer         i

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
#if XPOTS1
      close(ipot1)
#endif  /* XPOTS1 */
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
#if XPOTS2
      close(ipot2)
#endif  /* XPOTS2 */
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
#if XPOTS3
      close(ipot3)
#endif  /* XPOTS3 */
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
#if XPOTS4
      close(ipot4)
#endif  /* XPOTS4 */
      return
      end
