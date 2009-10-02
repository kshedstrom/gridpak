      program coast1

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!
!  This program converts a coastline given in Latitude-Longitude pairs
!  into meters on a plane for use by the grid generation program.
!  The plane projection is any one of the ten provided by the
!  NCAR Graphics routines in EZMAP.  See the EZMAP manual for a
!  description of its (many) subroutines and their parameters.
!
!  backward is a logical variable which is true if you read in the
!  data in the direction opposite to the direction it should be
!  written out for the grid program.
! *******************************************************************

#include "coast.h"
#include "proj.h"

!  Earth angular frequency and radius of sphere of equal volume
      real            Aomega, REarth
      parameter (     Aomega=7.292115e-5 , REarth=6.3708e6  )

      real            s(IBIG), u(IBIG), v(IBIG), x(IBIG), y(IBIG)
      real            xwrk(IBIG), ywrk(IBIG)
      character*21    ident
      real            DTOR, RTOD
      integer         i, nn
      real            u1, v1, u2, v2, udeg, s1

      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

      call opngks
      call pcseti('QUALITY',1)

!  The map of the coastline database

      write(ident,100)
  100 format(' EZMAP Coast Database')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ident,.012,0.,0.)
      call mapsti('GR',JGRD)
      call mapstc('OU','PS')
      call maproj(JPRJ,PLAT,PLONG,ROTA)
      call mapset(JLTS,P1,P2,P3,P4)
      call mapdrw
      call frame

!  The map of your coastline

      call cstiput(nn)
      write(ident,110)
  110 format('Digitized Coastline')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ident,.012,0.,0.)
      call mapstc('OU','NO')
      call mapset(JLTS,P1,P2,P3,P4)
      call mapdrw

      call mapit(clat(1),clong(1),0)
      do i=2,nn
        call mapit(clat(i),clong(i),2)
      enddo
      call mapiq
      call frame

      call maptrn(PLAT+.5,PLONG,u2,v2)
      call maptrn(PLAT-.5,PLONG,u1,v1)
      udeg = sqrt((u2-u1)*(u2-u1) + (v2-v1)*(v2-v1))
      do i=1,nn
        call maptrn(clat(i),clong(i),u(i),v(i))
        x(i) = u(i)*DTOR*REarth/udeg + XOFF
        y(i) = v(i)*DTOR*REarth/udeg + YOFF
      enddo

      open(unit=60,form='FORMATTED')
      write(60,*) nn
      if (backward) then
        do i=nn,1,-1
          write(60,*) x(i),' ',y(i)
        enddo
      else
        do i=1,nn
          write(60,*) x(i),' ',y(i)
        enddo
      end if
      close(60)

!  Everything up to this point is generic.  What follows should
!  use the same interpolation scheme that the grid program uses.

!  The map of your coastline, interpolated by cubic splines

      write(ident,120)
  120 format('Spline Smoothed Coast')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ident,.012,0.,0.)
      call mapset(JLTS,P1,P2,P3,P4)
      call mapdrw

      do i=1,nn
        s(i) = (i-1.)/(nn-1.)
      enddo
      call spline(s,u,nn,1.e30,1.e30,xwrk)
      call spline(s,v,nn,0.,0.,ywrk)

      call frstd(u(1),v(1))
      do i=1,400
        s1 = i/400.
        call splint(s,u,xwrk,nn,s1,u1)
        call splint(s,v,ywrk,nn,s1,v1)
        call vectd(u1,v1)
      enddo
      call mapiq
      call frame

      call clsgks
#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
      end

!*********************************************************************

      subroutine cstiput(nn)
      integer    nn
#include "coast.h"

!  I want a "while (scanf() != EOF)", but I have to fake it.
      do nn=1,IBIG
        read(5,*,err=120,end=110) clat(nn),clong(nn)
      enddo

!  If you get here then you have not reached the EOF mark
      print *,'Array size too small for the number of points'
      print *,'in coastline data file.  Only the first ',IBIG
      print *,'points will be used.'

  110 continue
      nn = nn-1
      print *,'Number of points read = ',nn
      return
  120 print *,'read error in cstiput'
#if NO_EXIT
      stop
#else
      call exit(1)
#endif  /* NO_EXIT */
      end

! **********************************************************************

      subroutine spline(x,y,N,yp1,ypn,y2)
      integer    N
      real       x(N), y(N), yp1, ypn, y2(N)

!   The following two subroutines are used to perform a cubic spline
!   interpolation.  The routines are taken from  Press,W.H., B.P.
!   Flannery, S.A.Teukolsky and W.T.Vetterling: "Numerical Recipes,
!   the Art of Scientific Computing" Cambridge University Press, 1986.

#include "coast.h"
      real       u(IBIG), sig, p, qn, un
      integer    i, k

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 100 i=2,N-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  100 continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(N)-x(N-1)))*(ypn-(y(N)-y(N-1))/(x(N)-x(N-1)))
      endif
      y2(N)=(un-qn*u(N-1))/(qn*y2(N-1)+1.)
      do 110 k=N-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
  110 continue
      return
      end

! **********************************************************************

      subroutine splint(xa,ya,y2a,N,x,y)
      integer    N
      real       xa(N), ya(N), y2a(N), x, y

      integer    klo, khi, k
      real       h, a, b

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
      if (h.eq.0.) then
        print *,'problem in splint'
#if NO_EXIT
        stop
#else
        call exit(1)
#endif  /* NO_EXIT */
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ &
     &        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
