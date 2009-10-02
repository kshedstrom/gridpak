      program maptst

c  *** In gridpak version 5.4 ***** October 18, 2001 ****************
c  Kate Hedstrom  (kate@arsc.edu)
c  John Wilkin  (wilkin@imcs.rutgers.edu)
c  ******************************************************************
c
c  This is a program to convert from the coordinate system in coast.f
c  to back to latitude and longitude.  If you change the parameters
c  in coast.f you have to change them here as well, including HINIT
c  and DH.  This program writes out the four lat,long grids on unit 62.
c  ******************************************************************

#include "proj.h"
#include "gridparam.h"
      parameter (     REarth = 6.3708e6 )
      character*40    ident
      character*80    gridid

      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

      call opngks

c  Initialize the mapping variables

      call pcseti('QUALITY',1)
#include "gridid.h"
      write(ident,100) gridid(1:40)
  100 format(40a)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ident,.012,0.,0.)
      call mapsti('GR',JGRD)
      call mapstc('OU','PS')
      call maproj(JPRJ,PLAT,PLONG,ROTA)
      call mapset(JLTS,P1,P2,P3,P4)
      call mapdrw
      call maptrn(PLAT+.5,PLONG,U2,V2)
      call maptrn(PLAT-.5,PLONG,U1,V1)
      udeg = sqrt((U2-U1)*(U2-U1) + (V2-V1)*(V2-V1))

c  Calculate the lat,long of the psi points

      hlat = HINIT
      do i=1,100
        print *,'Lat: '
	read *,rmlat
	if (rmlat .gt. 90) goto 120
	print *,'Long: '
	read *,rmlon
	call maptrn(rmlat,rmlon,u,v)
	print *, 'u,v ',u,v
        call mapinv(u,v,psilat,psilon,hlat-DH,hlat+DH)
	print *, 'lat, long ', psilat, psilon
      enddo
  120 continue
      call frame
      call clsgks

      call exit(0)
  250 call crash('read or write error while outputing solution', 0)
      call exit(1)
      end

C----------------------------------------------------------------------

      subroutine crash(icrash,ierr)
      character*(*) icrash
      print *,icrash
      if (ierr .ne. 0) print *,'ierr = ',ierr
      call exit(1)
      return
      end

C----------------------------------------------------------------------

      subroutine fx(A,f,dfdx)
      common /MAPKH1/ r
      save /MAPKH1/

      if (A .LE. .0001) then
        f = sin(A)/2. - r
        dfdx = cos(A)/2.
      else
        f = (1-cos(A))/sin(A) - r
        dfdx = (1-cos(A))/(sin(A)*sin(A))
      end if
      return
      end

C----------------------------------------------------------------------

      subroutine funcd1(lat,f,dfdx)
      real lat
      COMMON /MAPCM1/ IPRJ,SINO,COSO,SINR,COSR,PHOC
      SAVE /MAPCM1/
      common /MAPKH2/ A,B
      save /MAPKH2/

      f = sin(lat)*coso - sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat)))*cos(lat)*sino - cos(B)*sin(A)
      dfdx = cos(lat)*coso + sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat)))*sin(lat)*sino
     &   + sin(lat)*sino/sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat))) *(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat))
      return
      end

C----------------------------------------------------------------------

      subroutine funcd2(lat,f,dfdx)
      real lat
      COMMON /MAPCM1/ IPRJ,SINO,COSO,SINR,COSR,PHOC
      SAVE /MAPCM1/
      common /MAPKH2/ A,B
      save /MAPKH2/

      f = sin(lat)*sino + sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat)))*cos(lat)*coso - cos(A)
      dfdx = cos(lat)*sino - sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat)))*sin(lat)*coso
     &   - sin(lat)*coso/sqrt(1-(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat))) *(sin(A)*sin(A)*sin(B)*sin(B)
     &       /cos(lat)/cos(lat))
      return
      end

C----------------------------------------------------------------------

      subroutine gx(A,f,dfdx)
      common /MAPKH1/ r
      save /MAPKH1/

      f = (1+cos(A))/sin(A) - r
      dfdx = (1+cos(A))/(sin(A)*sin(A))
      return
      end

C----------------------------------------------------------------------

      subroutine mapinv (u,v,rlat,rlon,MINLAT,MAXLAT)
      real MINLAT,MAXLAT

c Inverse of MAPTRN.  Only works for some map projections -
c IPRJ = 1,2,4,5,6,11.
c IPRJ corresponds to the projections as follows:

C       IPRJ	PROJECTION                      TYPE
c	 1	Lambert conformal conic         conic
c	 2	Stereographic                   azimuthal
c	 3	Orthographic or Satellite-view
c	 4	Lambert equal area
c	 5	Gnomonic
c	 6	Azimuthal equidistant
c	 7	Cylindrical equidistant         cylindrical
c	 8	Mercator
c	 9	Mollweide
c	10	Fast-path cylindrical equidistant
c	11	Fast-path Mercator
c	12	Fast-path Mollweide

C Declare required common blocks.  See MAPBD for descriptions of these
C common blocks and the variables in them.

      parameter (pi=3.1415926536)
      common /MAPCM1/ IPRJ,SINO,COSO,SINR,COSR,PHOC
      save /MAPCM1/
      common /MAPCMB/ IIER
      save /MAPCMB/
      common /MAPCM4/ INTF,JPRJ,PHIA,PHIO,ROTA,ILTS,PLA1,PLA2,PLA3,PLA4,
     &                PLB1,PLB2,PLB3,PLB4,PLTR,GRID,IDSH,IDOT,LBLF,PRMF,
     &                ELPF,XLOW,XROW,YBOW,YTOW,IDTL,GRDR,SRCH,ILCW
      logical         INTF,LBLF,PRMF,ELPF
      save /MAPCM4/
      common /MAPKH1/ r
      save /MAPKH1/
      common /MAPKH2/ A,B
      save /MAPKH2/
      external fx,gx,funcd1,funcd2
      real lat

C Define required constants.  DTOR is pi over 180, DTRH is half of DTOR
C or pi over 360, and TOPI is 2 over pi.

      data DTOR / .017453292519943 /
      data DTRH / .008726646259971 /
      data RTOD / 57.2957795130823 /
      data TOPI / .636619772367581 /
      data PIOT / 1.570796327 /

      if (iprj .eq. 11) then
c  Fast-path Mercator.

        print *, 'PHOC ', PHOC, PHIO, PHIA
        ph = u + PHIO*DTOR
        lat = -PIOT + 2*atan(exp(v))

      else if (iprj .eq. 1) then
c  Lambert Conformal Conic

        tlon = atan(-sino*u/v)
        r = sqrt(v*v + u*u)
        rlon = RTOD*tlon/coso
        chi = 2*RTOD*atan(r**(1/coso))
        rlat = (90-chi)/sino
c  Clip to PHOC +- 180 degrees
        rlon = rlon-sign(180.,rlon+180.)+sign(180.,180-rlon) + PHOC
        return

      else
c Azimuthal projections.

c  Find the radius and angle of the point.

        r = sqrt(u*u + v*v)
        B = atan2(u,v) - ROTA*DTOR

        if (iprj .eq. 2) then
C Stereographic.

          if (r .lt. .0001) then
            A = asin(2.*r)
          else
            A = rtsafe(fx,.00001,pi-.00001,.00001)
          end if

        else if (iprj .eq. 4) then
C Lambert equal area.

          r = sqrt(4/(r*r) - 1)
          A = rtsafe(gx,.0001,pi-.00001,.00001)
          if (abs(COS(A)+1.).lt.1.e-6) go to 100

        else if (iprj .eq. 5) then
C Gnomonic.

          A = atan(R)

        else if (iprj .eq. 6) then
C Azimuthal equidistant.

          A=R

        else
C Projection not allowed
          goto 110
        end if

C rlon is a longitude, in degrees, between -180. and +180.,
C inclusive, and rlat is a latitude, in degrees.  ph and lat are
C phase (lon) and latitude in radians.

c Look for North Pole.  Change if you need both poles (sorry).

c      if ((abs(B) .le. .0001) .and. abs(SINO - cos(A)) .le .0001) then
c          lat = 1./TOPI
c          ph = 0.
c       else
        lat = rtsafe(funcd1,MINLAT*DTOR,MAXLAT*DTOR,.00001)
        if (lat .eq. 1.e20) then
          lat = rtsafe(funcd2,MINLAT*DTOR,MAXLAT*DTOR,.00001)
          if (lat .eq. 1.e20)
     &      call crash('We have a problem in mapinv', 0)
        end if
        ph = asin(sin(B)*sin(A)/cos(lat))
c       end if

c  End of Spherical projections
      end if

      rlat = RTOD*lat
      tmp = RTOD*ph
c  Clip to PHOC +- 180 degrees
      rlon = tmp - sign(180.,tmp+180.) + sign(180.,180.-tmp) + PHOC

      return

C Projection of point is invisible or undefined.

  100 rlat = 90.
      rlon = 361.
      return

C Error exit.

  110 if (IIER.ne.0) go to 100
      IIER=16
      call seter (' mapinv - attempt to use non-existent projection',
     &              IIER,1)
      go to 100

      end

C----------------------------------------------------------------------

      FUNCTION RTSAFE(FUNCD,X1,X2,XACC)
      PARAMETER (MAXIT=100)
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
c     IF (FL*FH .GE. 0.) PAUSE 'root must be bracketed'
      IF (FL*FH .GE. 0.) then
c  bogus value - flag that there was a problem
        rtsafe = 1.e20
        return
      end if
      IF (FL .LT. 0.) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      ENDIF
      RTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(RTSAFE,F,DF)
      DO 100 J=1,MAXIT
        IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.
     &      .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          RTSAFE=XL+DX
          IF(XL.EQ.RTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=RTSAFE
          RTSAFE=RTSAFE-DX
          IF(TEMP.EQ.RTSAFE)RETURN
        ENDIF
        IF(ABS(DX).LT.XACC) RETURN
        CALL FUNCD(RTSAFE,F,DF)
        IF(F.LT.0.) THEN
          XL=RTSAFE
          FL=F
        ELSE
          XH=RTSAFE
          FH=F
        ENDIF
  100 CONTINUE
      PAUSE 'RTSAFE exceeding maximum iterations'
      RETURN
      END

C----------------------------------------------------------------------
