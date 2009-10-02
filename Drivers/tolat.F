      program tolat

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!  This is a program to convert from the coordinate system in coast.f
!  to back to latitude and longitude.  If you change the parameters
!  in coast.f you have to change them here as well, including HINIT
!  and DH.  This program writes out the four lat,long grids to the
!  grid netCDF file.
!  ******************************************************************

#include "griddefs.h"
#include "proj.h"
#include "bathy.h"
#include "ncgrid.h"
      BIGREAL         REarth
      parameter (     REarth = 6.3708e6 )
      integer         i, j
      BIGREAL         minlat, minlon, maxlat, maxlon
      character*40    ident
      logical         dateline
      BIGREAL         DTOR, RTOD
      BIGREAL         udeg, vmin, vmax, addphase
      real            U1, U2, V1, V2, u, v, lat, lon

      data DTOR / .017453292519943 /
      data RTOD / 57.2957795130823 /

#include "gridid.h"

      call get_xy
#ifndef PLOTS
      call crash('This program depends on NCAR graphics', 1)
#endif  /* PLOTS */
      call start_plot

!  Initialize the mapping variables

      call pcseti('QUALITY',1)
      write(ident,100) gridid(1:40)
  100 format(40a)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,ident,.012,0.,0.)
      call mapsti('GR',JGRD)
      call mapstc('OU','NO')
      call maproj(JPRJ,PLAT,PLONG,ROTA)
      call mapset(JLTS,P1,P2,P3,P4)
      call mapdrw
      call plotit(0,0,0)
      call drawcoast
      call mapiq
      call maptrn(PLAT+.5,PLONG,U2,V2)
      call maptrn(PLAT-.5,PLONG,U1,V1)
      udeg = sqrt((U2-U1)*(U2-U1) + (V2-V1)*(V2-V1))

!  Calculate the lat,long of the psi points

      do i=1,L
        do j=1,M
          u = (xp(i,j)-XOFF)*udeg*RTOD/REarth
          v = (yp(i,j)-YOFF)*udeg*RTOD/REarth
          call maptri(u,v,lat,lon)
          lat_psi(i,j) = lat
          lon_psi(i,j) = lon
        enddo
      enddo

!  Plot the grid in lat,long space

#if DBLEPREC
      do j = 1,M
        call mapit(sngl(lat_psi(1,j)),sngl(lon_psi(1,j)),0)
        do i = 2,L
          call mapit(sngl(lat_psi(i,j)),sngl(lon_psi(i,j)),2)
        enddo
      enddo
      do i = 1,L
        call mapit(sngl(lat_psi(i,1)),sngl(lon_psi(i,1)),0)
        do j = 2,M
          call mapit(sngl(lat_psi(i,j)),sngl(lon_psi(i,j)),2)
        enddo
      enddo
#else
      do j = 1,M
        call mapit(lat_psi(1,j),lon_psi(1,j),0)
        do i = 2,L
          call mapit(lat_psi(i,j),lon_psi(i,j),2)
        enddo
      enddo
      do i = 1,L
        call mapit(lat_psi(i,1),lon_psi(i,1),0)
        do j = 2,M
          call mapit(lat_psi(i,j),lon_psi(i,j),2)
        enddo
      enddo
#endif  /* DBLEPREC */
      call mapiq
      call frame
      call end_plot

!  Calculate the lat,long of the u points

      do i=1,L
        do j=0,M
          u = (xu(i,j)-XOFF)*udeg*RTOD/REarth
          v = (yu(i,j)-YOFF)*udeg*RTOD/REarth
          call maptri(u,v,lat,lon)
          lat_u(i,j) = lat
          lon_u(i,j) = lon
        enddo
      enddo

!  Calculate the lat,long of the v points

      do i=0,L
        do j=1,M
          u = (xv(i,j)-XOFF)*udeg*RTOD/REarth
          v = (yv(i,j)-YOFF)*udeg*RTOD/REarth
          call maptri(u,v,lat,lon)
          lat_v(i,j) = lat
          lon_v(i,j) = lon
        enddo
      enddo

!  Calculate the lat,long of the rho points

      do i=0,L
        do j=0,M
          u = (xr(i,j)-XOFF)*udeg*RTOD/REarth
          v = (yr(i,j)-YOFF)*udeg*RTOD/REarth
          call maptri(u,v,lat,lon)
          lat_rho(i,j) = lat
          lon_rho(i,j) = lon
        enddo
        print *,i,lat_rho(i,0),lon_rho(i,0),lat_rho(i,M),lon_rho(i,M)
      enddo

      maxlat = vmax(lat_rho,Lp*Mp)
      minlat = vmin(lat_rho,Lp*Mp)
      maxlon = vmax(lon_rho,Lp*Mp)
      minlon = vmin(lon_rho,Lp*Mp)
      print *,'Latitude range ', minlat, maxlat
      print *,'Longitude range ', minlon, maxlon
      if (maxlat .eq. 1.e12 .or. maxlon .eq. 1.e12) then
        call crash('need to modify proj.h plot range', 1)
      end if

!  Write out the lat,long information

!  First check that grid doesn't span the international dateline
!  If so this can mess up things in topo and sphere
!  Beware: This test might not trap some nasty cases

      dateline = .false.

      addphase = 0.0
      do j = 1,M
        do i = 2,L
          if (abs(lon_psi(i,j)-lon_psi(i-1,j)) .gt. 270.) then
            if (lon_psi(i,j) .lt. lon_psi(i-1,j)) then
              addphase = 360.0
            else
              addphase = -360.0
            endif
            dateline  = .true.
          else
            addphase = 0.0
          endif
          lon_psi(i,j) = lon_psi(i,j) + addphase
        enddo
      enddo

      addphase = 0.0
      do j = 0,M
        do i = 1,L
          if (abs(lon_rho(i,j)-lon_rho(i-1,j)) .gt. 270.) then
            if (lon_rho(i,j) .lt. lon_rho(i-1,j)) then
              addphase = 360.0
            else
              addphase = -360.0
            endif
            dateline  = .true.
          else
            addphase = 0.0
          endif
          lon_rho(i,j) = lon_rho(i,j) + addphase
        enddo
      enddo

      addphase = 0.0
      do j = 0,M
        do i = 2,L
          if (abs(lon_u(i,j)-lon_u(i-1,j)) .gt. 270.) then
            if (lon_u(i,j) .lt. lon_u(i-1,j)) then
              addphase = 360.0
            else
              addphase = -360.0
            endif
            dateline  = .true.
          else
            addphase = 0.0
          endif
          lon_u(i,j) = lon_u(i,j) + addphase
        enddo
      enddo

      addphase = 0.0
      do j = 1,M
        do i = 1,L
          if (abs(lon_v(i,j)-lon_v(i-1,j)) .gt. 270.) then
            if (lon_v(i,j) .lt. lon_v(i-1,j)) then
              addphase = 360.0
            else
              addphase = -360.0
            endif
            dateline  = .true.
          else
            addphase = 0.0
          endif
          lon_v(i,j) = lon_v(i,j) + addphase
        enddo
      enddo
      if(dateline) then
        print *,'tolat: grid crosses Dateline (?)'
        print *,'       recovery attempted... check'
        maxlat = vmax(lat_rho,Lp*Mp)
        minlat = vmin(lat_rho,Lp*Mp)
        maxlon = vmax(lon_rho,Lp*Mp)
        minlon = vmin(lon_rho,Lp*Mp)
        print *,'Latitude range ', minlat, maxlat
        print *,'Longitude range ', minlon, maxlon
      endif

      call wrt_lat

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
  150 call crash('read or write error while outputing solution',0)
      end

!----------------------------------------------------------------------

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
