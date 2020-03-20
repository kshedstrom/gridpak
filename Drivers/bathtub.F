      program bathtub

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!  ******************************************************************
!
!  This program reads in the bathymetric data as well as the lat,long
!  points at which depth is carried in SPEM.  It interpolates the
!  bathymetry to those grid points by a bilinear fit from the four
!  nearest neighbors.  The bathymetry data used is the ETOPO5
!  bathymetry read from the netCDF file.  It is stored at
!  5 minute intervals (1/12 degrees = DELTA_X, DELTA_Y).
!
!  ILAT is # of latitude grid points;
!  ILON is # of longitude grid points.
!  NSUB decimates the input data -- solely for plotting purposes to
!  avoid overflowing conpack work arrays.
!********************************************************************

#include "griddefs.h"
#include "bathy.h"
#include "ncgrid.h"
      integer         Lm2, Mm2
      parameter (     Lm2=L-2, Mm2=M-2   )
      integer         ILON, ILAT
      BIGREAL         DELTA_X, DELTA_Y
      logical         evenflag
#if NJB
      parameter     ( evenflag = .true. )
      parameter     ( ILAT=1319, ILON=2600 )
      parameter     ( DELTA_X=1/120.d0, DELTA_Y=1/120.d0 )
#elif ETOPO5
      parameter     ( evenflag = .true. )
      parameter     ( ILON=360*12+1, ILAT=180*12+1 )
      parameter     ( DELTA_X=1/12.d0, DELTA_Y=1/12.d0 )
#elif ETOPO2
      parameter     ( evenflag = .true. )
      parameter     ( ILON=360*30+1, ILAT=180*30+1 )
      parameter     ( DELTA_X=1/30.d0, DELTA_Y=1/30.d0 )
#elif ETOPO1
      parameter     ( evenflag = .true. )
      parameter     ( ILON=360*60+1, ILAT=180*60+1 )
      parameter     ( DELTA_X=1/60.d0, DELTA_Y=1/60.d0 )
#elif SRTM30
      parameter     ( evenflag = .true. )
      parameter     ( ILON=360*120+1, ILAT=180*120 )
      parameter     ( DELTA_X=1/120.d0, DELTA_Y=1/120.d0 )
#elif IBCAO
      parameter     ( evenflag = .true. )
      parameter     ( ILON=43201, ILAT=3121 )
      parameter     ( DELTA_X=1/120.d0, DELTA_Y=1/120.d0 )
#elif GEBCO
      parameter     ( evenflag = .true. )
      parameter     ( ILON=360*60+1, ILAT=180*60+1 )
      parameter     ( DELTA_X=1/60.d0, DELTA_Y=1/60.d0 )
#elif ARDEM
      parameter     ( evenflag = .true. )
      parameter     ( ILON=6601, ILAT=3601)
      parameter     ( DELTA_X=1/60.d0, DELTA_Y=1/120.d0)
#elif PALAU_BATHY
      parameter     ( evenflag = .false. )
      parameter     ( ILON=2850, ILAT=2409)
#else
!  Smith and Sandwell subset
      parameter     ( evenflag = .false. )
      parameter     ( ILON=721, ILAT=918)

#endif
#if IBCAO
      real, allocatable      ::       bathy(:,:)
#else
      integer*2, allocatable ::       bathy(:,:)
#if ETOPO5 || ETOPO2 || SRTM30
      integer*2, allocatable ::       bathy2(:,:)
#endif
#endif
      BIGREAL         lat(ILAT), lon(ILON)
      common /ll/     lat, lon
      BIGREAL         corner(2,2), dx, dy
      BIGREAL         minlat, minlon, maxlat, maxlon

      integer   NSUB, LONSUB, LATSUB
      parameter(NSUB = ILON/250+ILAT/250)
      parameter(LONSUB=(ILON-1)/NSUB, LATSUB=(ILAT-1)/NSUB )
      real            tmp(LONSUB,LATSUB), htmp(0:L,0:M)
      character*24    lbl
      BIGREAL         alatn, alats, alonw, alone
      BIGREAL         depthmxr, depthmnr,                               &
     &                htmp1, htmp2, x1, y1
      BIGREAL         vmin, vmax
#if IBCAO
      real            vsimin, vsimax
#else
      integer*2       vsimin, vsimax
#endif
      integer         i, j, i1, j1
      BIGREAL         zero, one
      parameter     ( zero=0.d0, one=1.d0  )
      integer         find_i, find_j
      logical         mask_flag
#if IMG_AVG
      integer         i2, j2, ii, jj, np
      BIGREAL         aa, bb, cc, dd
#endif /* IMG_AVG */

#include "gridid.h"

      allocate(bathy(ILON,ILAT))
#if ETOPO5 || ETOPO2 || SRTM30
      allocate(bathy2(ILON-1,ILAT))
#endif
!  depths less than depthmin will be set to depthmin
      print *, 'minimum depth (negative is land elevation)'
      read(5,*) depthmin
!  depths greater than depthmax will be set to depthmax
      print *, 'maximum depth'
      read(5,*) depthmax
      print *, 'apply mask? (T/F)'
      read(5,*) mask_flag

#if ETOPO5 || ETOPO2 || SRTM30
      call extract(lon,lat,bathy2,ILON-1,ILAT)

!  Longitudes between 0 and 360
      lon(ILON) = lon(1)+360.
# if ETOPO5
      do i=1,ILON
        if (lon(i) .lt. 0.) lon(i) = lon(i) + 360.
      enddo
# endif

!  Copy value at one end of the Earth to other end
      do i=1,ILON-1
        do j=1,ILAT
           bathy(i,j) = bathy2(i,j)
        enddo
      enddo
      do j=1,ILAT
        bathy(ILON,j) = bathy(1,j)
      enddo
#else
      call extract(lon,lat,bathy,ILON,ILAT)
      print *, 'ILON, ILAT', ILON, ILAT

!  Longitudes between 0 and 360
# if !ETOPO1 && !IBCAO
      do i=1,ILON
        if (lon(i) .lt. 0.) lon(i) = lon(i) + 360.
      enddo
# endif
#endif   /* ETOPO5 || ETOPO2 */

!  Our convention has been negative for land elevations
      do j=1,ILAT
        do i=1,ILON
          bathy(i,j) = -bathy(i,j)
        enddo
      enddo

      alatn = lat(ILAT)
      alats = lat(1)
      alonw = vmin(lon,ILON)
      alone = vmax(lon,ILON)

      depthmxr = vsimax(bathy,ILON*ILAT)
      print *,' Maximum depth read is ',depthmxr
      depthmnr = vsimin(bathy,ILON*ILAT)
      print *,' Minimum depth read is ',depthmnr

      call get_lat

      minlat = vmin(lat_rho,Lp*Mp)
      minlon = vmin(lon_rho,Lp*Mp)
      maxlat = vmax(lat_rho,Lp*Mp)
      maxlon = vmax(lon_rho,Lp*Mp)

#if !ETOPO1 && !IBCAO
#if ETOPO5
      if (minlon .lt. 0.) then
#endif /* ETOPO5 */
        do j=0,M
          do i=0,L
            if (lon_rho(i,j) .lt. 0.) lon_rho(i,j) = lon_rho(i,j) + 360.
          enddo
        enddo
        minlon = vmin(lon_rho,Lp*Mp)
        maxlon = vmax(lon_rho,Lp*Mp)
#if ETOPO5
      endif
#endif /* ETOPO5 */
#endif /* ETOPO1 */

      print *,'latitude range of grid ',minlat,maxlat
      print *,'latitude range of data ',alats,alatn
      print *,'longitute range of grid ',minlon,maxlon
      print *,'longitute range of data ',alonw,alone
#if !defined IBCAO && !defined ARDEM
      if (minlat .lt. alats .or. maxlat .gt. alatn) then
        print *,'Insufficient latitude range of data'
        print *,'try again ... '
#if NO_EXIT
        stop
#else
        call exit(1)
#endif  /* NO_EXIT */
      end if
      if (minlon .lt. alonw .or. maxlon .gt. alone) then
        print *,'Insufficient longitute range of data'
        print *,'try again ... '
#if NO_EXIT
        stop
#else
        call exit(1)
#endif  /* NO_EXIT */
      end if
#endif  /* !IBCAO */

!  h at lat_rho, lon_rho, before h is clipped

#if IMG_AVG
!  h at lat_rho, lon_rho

      do j=1,Mm
        do i=1,Lm
          aa = lat_v(i,j)
          bb = lat_v(i,j+1)
          cc = lon_u(i,j)
          dd = lon_u(i+1,j)

          call locate(lon,ILON,cc,i1)
          call locate(lon,ILON,dd,i2)
          call locate(lat,ILAT,aa,j1)
          call locate(lat,ILAT,bb,j2)

          np=0
          h(i,j) = 0.
          do jj=j1+1,j2
            do ii=i1+1,i2
              h(i,j)=h(i,j)+bathy(ii,jj)
              if(bathy(ii,jj).ne.0.) np=np+1
            enddo
          enddo
          if (np .ne. 0) h(i,j)=h(i,j)/float(np)
        enddo
      enddo

! Valeurs au bord
      do i=1,Lm
        h(i,0)=h(i,1)
        h(i,M)=h(i,Mm)
      enddo
      do j=0,M
        h(0,j)=h(1,j)
        h(L,j)=h(Lm,j)
      enddo
#else
      if (evenflag) then
        do i=0,L
          do j=0,M
            i1 = int((lon_rho(i,j) - alonw)/DELTA_X) + 1
            j1 = int((lat_rho(i,j) - alats)/DELTA_Y) + 1
#if defined IBCAO
            h(i,j) = 0.0
            if (j1 < 1) cycle
#elif defined ARDEM
            h(i,j) = 0.0
            if (j1 < 1 .or. i1 < 1) cycle
            if (j1 >= ILAT .or. i1 >= ILON) cycle
#endif
            x1 = lon_rho(i,j) - (alonw + (i1-1.)*DELTA_X)
            y1 = lat_rho(i,j) - (alats + (j1-1.)*DELTA_Y)
!           print *, 'Interpolating', i, j, i1, j1, x1, y1
            corner(1,1) = bathy(i1,j1)
            corner(1,2) = bathy(i1,j1+1)
            corner(2,1) = bathy(i1+1,j1)
            corner(2,2) = bathy(i1+1,j1+1)
            htmp1 = corner(1,1) + (corner(2,1)-corner(1,1))*x1/DELTA_X
            htmp2 = corner(1,2) + (corner(2,2)-corner(1,2))*x1/DELTA_X
            h(i,j) = htmp1 + (htmp2-htmp1)*y1/DELTA_Y
          enddo
        enddo
      else
        do i=0,L
          do j=0,M
            i1 = find_i(lon_rho(i,j),ILON,lon)
            j1 = find_j(lat_rho(i,j),ILAT,lat)
            x1 = lon_rho(i,j) - lon(i1)
            y1 = lat_rho(i,j) - lat(j1)
            corner(1,1) = bathy(i1,j1)
            corner(1,2) = bathy(i1,j1+1)
            corner(2,1) = bathy(i1+1,j1)
            corner(2,2) = bathy(i1+1,j1+1)
            dx = lon(i1+1) - lon(i1)
            dy = lat(j1+1) - lat(j1)
            htmp1 = corner(1,1) + (corner(2,1)-corner(1,1))*x1/dx
            htmp2 = corner(1,2) + (corner(2,2)-corner(1,2))*x1/dx
            h(i,j) = htmp1 + (htmp2-htmp1)*y1/dy
          enddo
        enddo
      endif
#endif  /* IMG_AVG */

!  calculate mask from bathymetry
      if (mask_flag) then
        do i=0,L
          do j=0,M
#ifdef PALAU_BATHY
            if (h(i,j) .le. 50.) then
#else
            if (h(i,j) .le. 0.) then
#endif
              mask_rho(i,j) = zero
            else
              mask_rho(i,j) = one
            endif
          enddo
        enddo
        call uv_mask

!  Write mask
        call wrt_mask
      endif
!  Write unclipped bathymetry
      call wrt_h

!  Now clip bathy values before interpolating h again
      if (depthmnr .lt. depthmin) then
        do j=1,ILAT
          do i=1,ILON
            bathy(i,j) = max(depthmin,bathy(i,j))
          enddo
        enddo
        print *,' There were depths less than ',depthmin
        print *,' These were set to ',depthmin
      endif
      if (depthmxr .gt. depthmax) then
        do j=1,ILAT
          do i=1,ILON
            bathy(i,j) = min(depthmax,bathy(i,j))
          enddo
        enddo
        print *,' There were depths greater than ',depthmax
        print *,' These were set to ',depthmax
      endif

!  h at lat_rho, lon_rho, take 2
#if IMG_AVG
!  h at lat_rho, lon_rho
      do j=1,Mm
        do i=1,Lm
          aa = lat_v(i,j)
          bb = lat_v(i,j+1)
          cc = lon_u(i,j)
          dd = lon_u(i+1,j)

          call locate(lon,ILON,cc,i1)
          call locate(lon,ILON,dd,i2)
          call locate(lat,ILAT,aa,j1)
          call locate(lat,ILAT,bb,j2)

          np=0
          h(i,j) = 0.
          do jj=j1+1,j2
            do ii=i1+1,i2
              h(i,j)=h(i,j)+bathy(ii,jj)
              if(bathy(ii,jj).ne.0.) np=np+1
            enddo
          enddo
          if (np .ne. 0) h(i,j)=h(i,j)/float(np)
        enddo
      enddo

! Valeurs au bord
      do i=1,Lm
        h(i,0)=h(i,1)
        h(i,M)=h(i,Mm)
      enddo
      do j=0,M
        h(0,j)=h(1,j)
        h(L,j)=h(Lm,j)
      enddo
#else
      if (evenflag) then
        do i=0,L
          do j=0,M
            i1 = int((lon_rho(i,j) - alonw)/DELTA_X) + 1
            j1 = int((lat_rho(i,j) - alats)/DELTA_Y) + 1
#if defined IBCAO
            h(i,j) = 0.0
            if (j1 < 1) cycle
#elif defined ARDEM
            h(i,j) = 0.0
            if (j1 < 1 .or. i1 < 1) cycle
            if (j1 >= ILAT .or. i1 >= ILON) cycle
#endif
            x1 = lon_rho(i,j) - (alonw + (i1-1.)*DELTA_X)
            y1 = lat_rho(i,j) - (alats + (j1-1.)*DELTA_Y)
            corner(1,1) = bathy(i1,j1)
            corner(1,2) = bathy(i1,j1+1)
            corner(2,1) = bathy(i1+1,j1)
            corner(2,2) = bathy(i1+1,j1+1)
            htmp1 = corner(1,1) + (corner(2,1)-corner(1,1))*x1/DELTA_X
            htmp2 = corner(1,2) + (corner(2,2)-corner(1,2))*x1/DELTA_X
            h(i,j) = htmp1 + (htmp2-htmp1)*y1/DELTA_Y
          enddo
        enddo
      else
        do i=0,L
          do j=0,M
            i1 = find_i(lon_rho(i,j),ILON,lon)
            j1 = find_j(lat_rho(i,j),ILAT,lat)
            x1 = lon_rho(i,j) - lon(i1)
            y1 = lat_rho(i,j) - lat(j1)
            corner(1,1) = bathy(i1,j1)
            corner(1,2) = bathy(i1,j1+1)
            corner(2,1) = bathy(i1+1,j1)
            corner(2,2) = bathy(i1+1,j1+1)
            dx = lon(i1+1) - lon(i1)
            dy = lat(j1+1) - lat(j1)
            htmp1 = corner(1,1) + (corner(2,1)-corner(1,1))*x1/dx
            htmp2 = corner(1,2) + (corner(2,2)-corner(1,2))*x1/dx
            h(i,j) = htmp1 + (htmp2-htmp1)*y1/dy
          enddo
        enddo
      endif
#endif  /* IMG_AVG */

!  write clipped bathymetry
      call wrt_h

#if PLOTS
!  Plot the original topography
      call start_plot

!  set foreground to black, background to white
      call gscr(1,0,1.,1.,1.)
      call gscr(1,1,0.,0.,0.)

      call pcseti('QUALITY',1)
      write(lbl,200)
  200 format('Bottom Topography')
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,lbl,.012,0.,0.)
      do j=1,LATSUB
        do i=1,LONSUB
          tmp(i,j)=bathy(NSUB*i,NSUB*j)
        enddo
      enddo
      call set(.05,.95,.25,.95,1.,float(LONSUB),1.,float(LATSUB),1)
      call perim(1,ILON/12,1,ILAT/12)
      call cpseti('SET',0)
      call cpshift(tmp,LONSUB,LONSUB,LATSUB,0.,.false.,1.)
      call frame
!  set foreground to black, background to white
      call gscr(1,0,1.,1.,1.)
      call gscr(1,1,0.,0.,0.)
      call cpsfill(tmp,LONSUB,LONSUB,LATSUB,5,2,.false.,.false.)
      call frame

!  Plot the topography on the grid

      do j=0,M
        do i=0,L
          htmp(i,j)=h(i,j)
        enddo
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call plchhq(0.5,0.98,lbl,.012,0.,0.)
      call set(.05,.95,.05,.95,1.,float(Lp),1.,float(Mp),1)
      call perim(1,L,1,M)
      call cpcnrc(htmp,Lp,Lp,Mp,0.,0.,0.,1,-1,-682)
      call frame
      call end_plot
#endif  /* PLOTS */

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
      end

! ******************************************************************

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

! ***********************************************************

#if IBCAO
      real function vsimax(vect,N)
      integer   N
      real vect(N)
      real tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = max(tmp,vect(i))
      enddo
      vsimax = tmp
      return
      end

      real function vsimin(vect,N)
      integer   N
      real vect(N)
      real tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = min(tmp,vect(i))
      enddo
      vsimin = tmp
      return
      end
#else
      integer*2 function vsimax(vect,N)
      integer   N
      integer*2 vect(N)
      integer*2 tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = max(tmp,vect(i))
      enddo
      vsimax = tmp
      return
      end

      integer*2 function vsimin(vect,N)
      integer   N
      integer*2 vect(N)
      integer*2 tmp
      integer   i

      tmp = vect(1)
      do i=2,N
        tmp = min(tmp,vect(i))
      enddo
      vsimin = tmp
      return
      end
#endif

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine locate(xx,n,x,j)

      integer j,n
      BIGREAL xx(n)
      BIGREAL x
      integer jl,jm,ju

      jl=0
      ju=n+1
 10   if (ju-jl .gt. 1) then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
        goto 10
      endif
      j=jl
      return
      end

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      function find_i(lon_rho,ILON,lon)
      integer         find_i
      BIGREAL         lon_rho
      integer         ILON
      BIGREAL         lon(ILON)
      integer         mid, high, low

      low = 1
      high = ILON
      mid = (high+low)/2

      do while (.true.)
        if (lon_rho .gt. lon(mid)) then
          low = mid
          mid = (high+low)/2
          if (high - low .eq. 1) then
            find_i = low
            return
          endif
        else
          high = mid
          mid = (high+low)/2
          if (high - low .eq. 1) then
            find_i = low
            return
          endif
        endif
      enddo

      end

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      function find_j(lat_rho,ILAT,lat)
      integer         find_j
      BIGREAL         lat_rho
      integer         ILAT
      BIGREAL         lat(ILAT)
      integer         mid, high, low

      low = 1
      high = ILAT
      mid = (high+low)/2

      do while (.true.)
        if (lat_rho .gt. lat(mid)) then
          low = mid
          mid = (high+low)/2
          if (high - low .eq. 1) then
            find_j = low
            return
          endif
        else
          high = mid
          mid = (high+low)/2
          if (high - low .eq. 1) then
            find_j = low
            return
          endif
        endif
      enddo

      end
