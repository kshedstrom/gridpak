#include "griddefs.h"

      subroutine extract(lon,lat,topo,im,jm)
!
!=======================================================================
!  Define global data.
!=======================================================================
!
      use netcdf
!
!=======================================================================
!  Define local data.
!=======================================================================
!
      integer      ncinpid, rcode, i, j
      character*60 fname
      integer      im, jm, stdout, varid
      integer      count(2), start(2)
      integer*2    topo(im,jm)
      real         c0, c5, c60
      BIGREAL      lat(jm), lon(im)
      parameter (  c0=0.0,c5=5.0,c60=60.0,stdout=6)
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
#if ETOPO5
      call getenv('ETOPO5',fname)
#elif ETOPO2
      call getenv('ETOPO2',fname)
#elif ETOPO1
      call getenv('ETOPO1',fname)
#elif GEBCO
      call getenv('GEBCO',fname)
#elif IBCAO
      call getenv('IBCAO',fname)
#else
      call getenv('BATHY_FILE',fname)
#endif
      rcode = nf90_open(fname,nf90_nowrite,ncinpid)
      if (rcode.ne.0) then
        write(stdout,904) fname
        goto 10
      endif
!
!-----------------------------------------------------------------------
!  Extract longitudes of the requested area.
!-----------------------------------------------------------------------
!
      start(1)=1
      count(1)=im
      rcode = nf90_inq_varid(ncinpid,'lon',varid)
      if (rcode.eq.0) then
        rcode = nf90_get_var(ncinpid,varid,lon)
        if (rcode.ne.0) then
          write(stdout,900) 'lon'
          goto 10
        endif
      else
        write(stdout,901) 'lon'
      endif
!
!-----------------------------------------------------------------------
!  Extract latitudes of the requested area.
!-----------------------------------------------------------------------
!
      rcode = nf90_inq_varid(ncinpid,'lat',varid)
      if (rcode.eq.0) then
        rcode = nf90_get_var(ncinpid,varid,lat)
        if (rcode.ne.0) then
          write(stdout,900) 'lat'
          goto 10
        endif
      else
        write(stdout,901) 'lat'
      endif
      print *, 'im ', im, jm, lon(1)
!
!-----------------------------------------------------------------------
!  Extract bathmetry
!-----------------------------------------------------------------------
!
      rcode=nf90_inq_varid(ncinpid,'z',varid)
      if (rcode.eq.0) then
          rcode = nf90_get_var(ncinpid,varid,topo)
          if(rcode.ne.0) then
            write(stdout,900) 'z'
            goto 10
          endif
      else
        write(stdout,901) 'z'
      endif
!
      goto 20
  10  write(stdout,903)
      print *, trim(nf90_strerror(rcode))
      stop
  20  continue
 900  format(/' EXTRACT - error while reading variable: ',a)
 901  format(/' EXTRACT - cannot find variable: ',a)
 902  format(/' EXTRACT - error while writing variable: ',a)
 903  format(/' EXTRACT - terminated abnormally.')
 904  format(/' EXTRACT - error while opening file: ',a)
      return
      end
