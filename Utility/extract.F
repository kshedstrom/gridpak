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
#elif GEBCO
      call getenv('GEBCO',fname)
#else
      call getenv('BATHY_FILE',fname)
#endif  /* ETOPO5 */
      rcode = nf90_open(fname,nf90_nowrite,ncinpid)
!
!-----------------------------------------------------------------------
!  Extract longitudes of the requested area.
!-----------------------------------------------------------------------
!
      start(1)=1
      count(1)=im
      rcode = nf90_inq_varid(ncinpid,'topo_lon',varid)
      if (rcode.eq.0) then
        rcode = nf90_get_var(ncinpid,varid,lon)
        if (rcode.ne.0) then
          write(stdout,900) 'topo_lon'
          goto 10
        endif
      else
        write(stdout,901) 'topo_lon'
      endif
!
!-----------------------------------------------------------------------
!  Extract latitudes of the requested area.
!-----------------------------------------------------------------------
!
      rcode = nf90_inq_varid(ncinpid,'topo_lat',varid)
      if (rcode.eq.0) then
        rcode = nf90_get_var(ncinpid,varid,lat)
        if (rcode.ne.0) then
          write(stdout,900) 'topo_lat'
          goto 10
        endif
      else
        write(stdout,901) 'topo_lat'
      endif
      print *, 'im ', im, jm, lon(1)
!
!-----------------------------------------------------------------------
!  Extract bathmetry
!-----------------------------------------------------------------------
!
      rcode=nf90_inq_varid(ncinpid,'topo',varid)
      if (rcode.eq.0) then
          rcode = nf90_get_var(ncinpid,varid,topo)
          if(rcode.ne.0) then
            write(stdout,900) 'topo'
            goto 10
          endif
      else
        write(stdout,901) 'topo'
      endif
!
      goto 20
  10  write(stdout,903)
      stop
  20  continue
 900  format(/' EXTRACT - error while reading variable: ',a)
 901  format(/' EXTRACT - cannot find variable: ',a)
 902  format(/' EXTRACT - error while writing variable: ',a)
 903  format(/' EXTRACT - terminated abnormally.')
      return
      end
