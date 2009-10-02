      program lectgrid
#include "griddefs.h"
#include "bathy.h"
#include "ncgrid.h"
      integer         i, j
      
#include "gridid.h"

      call get_lat
      
      print*,'rholat(',L,',',M,')=',lat_rho(L,M)
      print*,'rholon(0,0)=',lon_rho(0,0)
      print*,'psilat(',L,',',M,')=',lat_psi(L,M)
      print*,'psilon(1,1)=',lon_psi(1,1)
      
      open(1,file='rlon.ce',form='formatted')
      write(1,110)
      do i=0,L
         write(1,210) i,lon_rho(i,M)
      enddo
      do j=0,M
         write(1,211) j,lat_rho(L,j)
      enddo
      close(1)
      
100   format('Latitude (degre) des points T : " rholat(j) j=1,353"')
110   format('Longitude (degre) des points T : "rholon(i) i=0,353"')
200   format(7(x,f10.5))
210   format('i=',i3,2x,'rholon=',f11.4)
211   format('i=',i3,2x,'rholat=',f11.4)
220   continue
      
      call crash(' ', 0)
      end

c ******************************************************************

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

