#include "griddefs.h"
      program bloc_data
c       
#include "bathy.h"
#include "ncgrid.h"
      integer   incri, incrj, ibloc, jbloc
      parameter(incri=22,incrj=33)
      parameter(ibloc=4,jbloc=2)
      BIGREAL   tmp(0:L,0:M)
      integer   jmin, jmax, imin, imax, nblocj
      integer   i, j, nbloc

#include "gridid.h"

c Lecture de la bathymetrie            
      call get_h
      
c Division par 10 de la bathymetrie
      do i=0,L
         do j=0,M
            tmp(i,j) = h(i,j)/10.
         enddo
      enddo 
           
c Ecriture des blocs
      open(15,file='Bloc_data.h100m.br',form='formatted')
cm      open(15,file='Bloc_mask.etopo5',form='formatted')
      jmin=0
      jmax=incrj
      nblocj=1
      do while(jmin .lt. jbloc*incrj)
         nbloc=nblocj
         imin=0
         imax=incri  
         do while (imin .lt. ibloc*incri)
           write(15,99)nbloc 
           write(15,100)imin,imax
           write(15,101)jmin,jmax
           write(15,102)(i,i=imin,imax)
           write(15,*)'  '
           do j = jmax, jmin,-1
             write(15,990) j,(tmp(i,j),i=imin,imax)
           enddo
           write(15,*)'  '
           imin=imax
           imax=imax+incri
           nbloc=nbloc+1
           if (imax.gt.L) imax=L
         enddo
         jmin=jmax
         jmax=jmax+incrj
         nblocj=nbloc
         if (jmax.gt.M) jmax=M 
      enddo
      close(15)
99    format('Donnees du bloc No',i3)       
100   format('imin=',i3,x,'imax=',i3) 
101   format('jmin=',i3,x,'jmax=',i3) 
102   format(2x,31(x,i4))                  
990   format(x,i3,31(x,f4.0)) 
cm ---------------------------                
c      
      call crash(' ', 0)
      end

c ********************************************************************

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
