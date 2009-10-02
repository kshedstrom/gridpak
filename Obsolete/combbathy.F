#include "griddefs.h"
      program combl_bathy

c Ce programme permet de combler :
c - le Gulf du mexique a depthgm
c - la mer des antilles a depthgm
c - le pacifique a depthmin
       
#include "bathy.h"
#include "ncgrid.h"
      real       depthmin, depthgm
      parameter (depthmin=200.,depthgm=3000.)
!* Mediterranee
      integer    imedmin1, imedmax1, jmedmin1, jmedmax1
      parameter (imedmin1=183,imedmax1=L,jmedmin1=48,jmedmax1=70)
!* Pacifique
      integer    ipmin, ipmax, jpmin, jpmax
      parameter (ipmin=0,ipmax=41,jpmin=0,jpmax=78)
!* Gulf du mexique
      integer    igmmin, igmmax, jgmmin, jgmmax
      parameter (igmmin=0,igmmax=35,jgmmin=28,jgmmax=47)
!* Mer des Antilles 1
      integer    imamin1, imamax1, jmamin1, jmamax1
      parameter (imamin1=16,imamax1=55,jmamin1=8,jmamax1=27)
!* Mer des Antilles 2
      integer    imamin2, imamax2, jmamin2, jmamax2
      parameter (imamin2=56,imamax2=69,jmamin2=8,jmamax2=23)
      integer    i, j, imax0, jmax0, imaxj
      logical    pacific_auto

#include "gridid.h"

      pacific_auto=.true.
      
c Lecture des donnees.            
      call get_h
      
c GULF DU MEXIQUE a mettre a 3000 m maxi
c lat : 20N  ---> 30N 
c lon : 99.8W ---> 80W 
        do i=igmmin,igmmax
          do j=jgmmin,jgmmax
            if (h(i,j).gt.depthgm) h(i,j)=depthgm
          enddo
        enddo
      
c MER DES ANTILLES a mettre a 3000 m maxi
c premier cadran
c lat : 10N ---> 20N 
c lon : 90W ---> 70W 
        do i=imamin1,imamax1
          do j=jmamin1,jmamax1
            if (h(i,j).gt.depthgm) h(i,j)=depthgm
          enddo
        enddo
c deuxieme cadran
c lat : 10N ---> 18N 
c lon : 70W ---> 63W 
        do i=imamin2,imamax2
          do j=jmamin2,jmamax2
            if (h(i,j).gt.depthgm) h(i,j)=depthgm
          enddo
        enddo 
        
c DOMAINE DE LA MEDITERRANNEE a mettre a depthmin
c premier cadran
c lat : 30N ---> 40N 
c lon : 6W ---> 0W 
        do i=imedmin1,imedmax1
          do j=jmedmin1,jmedmax1
            if (h(i,j).gt.depthmin) h(i,j)=depthmin
          enddo
        enddo 
c deuxieme cadran
c lat : 30N ---> 45N 
c lon : 0W ---> 17.8W 
!       do i=imedmin2,imedmax2
!         do j=jmedmin2,jmedmax2
!           if (h(i,j).gt.depthmin) h(i,j)=depthmin
!         enddo
!       enddo 
        
c DOMAINE DU PACIFIQUE a mettre a 200 m (ou 100m)
c lat : 17.7S ---> 8N 
c lon : 99.8W ---> 60W          
      if (pacific_auto.eq..false.) then 
           print*,'pacific_auto est faux '
        do i=ipmin,ipmax
          do j=jpmin,jpmax
            h(i,j)=depthmin
          enddo
        enddo

      else
      
	do j=0,M
	   if(h(0,j).le.depthmin) then
	     jmax0=j
 	     goto 11
 	   endif
 	enddo
  11	print *,'jmax0=',jmax0
	do i=0,L
	   if(h(i,0).le.depthmin) then
	     imax0=i
 	     goto 13
 	   endif
        enddo
  13	print *,'imax0=',imax0  	
	do j=0,jmax0
	   do i=0,imax0
	      if (h(i,j).le.depthmin) then
	         imaxj = i
c	         print*,'imaxj = ',imaxj
	         goto 22
	      endif
	   enddo
  22	   do i=0,imaxj
              h(i,j)=depthmin
 	   enddo
        enddo
        
      endif      ! pacific_auto
      
c Sortie de la bathymetrie comblee
      call wrt_h

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
