      subroutine uv_mask

!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
! *******************************************************************

#include "griddefs.h"
#include "bathy.h"
#include "ncgrid.h"
      integer         i, j
      BIGREAL         a1, a2, av2

      av2(a1,a2) = .5*(a1+a2)

#include "gridid.h"

!  Compute u, v masks.
      do i=1,L
        do j=0,M
          mask_u(i,j) = mask_rho(i,j) * mask_rho(i-1,j)
        enddo
      enddo
      do i=0,L
        do j=1,M
          mask_v(i,j) = mask_rho(i,j) * mask_rho(i,j-1)
        enddo
      enddo

!  Compute psi mask.
      do i=1,L
        do j=1,M
          mask_psi(i,j) = mask_rho(i,j) * mask_rho(i-1,j) &
     &                  * mask_rho(i,j-1) * mask_rho(i-1,j-1)
        enddo
      enddo

      return
      end
