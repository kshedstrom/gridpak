#include "griddefs.h"
#include "bathy.h"
      integer         ITMAX, IBIG
      parameter (     ITMAX=8, IBIG=400    )
!  ITMAX is the number of iterations to perform
!  IBIG is the largest number of points to be read in for one
!  boundary.
!
!  original distribution of x,y points is preserved on boundary kb1
!  and kb2:
      integer         kb1, kb2
      parameter (     kb1 = 3, kb2 = 2   )

      integer         L2, M2, L2big, M2big, nwrk
      integer         N, N1, N2, N3, N4
      parameter (     L2=2*(L-1), M2=2*(M-1)   )
      parameter (     L2big=2*Lm, M2big=2*Mm   )
      parameter (     N1=M2, N2=M2+L2, N3=M2+L2+M2, &
     &                N4=M2+L2+M2+L2,  N=N4    )
      integer         KK
      parameter (     KK = 9   )
      parameter (     nwrk = 2*(KK-2)*(2**(KK+1)) + KK + 10*M2big + &
     &                       12*L2big + 27    )
      BIGREAL         sxi(0:L2big), seta(0:M2big)
      common / xiej / sxi, seta
      BIGREAL         x1spl(IBIG),x2spl(IBIG),x3spl(IBIG),x4spl(IBIG), &
     &                y1spl(IBIG),y2spl(IBIG),y3spl(IBIG),y4spl(IBIG), &
     &                s1spl(IBIG),s2spl(IBIG),s3spl(IBIG),s4spl(IBIG), &
     &                b1spl(IBIG),b2spl(IBIG),b3spl(IBIG),b4spl(IBIG), &
     &                c1spl(IBIG),c2spl(IBIG),c3spl(IBIG),c4spl(IBIG)
      integer         nb1pts,nb2pts,nb3pts,nb4pts
      common / bdata/ x1spl, x2spl, x3spl, x4spl, &
     &                y1spl, y2spl, y3spl, y4spl, &
     &                s1spl, s2spl, s3spl, s4spl, &
     &                b1spl, b2spl, b3spl, b4spl, &
     &                c1spl, c2spl, c3spl, c4spl, &
     &                nb1pts, nb2pts, nb3pts, nb4pts
!  The boundary values are read from stdin for edges which have
!  rbnd true.  For boundaries which are read in, the grid spacing
!  is proportional to distance along the boundary if even? is true.
!  Otherwise, it is proportional to the spacing of the supplied
!  boundary points.
      logical         rbnd1, rbnd2, rbnd3, rbnd4, &
     &                even1, even2, even3, even4
      parameter   (  rbnd1=.true., rbnd2=.true., &
     &               rbnd3=.true., rbnd4=.true., &
     &               even1=.true., even2=.true., &
     &               even3=.true., even4=.true.  )

!  The following are used when you need to fit a boundary with
!  bumps on opposite sides and need to make intermediate partial
!  grids.  Set pleft1 or pbot1 to true to print out the boundaries
!  of a partial left or bottom grid.  Set pleft2 or pbot2 to true
!  to print out the new left or bottom boundary.  Lmiddle or Mmiddle
!  gives the position of the interior boundary for the intermediate
!  grid.  The boundaries are written out to iout1 or iout2.
!
!  Don't forget to adjust the evenx flags, kb1 and kb2 accordingly.
      logical         pleft1, pleft2, pbot1, pbot2
      integer         Lmiddle, Mmiddle, iout1, iout2
      parameter   (   pleft1=.false., pleft2=.false., &
     &                pbot1=.false., pbot2=.false., &
     &                Lmiddle=11, Mmiddle=15, &
     &                iout1=13, iout2=14             )

!  These variables are used for writing out a subset of the psi points
!  to be used in generating a nested domain.
      logical        subset
      integer        Lwest, Least, Msouth, Mnorth, iout3
      parameter    ( subset = .false., Lwest = 14, Least = 17, &
     &               Msouth = 6, Mnorth = 14, iout3 = 15    )

!  xpots unit numbers
      integer         ipot1, ipot2, ipot3, ipot4
      parameter   (   ipot1=41, ipot2=42, ipot3=43, ipot4=44   )
