#include "gridparam.h"
      integer         L, M, Lp, Mp, L2d
      parameter (     L=Lm+1, M=Mm+1, Lp=Lm+2,  Mp=Mm+2  )
      parameter (     L2d=Lp*Mp                               )
!  lcflag tells which color bar to use, 1 - 6 (so far).  If negative
!  then black outlines are drawn over the color contour regions.
!  5 is shades of grey, 2 is John's Gebco chart attempt.
      integer         lcflag
      parameter (     lcflag=-2                               )
      real            x_v(0:Lm+3,0:Mm+3), y_v(0:Lm+3,0:Mm+3)
      common /xxyys/  x_v, y_v
      BIGREAL         xp(L,M), yp(L,M), xr(0:L,0:M), yr(0:L,0:M), &
     &                xu(L,0:M), yu(L,0:M), xv(0:L,M), yv(0:L,M), &
     &                xl, el
      real            xmin, ymin, xmax, ymax
      common /grdpts/ xp, yp, xr, yr, xu, yu, xv, yv, xl, el, &
     &                xmin, ymin, xmax, ymax
      BIGREAL         f(0:L,0:M), h(0:L,0:M)
      common /parm/   f, h
      BIGREAL         pm(0:L,0:M), pn(0:L,0:M), &
     &                dndx(0:L,0:M), dmde(0:L,0:M)
      common /pmpn/   pm, pn, dndx, dmde
      BIGREAL         lat_psi(L,M), lon_psi(L,M), &
     &                lat_rho(0:L,0:M), lon_rho(0:L,0:M), &
     &                lat_u(L,0:M), lon_u(L,0:M), &
     &                lat_v(0:L,M), lon_v(0:L,M)
      common /latlon/ lat_psi, lon_psi, lat_rho, lon_rho, &
     &                lat_u, lon_u, lat_v, lon_v
      BIGREAL         mask_rho(0:L,0:M), mask_u(L,0:M), &
     &                mask_v(0:L,M), mask_psi(L,M)
      common /rmask/  mask_rho, mask_u, mask_v, mask_psi
      BIGREAL         angle(0:L,0:M)
      common /angles/ angle
      integer*2       depthmin, depthmax
      common /hmins/  depthmin, depthmax
      integer         spherical
      common /logic/  spherical
