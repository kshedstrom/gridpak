#include "griddefs.h"
! ***********************************************************

      BIGREAL function vmax(vect,N)
      integer   N
      BIGREAL   vect(N)
      BIGREAL   tmp
      integer   i

      tmp = vect(1)
      do 100 i=2,N
        tmp = max(tmp,vect(i))
  100 continue
      vmax = tmp
      return
      end

      BIGREAL function vmin(vect,N)
      integer   N
      BIGREAL   vect(N)
      BIGREAL   tmp
      integer   i

      tmp = vect(1)
      do 100 i=2,N
        tmp = min(tmp,vect(i))
  100 continue
      vmin = tmp
      return
      end

