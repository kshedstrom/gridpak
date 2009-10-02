
      subroutine checkdefs
#include "grid.h"
#include "ncgrid.h"

      integer     is

      CPPoptions=' '
      write(stdout,10)
  10  format(/,' Activated C-preprocessing Options:',/)
  20  format(2x,a,t20,a)

#if DCOMPLEX
      write(stdout,20) 'DCOMPLEX', &
     &                 'Double complex'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+10)=' DCOMPLEX,'
#endif /* DCOMPLEX */
#if DBLEPREC
      write(stdout,20) 'DBLEPREC', &
     &                 'Double precision'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+10)=' DBLEPREC,'
#endif /* DBLEPREC */
#if ELLIPSOID
      write(stdout,20) 'ELLIPSOID', &
     &                 'Distances computed with ellipsoidal terms'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+11)=' ELLIPSOID,'
#endif /* ELLIPSOID */
#if ETOPO5
      write(stdout,20) 'ETOPO5', &
     &                 'ETOPO5 bathymetry'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' ETOPO5,'
#endif /* ETOPO5 */
#if ETOPO2
      write(stdout,20) 'ETOPO2', &
     &                 'ETOPO2 bathymetry'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' ETOPO2,'
#endif /* ETOPO2 */
#if DRAW_COASTS
      write(stdout,20) 'DRAW_COASTS', &
     &                 'Draw the coastlines on some plots'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+13)=' DRAW_COASTS,'
#endif /* DRAW_COASTS */
#if IMG_AVG
      write(stdout,20) 'IMG_AVG', &
     &                 'Average bathymetry within a gridbox'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+9)=' IMG_AVG,'
#endif /* IMG_AVG */
#if KEEP_SHALLOW
      write(stdout,20) 'KEEP_SHALLOW', &
     &                 'Keep shallow areas shallow'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+14)=' KEEP_SHALLOW,'
#endif /* KEEP_SHALLOW */
#if PLOTS
      write(stdout,20) 'PLOTS', &
     &                 'Write out NCAR graphics plots'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+7)=' PLOTS,'
#endif /* PLOTS */
#if SYS_POTS
      write(stdout,20) 'SYS_POTS', &
     &                 'System call to run xpots (not working)'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+10)=' SYS_POTS,'
#endif /* SYS_POTS */
#if X_WIN
      write(stdout,20) 'X_WIN', &
     &                 'X windows instead of a metafile'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+7)=' X_WIN,'
#endif /* X_WIN */
#if XPOTS1
      write(stdout,20) 'XPOTS1', &
     &                 'Read xpots file number 1'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' XPOTS1,'
#endif /* XPOTS1 */
#if XPOTS2
      write(stdout,20) 'XPOTS2', &
     &                 'Read xpots file number 2'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' XPOTS2,'
#endif /* XPOTS2 */
#if XPOTS3
      write(stdout,20) 'XPOTS3', &
     &                 'Read xpots file number 3'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' XPOTS3,'
#endif /* XPOTS3 */
#if XPOTS4
      write(stdout,20) 'XPOTS4', &
     &                 'Read xpots file number 4'
      is=LEN_TRIM(CPPoptions)+1
      CPPoptions(is:is+8)=' XPOTS4,'
#endif /* XPOTS4 */

      return
      end
