! *********************************************************************
!  Copyright (c) 1991, 1993, 1995 Rutgers University
!
! *********************************************************************

      subroutine start_plot

#include "griddefs.h"
#if PLOTS

!  Open and activate the GKS "workstation"
      call gopks(6,0)
#if X_WIN
      call gopwk(1,2,8)
#else
      call gopwk(1,2,1)
#endif  /* X_WIN */
      call gacwk(1)

!  Set some plotting parameters
      call pcseti('QUALITY',1)      ! medium quality font
!     call pcseti('CD - complex/duplex',1)    ! complex or duplex font

      return
      end

      subroutine end_plot
      call gdawk(1)
      call gclwk(1)
      call gclks

#endif          /* PLOTS */
      return
      end
