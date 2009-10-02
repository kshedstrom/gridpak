      integer         L, M
      parameter (     L=72      , M=43      )
!
!  mud2 requires that these values satisfy
!
!     L = NXL*2**(NSTEP-1)+1
!     M = NYL*2**(NSTEP-1)+1
!
!  where NXL, NYL and NSTEP are integers.  Try to have NSTEP as large as
!  possible (see mud2 documentation).
!
!  subroutine 'factor' now takes care of passing NXL, NYL and NSTEP to mud2.
