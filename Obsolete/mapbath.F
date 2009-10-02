      program mapbath
! *******************************************************************
!  From Ichiro Fukumori, edited by Kate Hedstrom to change
!  input file
!
! Map out bathymetry using NCAR contour routine STLINE-DRLINE
!     This is accomplished by replacing NCAR's FRSTD, VECTD, LASTD
!     to user supplied routines that just print out latitude,
!     longitude data to file.  Note also that user supplied FX, FY
!     functions that map data location to coordinate is needed.
!     Output is written to file 59. (58 is scratch file)
!     Format is number of points followed by that number of lat,
!     long pairs.
!
! Bathymetry data is subset of etopo5 as done by bathers
!
! *******************************************************************

#include "griddefs.h"
#include "bath2.h"

! Bathymetry data
      integer*2       bathy(ILON,ILAT)
      real            bot(ILON,ILAT)

! Common for functions fx, fy
      integer         i, j
      real            depth

!---------------------------------------------------------------------

      call extract(lon,lat,bathy,ILON,ILAT)

! parameters describing grid of bot
      xlonmi = lon(1)
      xlonma = lon(ILON)
      ylatmi = lat(1)
      ylatma = lat(ILAT)

      do j=1,ILAT
        do i=1,ILON
          bot(i,j) = bathy(i,j)
        enddo
      enddo

      print *,'Print depth of contour (positive is above sea level): '
      read *,depth
!     depth = -depth

! Start NCAR graphics
      call opngks

! Contour
      call stline (bot,ILON,ILON,ILAT,depth)

! End NCAR
      call clsgks

#if NO_EXIT
      stop
#else
      call exit(0)
#endif  /* NO_EXIT */
 210  call crash('read or write error while outputing solution', 0)
 9999 call crash(' Error reading in topography file ', 1)
      end

!----------------------------------------------------------------------

      subroutine crash(icrash,ierr)
      character*(*) icrash
      integer       ierr

      print *,icrash
      if (ierr .gt. 1) print *,'ierr = ',ierr
#if NO_EXIT
      stop
#else
      call exit(ierr)
#endif  /* NO_EXIT */
      return
      end

!----------------------------------------------------------------------

      real function fx(x,y)
      real          x, y
!
#include "bath2.h"
      real          xlon, ylat
      integer       i, j
!
      if (evenflag) then
        xlon=xlonmi+(xlonma-xlonmi)*(x-1.)/(float(ILON)-1.)
        ylat=ylatmi+(ylatma-ylatmi)*(y-1.)/(float(ILAT)-1.)
      else
        i = x
        j = y
        if (i .eq. ILON) then
          xlon = lon(ILON)
        else
          xlon = lon(i) + (lon(i+1)-lon(i))*(x-i)
          ylat = lat(j) + (lat(j+1)-lat(j))*(y-j)
        endif
      endif
      fx=xlon
      return
      end
!----------------------------------------------------------------------
      real function fy(x,y)
      real          x, y
!
#include "bath2.h"
      real          xlon, ylat
      integer       i, j
!
      if (evenflag) then
        xlon=xlonmi+(xlonma-xlonmi)*(x-1.)/(float(ILON)-1.)
        ylat=ylatmi+(ylatma-ylatmi)*(y-1.)/(float(ILAT)-1.)
      else
        i = x
        j = y
        if (j .eq. ILAT) then
          ylat = lat(ILAT)
        else
          xlon = lon(i) + (lon(i+1)-lon(i))*(x-i)
          ylat = lat(j) + (lat(j+1)-lat(j))*(y-j)
        endif
      endif
      fy=ylat
      return
      end
!----------------------------------------------------------------------
      subroutine frstd (x,y)
      real       x, y
      write (58,5901) y,x
 5901 format (2f9.3)
      return
      end
      subroutine vectd (x,y)
      real       x, y
      write (58,5901) y,x
 5901 format (2f9.3)
      return
      end
      subroutine lastd
      write (58,5903)
 5903 format (' 999.99  999.99')
      return
      end
!----------------------------------------------------------------------
      subroutine drline (Z,L,MM,NN)
      save
      integer    L, MM, NN
      real       Z(L,NN)
!
! THIS ROUTINE TRACES A CONTOUR LINE WHEN GIVEN THE BEGINNING BY STLINE.
! TRANSFORMATIONS CAN BE ADDED BY DELETING THE STATEMENT FUNCTIONS FOR
! FX AND FY IN DRLINE AND MINMAX AND ADDING EXTERNAL FUNCTIONS.
! X=1. AT Z(1,J), X=FLOAT(M) AT Z(M,J). X TAKES ON NON-INTEGER VALUES.
! Y=1. AT Z(I,1), Y=FLOAT(N) AT Z(I,N). Y TAKES ON NON-INTEGER VALUES.
!
      integer         IX, IY, IDX, IDY, IS, ISS, NP, INX(8), INY(8), &
     &                IR(50000), NR
      real            CV
      common /conre2/ IX         ,IY         ,IDX        ,IDY        , &
     &                IS         ,ISS        ,NP         ,CV         , &
     &                INX        ,INY        ,IR          ,NR
      integer         IOFFP
      real            SPVAL
      common /conre1/ IOFFP      ,SPVAL
      integer         IXBITS     ,IYBITS
      common /conre3/ IXBITS     ,IYBITS
      logical         IPEN       ,IPENO
      integer         IXX, IYY, ishift, M, N, JUMP1, JUMP2, IX0, IY0, &
     &                IS0, IX2, IY2, ISUB, ISBIG, IX3, IY3, IX4, IY4, &
     &                IXYPAK
      real            P1, P2, X, Y, C, fx, fy, XOLD, YOLD
!
!
!     STATEMENT FUNCTIONS TO BE REPLACED IF NON-IDENTITY
!     TRANSFORMATIONS ARE DESIRED.
!
!      fX(X,Y) = X
!      fy(X,Y) = Y
!
      IXYPAK(IXX,IYY) = ishift(IXX,IYBITS)+IYY
      C(P1,P2) = (P1-CV)/(P1-P2)
!
      DATA IPEN,IPENO/.TRUE.,.TRUE./
!
      M = MM
      N = NN
      if (IOFFP .EQ. 0) GO TO 101
      assign 110 to JUMP1
      assign 115 to JUMP2
      go to 102
  101 assign 112 to JUMP1
      assign 117 to JUMP2
  102 IX0 = IX
      IY0 = IY
      IS0 = IS
      if (IOFFP .eq. 0) go to 103
      IX2 = IX+INX(IS)
      IY2 = IY+INY(IS)
      IPEN = Z(IX,IY).ne.SPVAL .and. Z(IX2,IY2).ne.SPVAL
      IPENO = IPEN
  103 if (IDX .eq. 0) go to 104
      Y = IY
      ISUB = IX+IDX
      X = C(Z(IX,IY),Z(ISUB,IY))*float(IDX)+float(IX)
      GO TO 105
  104 X = IX
      ISUB = IY+IDY
      Y = C(Z(IX,IY),Z(IX,ISUB))*float(IDY)+float(IY)
  105 call frstd (fx(X,Y),fy(X,Y))
  106 IS = IS+1
      if (IS .gt. 8) IS = IS-8
      IDX = INX(IS)
      IDY = INY(IS)
      IX2 = IX+IDX
      IY2 = IY+IDY
      if (ISS .ne. 0) go to 107
      if (IX2.gt.M .or. IY2.gt.N .or. IX2.lt.1 .or. IY2.lt.1) go to 120
  107 if (CV-Z(IX2,IY2)) 108,108,109
  108 IS = IS+4
      IX = IX2
      IY = IY2
      go to 106
  109 if (IS/2*2 .eq. IS) go to 106
      go to JUMP1,(110,112)
  110 ISBIG = IS+(8-IS)/6*8
      IX3 = IX+INX(ISBIG-1)
      IY3 = IY+INY(ISBIG-1)
      IX4 = IX+INX(ISBIG-2)
      IY4 = IY+INY(ISBIG-2)
      IPENO = IPEN
      if (ISS .ne. 0) go to 111
      if (IX3.gt.M .or. IY3.gt.N .or. IX3.LT.1 .or. IY3.LT.1) go to 120
      if (IX4.gt.M .or. IY4.gt.N .or. IX4.LT.1 .or. IY4.LT.1) go to 120
  111 IPEN = Z(IX,IY).ne.SPVAL .and. Z(IX2,IY2).ne.SPVAL .and. &
     &       Z(IX3,IY3).ne.SPVAL .and. Z(IX4,IY4).ne.SPVAL
  112 if (IDX .eq. 0) go to 113
      Y = IY
      ISUB = IX+IDX
      X = C(Z(IX,IY),Z(ISUB,IY))*float(IDX)+float(IX)
      go to 114
  113 X = IX
      ISUB = IY+IDY
      Y = C(Z(IX,IY),Z(IX,ISUB))*float(IDY)+float(IY)
  114 go to JUMP2,(115,117)
  115 if (.not.IPEN) go to 118
      if (IPENO) go to 116
!
! END OF LINE SEGMENT
!
      call lastd
      call frstd (fx(XOLD,YOLD),fy(XOLD,YOLD))
!
! CONTINUE LINE SEGMENT
!
  116 continue
  117 call vectd (fx(X,Y),fy(X,Y))
  118 XOLD = X
      YOLD = Y
      if (IS .ne. 1) go to 119
      NP = NP+1
      if (NP .gt. NR) go to 120
      IR(NP) = IXYPAK(IX,IY)
  119 if (ISS .eq. 0) go to 106
      if (IX.ne.IX0 .or. IY.ne.IY0 .or. IS.ne.IS0) go to 106
!
! END OF LINE
!
  120 call lastd
      return
      end
!----------------------------------------------------------------------
      subroutine stline (Z,LL,MM,NN,CONV)
      save
      integer    LL, MM, NN
      real       Z(LL,NN), CONV
!
! THIS ROUTINE FINDS THE BEGINNINGS OF ALL CONTOUR LINES AT LEVEL CONV.
! FIRST THE EDGES ARE SEARCHED FOR LINES INTERSECTING THE EDGE (OPEN
! LINES) THEN THE INTERIOR IS SEARCHED FOR LINES WHICH DO NOT INTERSECT
! THE EDGE (CLOSED LINES).  BEGINNINGS ARE STORED IN IR TO PREVENT RE-
! TRACING OF LINES.  IF IR IS FILLED, THE SEARCH IS STOPPED FOR THIS
! CONV.
!
      integer         IX, IY, IDX, IDY, IS, ISS, NP, INX(8), INY(8), &
     &                IR(50000), NR
      real            CV
      common /conre2/ IX         ,IY         ,IDX        ,IDY        , &
     &                IS         ,ISS        ,NP         ,CV         , &
     &                INX        ,INY        ,IR          ,NR
      integer         IXBITS     ,IYBITS
      common /conre3/ IXBITS     ,IYBITS
      integer         IXX, IYY, I, J, K, L, M, N, IP1, JP1, IXY, &
     &                IXYPAK, LTYPE, ishift, IUNIT, i1mach
      real            VXA, VXB, VYA, VYB, XA, XB, YA, YB, X, Y
!
!
!
!
!
!
      IXYPAK(IXX,IYY) = ishift(IXX,IYBITS)+IYY
!
      L = LL
      M = MM
      N = NN
      CV = CONV
      NP = 0
      ISS = 0
      do 102 IP1=2,M
         I = IP1-1
         if (Z(I,1).ge.CV .or. Z(IP1,1).lt.CV) go to 101
         IX = IP1
         IY = 1
         IDX = -1
         IDY = 0
         IS = 1
         call drline (Z,L,M,N)
  101    if (Z(IP1,N).ge.CV .or. Z(I,N).lt.CV) go to 102
         IX = I
         IY = N
         IDX = 1
         IDY = 0
         IS = 5
         call drline (Z,L,M,N)
  102 continue
      do 104 JP1=2,N
         J = JP1-1
         if (Z(M,J).ge.CV .or. Z(M,JP1).lt.CV) go to 103
         IX = M
         IY = JP1
         IDX = 0
         IDY = -1
         IS = 7
         call drline (Z,L,M,N)
  103    if (Z(1,JP1).ge.CV .or. Z(1,J).lt.CV) go to 104
         IX = 1
         IY = J
         IDX = 0
         IDY = 1
         IS = 3
         call drline (Z,L,M,N)
  104 continue
      ISS = 1
      do 108 JP1=3,N
         J = JP1-1
         do 107 IP1=2,M
            I = IP1-1
            if (Z(I,J).ge.CV .or. Z(IP1,J).lt.CV) go to 107
            IXY = IXYPAK(IP1,J)
            if (NP .eq. 0) go to 106
            do 105 K=1,NP
               if (IR(K) .eq. IXY) go to 107
  105       continue
  106       NP = NP+1
            if (NP .gt. NR) then
!
! THIS PRINTS AN ERROR MESSAGE IF THE LOCAL ARRAY IR IN SUBROUTINE
! STLINE HAS AN OVERFLOW
! THIS MESSAGE IS WRITTEN BOTH ON THE FRAME AND ON THE STANDARD ERROR
! UNIT
!
      IUNIT = i1mach(4)
      write(IUNIT,2000)
 2000 format( &
     &' WARNING FROM ROUTINE STLINE IN CONREC--WORK ARRAY OVERFLOW')
      call getset(VXA,VXB,VYA,VYB,XA,XB,YA,YB,LTYPE)
      Y = (YB - YA) / 2.
      X = (XB - XA) / 2.
      call pwrit(X,Y, &
     &'**WARNING--PICTURE INCOMPLETE**', &
     & 31,3,0,0)
      Y = Y * .7
      call pwrit(X,Y, &
     &'WORK ARRAY OVERFLOW IN STLINE', &
     & 29,3,0,0)
        return
        endif
            IR(NP) = IXY
            IX = IP1
            IY = J
            IDX = -1
            IDY = 0
            IS = 1
            call drline (Z,L,M,N)
  107    continue
  108 continue
      return
      end
      blockdata conbd
      integer         IOFFP
      real            SPVAL
      common /conre1/ IOFFP      ,SPVAL
      integer         IX, IY, IDX, IDY, IS, ISS, NP, INX(8), INY(8), &
     &                IR(50000), NR
      real            CV
      common /conre2/ IX         ,IY         ,IDX        ,IDY        , &
     &                IS         ,ISS        ,NP         ,CV         , &
     &                INX        ,INY        ,IR          ,NR
      integer         ISIZEL     ,ISIZEM     ,ISIZEP     ,NREP       , &
     &                NCRT       ,ILAB       ,NULBLL     ,IOFFD      , &
     &                IOFFM      ,ISOLID     ,NLA        ,NLM
      real            EXT, XLT, YBT, SIDE
      common /conre4/ ISIZEL     ,ISIZEM     ,ISIZEP     ,NREP       , &
     &                NCRT       ,ILAB       ,NULBLL     ,IOFFD      , &
     &                EXT        ,IOFFM      ,ISOLID     ,NLA        , &
     &                NLM        ,XLT        ,YBT        ,SIDE
      integer         IRECMJ     ,IRECMN     ,IRECTX
      common /RECINT/ IRECMJ     ,IRECMN     ,IRECTX
      data IOFFP,SPVAL/1,99999./
      data ISIZEL,ISIZEM,ISIZEP,NLA,NLM,XLT,YBT,SIDE,ISOLID,NREP,NCRT/ &
     &     0       ,  0,   0, 16, 40,.05,.05,  .9,  1023,   6,   4 /
      data EXT,IOFFD,NULBLL,IOFFM,ILAB/.25,0,2,1,0/
      data INX(1),INX(2),INX(3),INX(4),INX(5),INX(6),INX(7),INX(8)/ &
     &        -1 ,   -1 ,    0 ,    1 ,    1 ,    1 ,    0 ,   -1 /
      data INY(1),INY(2),INY(3),INY(4),INY(5),INY(6),INY(7),INY(8)/ &
     &         0 ,    1 ,    1 ,    1 ,    0 ,   -1 ,   -1 ,   -1 /
      data NR/50000/
      data IRECMJ,IRECMN,IRECTX/ 1 ,   1 ,   1/
!
! REVISION HISTORY---
!
! JANUARY 1980     ADDED REVISION HISTORY AND CHANGED LIBRARY NAME
!                  FROM CRAYLIB TO PORTLIB FOR MOVE TO PORTLIB
!
! MAY 1980         ARRAYS IWORK AND ENCSCR, PREVIOUSLY TOO SHORT FOR
!                  SHORT-WORD-LENGTH MACHINES, LENGTHENED.  SOME
!                  DOCUMENTATION CLARIFIED AND CORRECTED.
!
! JUNE 1984        CONVERTED TO FORTRAN 77 AND TO GKS
!
! JUNE 1985        ERROR HANDLING LINES ADDED; IF OVERFLOW HAPPENS TO
!                  WORK ARRAY IN STLINE, A WARNING MESSAGE IS WRITTEN
!                  BOTH ON PLOT FRAME AND ON STANDARD ERROR MESSAGE.
!-------------------------------------------------------------------
!
      end

