!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                        f i s h p a k                          *
!     *                                                               *
!     *                                                               *
!     *     a package of fortran subprograms for the solution of      *
!     *                                                               *
!     *      separable elliptic partial differential equations        *
!     *                                                               *
!     *                  (version 3.2 , november 1988)                *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        john adams, paul swarztrauber and roland sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the national center for atmospheric research          *
!     *                                                               *
!     *                boulder, colorado  (80307)  u.s.a.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the national science foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! package gnbnaux
!
! latest revision        november 1988
!
! purpose                to provide auxiliary routines for fishpak
!                        entries genbun and poistg.
!
! usage                  there are no user entries in this package.
!                        the routines in this package are not intended
!                        to be called by users, but rather by routines
!                        in packages genbun and poistg.
!
! special conditions     none
!
! i/o                    none
!
! precision              single
!
! required library       comf from fishpak
! files
!
! language               fortran
!
! history                written in 1979 by roland sweet of ncar's
!                        scientific computing division.  made available
!                        on ncar's public libraries in january, 1980.
!
! portability            fortran 77
! ********************************************************************
#include "griddefs.h"
      subroutine cosgen (n,ijump,fnum,fden,a)
      integer         n, ijump
      BIGREAL         fnum, fden, a(*)

!  local variables
      BIGREAL         pi, dum, pimach, pibyn, x, y
      integer         i, k, k1, k2, k3, k4, k5, np1
!
!
!     this subroutine computes required cosine values in ascending
!     order.  when ijump .gt. 1 the routine computes values
!
!        2*cos(j*pi/l) , j=1,2,...,l and j .ne. 0(mod n/ijump+1)
!
!     where l = ijump*(n/ijump+1).
!
!
!     when ijump = 1 it computes
!
!            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... ,n
!
!     where
!        fnum = 0.5, fden = 0.0,  for regular reduction values
!        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
!        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
!        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
!                                in poisn2 only.
!
!
      pi = pimach(dum)
      if (n .eq. 0) go to 105
      if (ijump .eq. 1) go to 103
      k3 = n/ijump+1
      k4 = k3-1
      pibyn = pi/float(n+ijump)
      do 102 k=1,ijump
         k1 = (k-1)*k3
         k5 = (k-1)*k4
         do 101 i=1,k4
            x = k1+i
            k2 = k5+i
            a(k2) = -2.*cos(x*pibyn)
  101    continue
  102 continue
      go to 105
  103 continue
      np1 = n+1
      y = pi/(float(n)+fden)
      do 104 i=1,n
         x = float(np1-i)-fnum
         a(i) = 2.*cos(x*y)
  104 continue
  105 continue
      return
      end

      subroutine merge (tcos,i1,m1,i2,m2,i3)
      integer         i1, m1, i2, m2, i3
      BIGREAL         tcos(*)
!
!     this subroutine merges two ascending strings of numbers in the
!     array tcos.  the first string is of length m1 and starts at
!     tcos(i1+1).  the second string is of length m2 and starts at
!     tcos(i2+1).  the merged string goes into tcos(i3+1).
!
!  local variables
      integer       j, k, l, m, j1, j2
      BIGREAL       x, y
!
      j1 = 1
      j2 = 1
      j = i3
      if (m1 .eq. 0) go to 107
      if (m2 .eq. 0) go to 104
  101 j = j+1
      l = j1+i1
      x = tcos(l)
      l = j2+i2
      y = tcos(l)
      if (x-y) 102,102,103
  102 tcos(j) = x
      j1 = j1+1
      if (j1 .gt. m1) go to 106
      go to 101
  103 tcos(j) = y
      j2 = j2+1
      if (j2 .le. m2) go to 101
      if (j1 .gt. m1) go to 109
  104 k = j-j1+1
      do 105 j=j1,m1
         m = k+j
         l = j+i1
         tcos(m) = tcos(l)
  105 continue
      go to 109
  106 continue
      if (j2 .gt. m2) go to 109
  107 k = j-j2+1
      do 108 j=j2,m2
         m = k+j
         l = j+i2
         tcos(m) = tcos(l)
  108 continue
  109 continue
      return
      end

      subroutine trix (idegbr,idegcr,m,a,b,c,y,tcos,d,w)
!
!     subroutine to solve a system of linear equations where the
!     coefficient matrix is a rational function in the matrix given by
!     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
!
      integer         idegbr, idegcr, m
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,y(*)       , &
     &                tcos(*)    ,d(*)       ,w(*)

!  local variables
      integer         i, k, l, ip, lint, mm1, ifb, ifc
      BIGREAL         x, xx, z

      mm1 = m-1
      ifb = idegbr+1
      ifc = idegcr+1
      l = ifb/ifc
      lint = 1
      do 108 k=1,idegbr
         x = tcos(k)
         if (k .ne. l) go to 102
         i = idegbr+lint
         xx = x-tcos(i)
         do 101 i=1,m
            w(i) = y(i)
            y(i) = xx*y(i)
  101    continue
  102    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y(1) = y(1)*z
         do 103 i=2,mm1
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
  103    continue
         z = b(m)-x-a(m)*d(mm1)
         if (z .ne. 0.) go to 104
         y(m) = 0.
         go to 105
  104    y(m) = (y(m)-a(m)*y(mm1))/z
  105    continue
         do 106 ip=1,mm1
            i = m-ip
            y(i) = y(i)-d(i)*y(i+1)
  106    continue
         if (k .ne. l) go to 108
         do 107 i=1,m
            y(i) = y(i)+w(i)
  107    continue
         lint = lint+1
         l = (lint*ifb)/ifc
  108 continue
      return
      end

      subroutine tri3 (m,a,b,c,k,y1,y2,y3,tcos,d,w1,w2,w3)
      integer         m
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,k(4)       , &
     &                tcos(*)    ,y1(*)      ,y2(*)      ,y3(*)      , &
     &                d(*)       ,w1(*)      ,w2(*)      ,w3(*)
!
!     subroutine to solve three linear systems whose common coefficient
!     matrix is a rational function in the matrix given by
!
!                  tridiagonal (...,a(i),b(i),c(i),...)
!
!  local variables
      integer          i, n, k1, k2, k3, k4, mm1, l1, l2, l3, k2k3k4, &
     &                 if1, if2, if3, if4, lint1, lint2, lint3, kint1, &
     &                 kint2, kint3, ip
      BIGREAL          x, xx, z

      mm1 = m-1
      k1 = k(1)
      k2 = k(2)
      k3 = k(3)
      k4 = k(4)
      if1 = k1+1
      if2 = k2+1
      if3 = k3+1
      if4 = k4+1
      k2k3k4 = k2+k3+k4
      if (k2k3k4 .eq. 0) go to 101
      l1 = if1/if2
      l2 = if1/if3
      l3 = if1/if4
      lint1 = 1
      lint2 = 1
      lint3 = 1
      kint1 = k1
      kint2 = kint1+k2
      kint3 = kint2+k3
  101 continue
      do 115 n=1,k1
         x = tcos(n)
         if (k2k3k4 .eq. 0) go to 107
         if (n .ne. l1) go to 103
         do 102 i=1,m
            w1(i) = y1(i)
  102    continue
  103    if (n .ne. l2) go to 105
         do 104 i=1,m
            w2(i) = y2(i)
  104    continue
  105    if (n .ne. l3) go to 107
         do 106 i=1,m
            w3(i) = y3(i)
  106    continue
  107    continue
         z = 1./(b(1)-x)
         d(1) = c(1)*z
         y1(1) = y1(1)*z
         y2(1) = y2(1)*z
         y3(1) = y3(1)*z
         do 108 i=2,m
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
  108    continue
         do 109 ip=1,mm1
            i = m-ip
            y1(i) = y1(i)-d(i)*y1(i+1)
            y2(i) = y2(i)-d(i)*y2(i+1)
            y3(i) = y3(i)-d(i)*y3(i+1)
  109    continue
         if (k2k3k4 .eq. 0) go to 115
         if (n .ne. l1) go to 111
         i = lint1+kint1
         xx = x-tcos(i)
         do 110 i=1,m
            y1(i) = xx*y1(i)+w1(i)
  110    continue
         lint1 = lint1+1
         l1 = (lint1*if1)/if2
  111    if (n .ne. l2) go to 113
         i = lint2+kint2
         xx = x-tcos(i)
         do 112 i=1,m
            y2(i) = xx*y2(i)+w2(i)
  112    continue
         lint2 = lint2+1
         l2 = (lint2*if1)/if3
  113    if (n .ne. l3) go to 115
         i = lint3+kint3
         xx = x-tcos(i)
         do 114 i=1,m
            y3(i) = xx*y3(i)+w3(i)
  114    continue
         lint3 = lint3+1
         l3 = (lint3*if1)/if4
  115 continue
      return
!
! revision history---
!
! september 1973    version 1
! april     1976    version 2
! january   1978    version 3
! december  1979    version 3.1
! october   1980    changed several divides of floating integers
!                   to integer divides to accomodate cray-1 arithmetic.
! february  1985    documentation upgrade
! november  1988    version 3.2, fortran 77 changes
! june      1993    BIGREAL stuff added
!-----------------------------------------------------------------------
      end

