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
! package sepaux         contains no user entry points.
!
! latest revision        november 1988
!
! purpose                this package contains auxiliary routines for
!                        ncar public software packages such as sepeli
!                        and sepx4.
!
! usage                  since this package contains no user entries,
!                        no usage instructions or argument descriptions
!                        are given here.
!
! special conditions     none
!
! i/o                    none
!
! precision              single
!
! required library       none
! files
!
! language               fortran
!
! history                developed in the late 1970's by john c. adams
!                        of ncar's scienttific computing division.
!
! portability            fortran 77
! **********************************************************************
#include "griddefs.h"
      subroutine seport (usol,idmn,zn,zm,pertrb)
      integer         idmn
      BIGREAL         usol(idmn,1)    ,zn(*)      ,zm(*)   , pertrb
!
!     this subroutine orthoganalizes the array usol with respect to
!     the constant array in a weighted least squares norm
!
      integer         kswx       ,kswy       ,k          ,l          , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns
      BIGREAL         ait        ,bit        ,cit        ,dit        , &
     &                dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
      common /splp/   kswx       ,kswy       ,k          ,l          , &
     &                ait        ,bit        ,cit        ,dit        , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns         ,dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4

!  local variables
      integer         i, j, ii, jj, istr, ifnl, jstr, jfnl
      BIGREAL         ete, ute

      istr = is
      ifnl = ms
      jstr = js
      jfnl = ns
!
!     compute weighted inner products
!
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
!
!     set perturbation parameter
!
      pertrb = ute/ete
!
!     subtract off constant pertrb
!
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end

      subroutine sepmin (usol,idmn,zn,zm,pertb)
      integer         idmn
      BIGREAL         usol(idmn,1)   ,zn(*)      ,zm(*), pertb
!
!     this subroutine orhtogonalizes the array usol with respect to
!     the constant array in a weighted least squares norm
!
      integer         kswx       ,kswy       ,k          ,l          , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns
      BIGREAL         ait        ,bit        ,cit        ,dit        , &
     &                dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
      common /splp/   kswx       ,kswy       ,k          ,l          , &
     &                ait        ,bit        ,cit        ,dit        , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns         ,dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4

!  local variables
      integer         i, ii, j, jj, istr, ifnl, jstr, jfnl
      BIGREAL         ete, ute, pertrb
!
!     entry at sepmin occurrs when the final solution is
!     to be minimized with respect to the weighted
!     least squares norm
!
      istr = 1
      ifnl = k
      jstr = 1
      jfnl = l
!
!     compute weighted inner products
!
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
!
!     set perturbation parameter
!
      pertrb = ute/ete
!
!     subtract off constant pertrb
!
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end

      subroutine septri (n,a,b,c,d,u,z)
!
!     this subroutine solves for a non-zero eigenvector corresponding
!     to the zero eigenvalue of the transpose of the rank
!     deficient one matrix with subdiagonal a, diagonal b, and
!     superdiagonal c , with a(1) in the (1,n) position, with
!     c(n) in the (n,1) position, and all other elements zero.
!
      integer         n
      BIGREAL         a(n)       ,b(n)       ,c(n)       ,d(n)       , &
     &                u(n)       ,z(n)

!  local variables
      integer         j, k, nm1, nm2
      BIGREAL         bn, v, den, an

      bn = b(n)
      d(1) = a(2)/b(1)
      v = a(1)
      u(1) = c(n)/b(1)
      nm2 = n-2
      do  10 j=2,nm2
         den = b(j)-c(j-1)*d(j-1)
         d(j) = a(j+1)/den
         u(j) = -c(j-1)*u(j-1)/den
         bn = bn-v*u(j-1)
         v = -v*d(j-1)
   10 continue
      den = b(n-1)-c(n-2)*d(n-2)
      d(n-1) = (a(n)-c(n-2)*u(n-2))/den
      an = c(n-1)-v*d(n-2)
      bn = bn-v*u(n-2)
      den = bn-an*d(n-1)
!
!     set last component equal to one
!
      z(n) = 1.0
      z(n-1) = -d(n-1)
      nm1 = n-1
      do  20 j=2,nm1
         k = n-j
         z(k) = -d(k)*z(k+1)-u(k)*z(n)
   20 continue
      return
      end

      subroutine sepdx (u,idmn,i,j,uxxx,uxxxx)
      integer         idmn, i, j
      BIGREAL         u(idmn,1), uxxx, uxxxx
!
!     this program computes second order finite difference
!     approximations to the third and fourth x
!     partial derivatives of u at the (i,j) mesh point
!
      integer         kswx       ,kswy       ,k          ,l          , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns
      BIGREAL         ait        ,bit        ,cit        ,dit        , &
     &                dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
      common /splp/   kswx       ,kswy       ,k          ,l          , &
     &                ait        ,bit        ,cit        ,dit        , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns         ,dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4

      if (i.gt.2 .and. i.lt.(k-1)) go to  50
      if (i .eq. 1) go to  10
      if (i .eq. 2) go to  30
      if (i .eq. k-1) go to  60
      if (i .eq. k) go to  80
!
!     compute partial derivative approximations at x=a
!
   10 if (kswx .eq. 1) go to  20
      uxxx = (-5.0*u(1,j)+18.0*u(2,j)-24.0*u(3,j)+14.0*u(4,j)- &
     &                                               3.0*u(5,j))/(tdlx3)
      uxxxx = (3.0*u(1,j)-14.0*u(2,j)+26.0*u(3,j)-24.0*u(4,j)+ &
     &                                      11.0*u(5,j)-2.0*u(6,j))/dlx4
      return
!
!     periodic at x=a
!
   20 uxxx = (-u(k-2,j)+2.0*u(k-1,j)-2.0*u(2,j)+u(3,j))/(tdlx3)
      uxxxx = (u(k-2,j)-4.0*u(k-1,j)+6.0*u(1,j)-4.0*u(2,j)+u(3,j))/dlx4
      return
!
!     compute partial derivative approximations at x=a+dlx
!
   30 if (kswx .eq. 1) go to  40
      uxxx = (-3.0*u(1,j)+10.0*u(2,j)-12.0*u(3,j)+6.0*u(4,j)-u(5,j))/ &
     &       tdlx3
      uxxxx = (2.0*u(1,j)-9.0*u(2,j)+16.0*u(3,j)-14.0*u(4,j)+6.0*u(5,j)- &
     &                                                      u(6,j))/dlx4
      return
!
!     periodic at x=a+dlx
!
   40 uxxx = (-u(k-1,j)+2.0*u(1,j)-2.0*u(3,j)+u(4,j))/(tdlx3)
      uxxxx = (u(k-1,j)-4.0*u(1,j)+6.0*u(2,j)-4.0*u(3,j)+u(4,j))/dlx4
      return
!
!     compute partial derivative approximations on the interior
!
   50 continue
      uxxx = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
      uxxxx = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))/ &
     &        dlx4
      return
!
!     compute partial derivative approximations at x=b-dlx
!
   60 if (kswx .eq. 1) go to  70
      uxxx = (u(k-4,j)-6.0*u(k-3,j)+12.0*u(k-2,j)-10.0*u(k-1,j)+ &
     &                                                 3.0*u(k,j))/tdlx3
      uxxxx = (-u(k-5,j)+6.0*u(k-4,j)-14.0*u(k-3,j)+16.0*u(k-2,j)- &
     &                                     9.0*u(k-1,j)+2.0*u(k,j))/dlx4
      return
!
!     periodic at x=b-dlx
!
   70 uxxx = (-u(k-3,j)+2.0*u(k-2,j)-2.0*u(1,j)+u(2,j))/tdlx3
      uxxxx = (u(k-3,j)-4.0*u(k-2,j)+6.0*u(k-1,j)-4.0*u(1,j)+u(2,j))/ &
     &        dlx4
      return
!
!     compute partial derivative approximations at x=b
!
   80 uxxx = -(3.0*u(k-4,j)-14.0*u(k-3,j)+24.0*u(k-2,j)-18.0*u(k-1,j)+ &
     &                                                 5.0*u(k,j))/tdlx3
      uxxxx = (-2.0*u(k-5,j)+11.0*u(k-4,j)-24.0*u(k-3,j)+26.0*u(k-2,j)- &
     &                                    14.0*u(k-1,j)+3.0*u(k,j))/dlx4
      return
      end

      subroutine sepdy (u,idmn,i,j,uyyy,uyyyy)
      integer         idmn, i, j
      BIGREAL         u(idmn,*), uyyy, uyyyy
!
!     this program computes second order finite difference
!     approximations to the third and fourth y
!     partial derivatives of u at the (i,j) mesh point
!
      integer         kswx       ,kswy       ,k          ,l          , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns
      BIGREAL         ait        ,bit        ,cit        ,dit        , &
     &                dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
      common /splp/   kswx       ,kswy       ,k          ,l          , &
     &                ait        ,bit        ,cit        ,dit        , &
     &                mit        ,nit        ,is         ,ms         , &
     &                js         ,ns         ,dlx        ,dly        , &
     &                tdlx3      ,tdly3      ,dlx4       ,dly4
      if (j.gt.2 .and. j.lt.(l-1)) go to  50
      if (j .eq. 1) go to  10
      if (j .eq. 2) go to  30
      if (j .eq. l-1) go to  60
      if (j .eq. l) go to  80
!
!     compute partial derivative approximations at y=c
!
   10 if (kswy .eq. 1) go to  20
      uyyy = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)- &
     &                                                 3.0*u(i,5))/tdly3
      uyyyy = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+ &
     &                                      11.0*u(i,5)-2.0*u(i,6))/dly4
      return
!
!     periodic at x=a
!
   20 uyyy = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uyyyy = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))/dly4
      return
!
!     compute partial derivative approximations at y=c+dly
!
   30 if (kswy .eq. 1) go to  40
      uyyy = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))/ &
     &       tdly3
      uyyyy = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+6.0*u(i,5)- &
     &                                                      u(i,6))/dly4
      return
!
!     periodic at y=c+dly
!
   40 uyyy = (-u(i,l-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/tdly3
      uyyyy = (u(i,l-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
      return
!
!     compute partial derivative approximations on the interior
!
   50 continue
      uyyy = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uyyyy = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/ &
     &        dly4
      return
!
!     compute partial derivative approximations at y=d-dly
!
   60 if (kswy .eq. 1) go to  70
      uyyy = (u(i,l-4)-6.0*u(i,l-3)+12.0*u(i,l-2)-10.0*u(i,l-1)+ &
     &                                                 3.0*u(i,l))/tdly3
      uyyyy = (-u(i,l-5)+6.0*u(i,l-4)-14.0*u(i,l-3)+16.0*u(i,l-2)- &
     &                                     9.0*u(i,l-1)+2.0*u(i,l))/dly4
      return
!
!     periodic at y=d-dly
!
   70 continue
      uyyy = (-u(i,l-3)+2.0*u(i,l-2)-2.0*u(i,1)+u(i,2))/tdly3
      uyyyy = (u(i,l-3)-4.0*u(i,l-2)+6.0*u(i,l-1)-4.0*u(i,1)+u(i,2))/ &
     &        dly4
      return
!
!     compute partial derivative approximations at y=d
!
   80 uyyy = -(3.0*u(i,l-4)-14.0*u(i,l-3)+24.0*u(i,l-2)-18.0*u(i,l-1)+ &
     &                                                 5.0*u(i,l))/tdly3
      uyyyy = (-2.0*u(i,l-5)+11.0*u(i,l-4)-24.0*u(i,l-3)+26.0*u(i,l-2)- &
     &                                    14.0*u(i,l-1)+3.0*u(i,l))/dly4
      return
!
! revision history---
!
! september 1973    version 1
! april     1976    version 2
! january   1978    version 3
! december  1979    version 3.1
! february  1985    documentation upgrade
! november  1988    version 3.2, fortran 77 changes
! june      1993    BIGREAL stuff added
!-----------------------------------------------------------------------
      end
