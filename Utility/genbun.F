      subroutine genbun (nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)
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
!
! dimension of           a(m),b(m),c(m),y(idimy,n),
!                        w(see parameter list)
! arguments
!
! latest revision        november 1988
!
! purpose                the name of this package is a mnemonic for the
!                        generalized buneman algorithm.
!
!                        it solves the real linear system of equations
!
!                        a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
!                        + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
!
!                        for i = 1,2,...,m  and  j = 1,2,...,n.
!
!                        indices i+1 and i-1 are evaluated modulo m,
!                        i.e., x(0,j) = x(m,j) and x(m+1,j) = x(1,j),
!                        and x(i,0) may equal 0, x(i,2), or x(i,n),
!                        and x(i,n+1) may equal 0, x(i,n-1), or x(i,1)
!                        depending on an input parameter.
!
! usage                  call genbun (nperod,n,mperod,m,a,b,c,idimy,y,
!                                     ierror,w)
!
! arguments
!
! on input               nperod
!
!                          indicates the values that x(i,0) and
!                          x(i,n+1) are assumed to have.
!
!                          = 0  if x(i,0) = x(i,n) and x(i,n+1) =
!                               x(i,1).
!                          = 1  if x(i,0) = x(i,n+1) = 0  .
!                          = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
!                          = 3  if x(i,0) = x(i,2) and x(i,n+1) =
!                               x(i,n-1).
!                          = 4  if x(i,0) = x(i,2) and x(i,n+1) = 0.
!
!                        n
!                          the number of unknowns in the j-direction.
!                          n must be greater than 2.
!
!                        mperod
!                          = 0 if a(1) and c(m) are not zero
!                          = 1 if a(1) = c(m) = 0
!
!                        m
!                          the number of unknowns in the i-direction.
!                          n must be greater than 2.
!
!                        a,b,c
!                          one-dimensional arrays of length m that
!                          specify the coefficients in the linear
!                          equations given above.  if mperod = 0
!                          the array elements must not depend upon
!                          the index i, but must be constant.
!                          specifically, the subroutine checks the
!                          following condition .
!
!                            a(i) = c(1)
!                            c(i) = c(1)
!                            b(i) = b(1)
!
!                          for i=1,2,...,m.
!
!                        idimy
!                          the row (or first) dimension of the
!                          two-dimensional array y as it appears
!                          in the program calling genbun.
!                          this parameter is used to specify the
!                          variable dimension of y.
!                          idimy must be at least m.
!
!                        y
!                          a two-dimensional complex array that
!                          specifies the values of the right side
!                          of the linear system of equations given
!                          above.
!                          y must be dimensioned at least m*n.
!
!                        w
!                          a one-dimensional array that must
!                          be provided by the user for work
!                          space.  w may require up to 4*n +
!                          (10 + int(log2(n)))*m locations.
!                          the actual number of locations used is
!                          computed by genbun and is returned in
!                          location w(1).
!
!
!  on output             y
!
!                          contains the solution x.
!
!                        ierror
!                          an error flag which indicates invalid
!                          input parameters  except for number
!                          zero, a solution is not attempted.
!
!                          = 0  no error.
!                          = 1  m .le. 2  .
!                          = 2  n .le. 2
!                          = 3  idimy .lt. m
!                          = 4  nperod .lt. 0 or nperod .gt. 4
!                          = 5  mperod .lt. 0 or mperod .gt. 1
!                          = 6  a(i) .ne. c(1) or c(i) .ne. c(1) or
!                               b(i) .ne. b(1) for
!                               some i=1,2,...,m.
!                          = 7  a(1) .ne. 0 or c(m) .ne. 0 and
!                                 mperod = 1
!
!                        w
!                          w(1) contains the required length of w.
!
! special conditons      none
!
! i/o                    none
!
! precision              single
!
! required library       comf and gnbnaux from fishpak
! files
!
! language               fortran
!
! history                written in 1979 by roland sweet of ncar's
!                        scientific computing division.  made available
!                        on ncar's public libraries in january, 1980.
!
! algorithm              the linear system is solved by a cyclic
!                        reduction algorithm described in the
!                        reference.
!
! portability            fortran 77 --
!                        the machine dependent constant pi is
!                        defined in function pimach.
!
! references             sweet, r., "a cyclic reduction algorithm for
!                        solving block tridiagonal systems of arbitrary
!                        dimensions," siam j. on numer. anal., 14 (1977)
!                        pp. 706-720.
!
! accuracy               this test was performed on a cdc 7600:
!
!                        a uniform random number generator was used
!                        to create a solution array x for the system
!                        given in the 'purpose' description above
!                        with
!                          a(i) = c(i) = -0.5*b(i) = 1, i=1,2,...,m
!
!                        and, when mperod = 1
!
!                          a(1) = c(m) = 0
!                          a(m) = c(1) = 2.
!
!                        the solution x was substituted into the
!                        given system  and, using double precision
!                        a right side y was computed.
!                        using this array y, subroutine genbun
!                        was called to produce approximate
!                        solution z.  then relative error
!                          e = max(abs(z(i,j)-x(i,j)))/
!                              max(abs(x(i,j)))
!                        was computed, where the two maxima are taken
!                        over i=1,2,...,m and j=1,...,n.
!
!                        the value of e is given in the table
!                        below for some typical values of m and n.
!
!                   m (=n)    mperod    nperod    t(msecs)    e
!                   ------    ------    ------    --------  ------
!
!                     31        0         0          36     6.e-14
!                     31        1         1          21     4.e-13
!                     31        1         3          41     3.e-13
!                     32        0         0          29     9.e-14
!                     32        1         1          32     3.e-13
!                     32        1         3          48     1.e-13
!                     33        0         0          36     9.e-14
!                     33        1         1          30     4.e-13
!                     33        1         3          34     1.e-13
!                     63        0         0         150     1.e-13
!                     63        1         1          91     1.e-12
!                     63        1         3         173     2.e-13
!                     64        0         0         122     1.e-13
!                     64        1         1         128     1.e-12
!                     64        1         3         199     6.e-13
!                     65        0         0         143     2.e-13
!                     65        1         1         120     1.e-12
!                     65        1         3         138     4.e-13
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#include "griddefs.h"
      integer       nperod, n, mperod, m, idimy, ierror
      BIGREAL       y(idimy,1)
      BIGREAL       w(*)       ,b(*)       ,a(*)       ,c(*)

!  local variables
      integer       i, j, k, np, mp, mh, mp1, iwba, iwbb, iwbc, &
     &              iwb2, iwb3, iww1, iww2, iww3, iwd, iwp, iwtcos, &
     &              ipstor, irev, modd, mhm1, mhmi, mhpi, nby2, &
     &              mskip
      BIGREAL       a1

!
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 102
      do 101 i=2,m
         if (a(i) .ne. c(1)) go to 103
         if (c(i) .ne. c(1)) go to 103
         if (b(i) .ne. b(1)) go to 103
  101 continue
      go to 104
  102 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
      go to 104
  103 ierror = 6
  104 if (ierror .ne. 0) return
      mp1 = m+1
      iwba = mp1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      mp = mperod+1
      np = nperod+1
      go to (114,107),mp
  107 go to (108,109,110,111,123),np
  108 call poisp2 (m,n,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2), &
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos), &
     &             w(iwp))
      go to 112
  109 call poisd2 (m,n,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iww1), &
     &             w(iwd),w(iwtcos),w(iwp))
      go to 112
  110 call poisn2 (m,n,1,2,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2), &
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos), &
     &             w(iwp))
      go to 112
  111 call poisn2 (m,n,1,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2), &
     &             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos), &
     &             w(iwp))
  112 ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 124
  113 go to (127,133),mp
  114 continue
!
!     reorder unknowns when mp =0
!
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 119 j=1,n
         do 115 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  115    continue
         w(mh) = 2.*y(mh,j)
         go to (117,116),modd
  116    w(m) = 2.*y(m,j)
  117    continue
         do 118 i=1,m
            y(i,j) = w(i)
  118    continue
  119 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.
      w(i) = 0.
      w(k+1) = 2.*w(k+1)
      go to (120,121),modd
  120 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 122
  121 w(iwbb-1) = w(k+1)
  122 continue
      go to 107
!
!     reverse columns when nperod = 4.
!
  123 irev = 1
      nby2 = n/2
  124 do 126 j=1,nby2
         mskip = n+1-j
         do 125 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  125    continue
  126 continue
      go to (110,113),irev
  127 continue
      do 132 j=1,n
         do 128 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  128    continue
         w(mh) = .5*y(mh,j)
         go to (130,129),modd
  129    w(m) = .5*y(m,j)
  130    continue
         do 131 i=1,m
            y(i,j) = w(i)
  131    continue
  132 continue
  133 continue
!
!     return storage requirements for w array.
!
      w(1) = ipstor+iwp-1
      return
      end

      subroutine poisd2 (mr,nr,istag,ba,bb,bc,q,idimq,b,w,d,tcos,p)
!
!     subroutine to solve poisson's equation for dirichlet boundary
!     conditions.
!
!     istag = 1 if the last diagonal block is the matrix a.
!     istag = 2 if the last diagonal block is the matrix a+i.
!
      integer         mr, nr, istag, idimq
      BIGREAL         q(idimq,1) ,ba(*)      ,bb(*)      ,bc(*)      , &
     &                tcos(*)    ,b(*)       ,d(*)       ,w(*)       , &
     &                p(*)
!  local variables
      integer         i, j, l, m, n, jsh, ip, ipstor, kr, irreg, &
     &                jstsav, lr, nun, jst, jsp, nodd, jp3, jm2, jp2, &
     &                jm1, jp1, noddpr, jm3, ip1, krpi, ideg, jdeg
      BIGREAL         fi, t
      external        merge
      BIGREAL         half, zero, one
      parameter    (  half = 0.5, zero = 0., one = 1.0 )

      m = mr
      n = nr
      jsh = 0
      fi = 1./float(istag)
      ip = -m
      ipstor = 0
      go to (101,102),istag
  101 kr = 0
      irreg = 1
      if (n .gt. 1) go to 106
      tcos(1) = 0.
      go to 103
  102 kr = 1
      jstsav = 1
      irreg = 2
      if (n .gt. 1) go to 106
      tcos(1) = -1.
  103 do 104 i=1,m
         b(i) = q(i,1)
  104 continue
      call trix (1,0,m,ba,bb,bc,b,tcos,d,w)
      do 105 i=1,m
         q(i,1) = b(i)
  105 continue
      go to 183
  106 lr = 0
      do 107 i=1,m
         p(i) = 0.
  107 continue
      nun = n
      jst = 1
      jsp = n
!
!     irreg = 1 when no irregularities have occurred, otherwise it is 2.
!
  108 l = 2*jst
      nodd = 2-2*((nun+1)/2)+nun
!
!     nodd = 1 when nun is odd, otherwise it is 2.
!
      go to (110,109),nodd
  109 jsp = jsp-l
      go to 111
  110 jsp = jsp-jst
      if (irreg .ne. 1) jsp = jsp-l
  111 continue
!
!     regular reduction
!
      call cosgen (jst,1,half,one,tcos)
      if (l .gt. jsp) go to 118
      do 117 j=l,jsp,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         jm3 = jm2-jsh
         jp3 = jp2+jsh
         if (jst .ne. 1) go to 113
         do 112 i=1,m
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  112    continue
         go to 115
  113    do 114 i=1,m
            t = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = t+q(i,j)-q(i,jm3)-q(i,jp3)
            q(i,j) = t
  114    continue
  115    continue
         call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
         do 116 i=1,m
            q(i,j) = q(i,j)+b(i)
  116    continue
  117 continue
!
!     reduction for last unknown
!
  118 go to (119,136),nodd
  119 go to (152,120),irreg
!
!     odd number of unknowns
!
  120 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (123,121),istag
  121 continue
      if (jst .ne. 1) go to 123
      do 122 i=1,m
         b(i) = q(i,j)
         q(i,j) = 0.
  122 continue
      go to 130
  123 go to (124,126),noddpr
  124 do 125 i=1,m
         ip1 = ip+i
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+p(ip1)+q(i,j)
  125 continue
      go to 128
  126 do 127 i=1,m
         b(i) = .5*(q(i,jm2)-q(i,jm1)-q(i,jm3))+q(i,jp2)-q(i,jp1)+q(i,j)
  127 continue
  128 do 129 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  129 continue
  130 call trix (jst,0,m,ba,bb,bc,b,tcos,d,w)
      ip = ip+m
      ipstor = max0(ipstor,ip+m)
      do 131 i=1,m
         ip1 = ip+i
         p(ip1) = q(i,j)+b(i)
         b(i) = q(i,jp2)+p(ip1)
  131 continue
      if (lr .ne. 0) go to 133
      do 132 i=1,jst
         krpi = kr+i
         tcos(krpi) = tcos(i)
  132 continue
      go to 134
  133 continue
      call cosgen (lr,jstsav,zero,fi,tcos(jst+1))
      call merge (tcos,0,jst,jst,lr,kr)
  134 continue
      call cosgen (kr,jstsav,zero,fi,tcos)
      call trix (kr,kr,m,ba,bb,bc,b,tcos,d,w)
      do 135 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+b(i)+p(ip1)
  135 continue
      lr = kr
      kr = kr+l
      go to 152
!
!     even number of unknowns
!
  136 jsp = jsp+l
      j = jsp
      jm1 = j-jsh
      jp1 = j+jsh
      jm2 = j-jst
      jp2 = j+jst
      jm3 = jm2-jsh
      go to (137,138),irreg
  137 continue
      jstsav = jst
      ideg = jst
      kr = l
      go to 139
  138 call cosgen (kr,jstsav,zero,fi,tcos)
      call cosgen (lr,jstsav,zero,fi,tcos(kr+1))
      ideg = kr
      kr = kr+jst
  139 if (jst .ne. 1) go to 141
      irreg = 2
      do 140 i=1,m
         b(i) = q(i,j)
         q(i,j) = q(i,jm2)
  140 continue
      go to 150
  141 do 142 i=1,m
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  142 continue
      go to (143,145),irreg
  143 do 144 i=1,m
         q(i,j) = q(i,jm2)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
  144 continue
      irreg = 2
      go to 150
  145 continue
      go to (146,148),noddpr
  146 do 147 i=1,m
         ip1 = ip+i
         q(i,j) = q(i,jm2)+p(ip1)
  147 continue
      ip = ip-m
      go to 150
  148 do 149 i=1,m
         q(i,j) = q(i,jm2)+q(i,j)-q(i,jm1)
  149 continue
  150 call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      do 151 i=1,m
         q(i,j) = q(i,j)+b(i)
  151 continue
  152 nun = nun/2
      noddpr = nodd
      jsh = jst
      jst = 2*jst
      if (nun .ge. 2) go to 108
!
!     start solution.
!
      j = jsp
      do 153 i=1,m
         b(i) = q(i,j)
  153 continue
      go to (154,155),irreg
  154 continue
      call cosgen (jst,1,half,zero,tcos)
      ideg = jst
      go to 156
  155 kr = lr+jst
      call cosgen (kr,jstsav,zero,fi,tcos)
      call cosgen (lr,jstsav,zero,fi,tcos(kr+1))
      ideg = kr
  156 continue
      call trix (ideg,lr,m,ba,bb,bc,b,tcos,d,w)
      jm1 = j-jsh
      jp1 = j+jsh
      go to (157,159),irreg
  157 do 158 i=1,m
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  158 continue
      go to 164
  159 go to (160,162),noddpr
  160 do 161 i=1,m
         ip1 = ip+i
         q(i,j) = p(ip1)+b(i)
  161 continue
      ip = ip-m
      go to 164
  162 do 163 i=1,m
         q(i,j) = q(i,j)-q(i,jm1)+b(i)
  163 continue
  164 continue
!
!     start back substitution.
!
      jst = jst/2
      jsh = jst/2
      nun = 2*nun
      if (nun .gt. n) go to 183
      do 182 j=jst,n,l
         jm1 = j-jsh
         jp1 = j+jsh
         jm2 = j-jst
         jp2 = j+jst
         if (j .gt. jst) go to 166
         do 165 i=1,m
            b(i) = q(i,j)+q(i,jp2)
  165    continue
         go to 170
  166    if (jp2 .le. n) go to 168
         do 167 i=1,m
            b(i) = q(i,j)+q(i,jm2)
  167    continue
         if (jst .lt. jstsav) irreg = 1
         go to (170,171),irreg
  168    do 169 i=1,m
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  169    continue
  170    continue
         call cosgen (jst,1,half,zero,tcos)
         ideg = jst
         jdeg = 0
         go to 172
  171    if (j+l .gt. n) lr = lr-jst
         kr = jst+lr
         call cosgen (kr,jstsav,zero,fi,tcos)
         call cosgen (lr,jstsav,zero,fi,tcos(kr+1))
         ideg = kr
         jdeg = lr
  172    continue
         call trix (ideg,jdeg,m,ba,bb,bc,b,tcos,d,w)
         if (jst .gt. 1) go to 174
         do 173 i=1,m
            q(i,j) = b(i)
  173    continue
         go to 182
  174    if (jp2 .gt. n) go to 177
  175    do 176 i=1,m
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
  176    continue
         go to 182
  177    go to (175,178),irreg
  178    if (j+jsh .gt. n) go to 180
         do 179 i=1,m
            ip1 = ip+i
            q(i,j) = b(i)+p(ip1)
  179    continue
         ip = ip-m
         go to 182
  180    do 181 i=1,m
            q(i,j) = b(i)+q(i,j)-q(i,jm1)
  181    continue
  182 continue
      l = l/2
      go to 164
  183 continue
!
!     return storage requirements for p vectors.
!
      w(1) = ipstor
      return
      end

      subroutine poisn2 (m,n,istag,mixbnd,a,bb,c,q,idimq,b,b2,b3,w,w2, &
     &                   w3,d,tcos,p)
!
!     subroutine to solve poisson's equation with neumann boundary
!     conditions.
!
!     istag = 1 if the last diagonal block is a.
!     istag = 2 if the last diagonal block is a-i.
!     mixbnd = 1 if have neumann boundary conditions at both boundaries.
!     mixbnd = 2 if have neumann boundary conditions at bottom and
!     dirichlet condition at top.  (for this case, must have istag = 1.)
!
      integer         m, n, istag, mixbnd, idimq
      BIGREAL         a(*)       ,bb(*)      ,c(*)       ,q(idimq,*) , &
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       , &
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    , &
     &                p(*)
!  local variables
      integer         i, j, k(4), k1, k2, k3, k4, mr, ipstor, i2r, jr, &
     &                nr, nlast, kr, lr, jstop, jp3, jp2, i2rby2, jm1, &
     &                jp1, jm2, jm3, nrodpr, ii, ip, i1, i2, jr2, &
     &                nlastp, jstep, nrod, jstart
      BIGREAL         fistag, fnum, fden, fi, t
      equivalence     (k(1),k1)  ,(k(2),k2)  ,(k(3),k3)  ,(k(4),k4)
      external        merge
      BIGREAL         half, zero, one
      parameter    (  half = 0.5, zero = 0., one = 1.0 )

      fistag = 3-istag
      fnum = 1./float(istag)
      fden = 0.5*float(istag-1)
      mr = m
      ip = -mr
      ipstor = 0
      i2r = 1
      jr = 2
      nr = n
      nlast = n
      kr = 1
      lr = 0
      go to (101,103),istag
  101 continue
      do 102 i=1,mr
         q(i,n) = .5*q(i,n)
  102 continue
      go to (103,104),mixbnd
  103 if (n .le. 3) go to 155
  104 continue
      jr = 2*i2r
      nrod = 1
      if ((nr/2)*2 .eq. nr) nrod = 0
      go to (105,106),mixbnd
  105 jstart = 1
      go to 107
  106 jstart = jr
      nrod = 1-nrod
  107 continue
      jstop = nlast-jr
      if (nrod .eq. 0) jstop = jstop-i2r
      call cosgen (i2r,1,half,zero,tcos)
      i2rby2 = i2r/2
      if (jstop .ge. jstart) go to 108
      j = jr
      go to 116
  108 continue
!
!     regular reduction.
!
      do 115 j=jstart,jstop,jr
         jp1 = j+i2rby2
         jp2 = j+i2r
         jp3 = jp2+i2rby2
         jm1 = j-i2rby2
         jm2 = j-i2r
         jm3 = jm2-i2rby2
         if (j .ne. 1) go to 109
         jm1 = jp1
         jm2 = jp2
         jm3 = jp3
  109    continue
         if (i2r .ne. 1) go to 111
         if (j .eq. 1) jm2 = jp2
         do 110 i=1,mr
            b(i) = 2.*q(i,j)
            q(i,j) = q(i,jm2)+q(i,jp2)
  110    continue
         go to 113
  111    continue
         do 112 i=1,mr
            fi = q(i,j)
            q(i,j) = q(i,j)-q(i,jm1)-q(i,jp1)+q(i,jm2)+q(i,jp2)
            b(i) = fi+q(i,j)-q(i,jm3)-q(i,jp3)
  112    continue
  113    continue
         call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
         do 114 i=1,mr
            q(i,j) = q(i,j)+b(i)
  114    continue
!
!     end of reduction for regular unknowns.
!
  115 continue
!
!     begin special reduction for last unknown.
!
      j = jstop+jr
  116 nlast = j
      jm1 = j-i2rby2
      jm2 = j-i2r
      jm3 = jm2-i2rby2
      if (nrod .eq. 0) go to 128
!
!     odd number of unknowns
!
      if (i2r .ne. 1) go to 118
      do 117 i=1,mr
         b(i) = fistag*q(i,j)
         q(i,j) = q(i,jm2)
  117 continue
      go to 126
  118 do 119 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  119 continue
      if (nrodpr .ne. 0) go to 121
      do 120 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)
  120 continue
      ip = ip-mr
      go to 123
  121 continue
      do 122 i=1,mr
         q(i,j) = q(i,j)-q(i,jm1)+q(i,jm2)
  122 continue
  123 if (lr .eq. 0) go to 124
      call cosgen (lr,1,half,fden,tcos(kr+1))
      go to 126
  124 continue
      do 125 i=1,mr
         b(i) = fistag*b(i)
  125 continue
  126 continue
      call cosgen (kr,1,half,fden,tcos)
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 127 i=1,mr
         q(i,j) = q(i,j)+b(i)
  127 continue
      kr = kr+i2r
      go to 151
  128 continue
!
!     even number of unknowns
!
      jp1 = j+i2rby2
      jp2 = j+i2r
      if (i2r .ne. 1) go to 135
      do 129 i=1,mr
         b(i) = q(i,j)
  129 continue
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      ip = 0
      ipstor = mr
      go to (133,130),istag
  130 do 131 i=1,mr
         p(i) = b(i)
         b(i) = b(i)+q(i,n)
  131 continue
      tcos(1) = 1.
      tcos(2) = 0.
      call trix (1,1,mr,a,bb,c,b,tcos,d,w)
      do 132 i=1,mr
         q(i,j) = q(i,jm2)+p(i)+b(i)
  132 continue
      go to 150
  133 continue
      do 134 i=1,mr
         p(i) = b(i)
         q(i,j) = q(i,jm2)+2.*q(i,jp2)+3.*b(i)
  134 continue
      go to 150
  135 continue
      do 136 i=1,mr
         b(i) = q(i,j)+.5*(q(i,jm2)-q(i,jm1)-q(i,jm3))
  136 continue
      if (nrodpr .ne. 0) go to 138
      do 137 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  137 continue
      go to 140
  138 continue
      do 139 i=1,mr
         b(i) = b(i)+q(i,jp2)-q(i,jp1)
  139 continue
  140 continue
      call trix (i2r,0,mr,a,bb,c,b,tcos,d,w)
      ip = ip+mr
      ipstor = max0(ipstor,ip+mr)
      do 141 i=1,mr
         ii = ip+i
         p(ii) = b(i)+.5*(q(i,j)-q(i,jm1)-q(i,jp1))
         b(i) = p(ii)+q(i,jp2)
  141 continue
      if (lr .eq. 0) go to 142
      call cosgen (lr,1,half,fden,tcos(i2r+1))
      call merge (tcos,0,i2r,i2r,lr,kr)
      go to 144
  142 do 143 i=1,i2r
         ii = kr+i
         tcos(ii) = tcos(i)
  143 continue
  144 call cosgen (kr,1,half,fden,tcos)
      if (lr .ne. 0) go to 145
      go to (146,145),istag
  145 continue
      call trix (kr,kr,mr,a,bb,c,b,tcos,d,w)
      go to 148
  146 continue
      do 147 i=1,mr
         b(i) = fistag*b(i)
  147 continue
  148 continue
      do 149 i=1,mr
         ii = ip+i
         q(i,j) = q(i,jm2)+p(ii)+b(i)
  149 continue
  150 continue
      lr = kr
      kr = kr+jr
  151 continue
      go to (152,153),mixbnd
  152 nr = (nlast-1)/jr+1
      if (nr .le. 3) go to 155
      go to 154
  153 nr = nlast/jr
      if (nr .le. 1) go to 192
  154 i2r = jr
      nrodpr = nrod
      go to 104
  155 continue
!
!      begin solution
!
      j = 1+jr
      jm1 = j-i2r
      jp1 = j+i2r
      jm2 = nlast-i2r
      if (nr .eq. 2) go to 184
      if (lr .ne. 0) go to 170
      if (n .ne. 3) go to 161
!
!     case n = 3.
!
      go to (156,168),istag
  156 continue
      do 157 i=1,mr
         b(i) = q(i,2)
  157 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 158 i=1,mr
         q(i,2) = b(i)
         b(i) = 4.*b(i)+q(i,1)+2.*q(i,3)
  158 continue
      tcos(1) = -2.
      tcos(2) = 2.
      i1 = 2
      i2 = 0
      call trix (i1,i2,mr,a,bb,c,b,tcos,d,w)
      do 159 i=1,mr
         q(i,2) = q(i,2)+b(i)
         b(i) = q(i,1)+2.*q(i,2)
  159 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 160 i=1,mr
         q(i,1) = b(i)
  160 continue
      jr = 1
      i2r = 0
      go to 194
!
!     case n = 2**p+1
!
  161 continue
      go to (162,170),istag
  162 continue
      do 163 i=1,mr
         b(i) = q(i,j)+.5*q(i,1)-q(i,jm1)+q(i,nlast)-q(i,jm2)
  163 continue
      call cosgen (jr,1,half,zero,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 164 i=1,mr
         q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))+b(i)
         b(i) = q(i,1)+2.*q(i,nlast)+4.*q(i,j)
  164 continue
      jr2 = 2*jr
      call cosgen (jr,1,zero,zero,tcos)
      do 165 i=1,jr
         i1 = jr+i
         i2 = jr+1-i
         tcos(i1) = -tcos(i2)
  165 continue
      call trix (jr2,0,mr,a,bb,c,b,tcos,d,w)
      do 166 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  166 continue
      call cosgen (jr,1,half,zero,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      do 167 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  167 continue
      go to 194
!
!     case of general n with nr = 3 .
!
  168 do 169 i=1,mr
         b(i) = q(i,2)
         q(i,2) = 0.
         b2(i) = q(i,3)
         b3(i) = q(i,1)
  169 continue
      jr = 1
      i2r = 0
      j = 2
      go to 177
  170 continue
      do 171 i=1,mr
         b(i) = .5*q(i,1)-q(i,jm1)+q(i,j)
  171 continue
      if (nrod .ne. 0) go to 173
      do 172 i=1,mr
         ii = ip+i
         b(i) = b(i)+p(ii)
  172 continue
      go to 175
  173 do 174 i=1,mr
         b(i) = b(i)+q(i,nlast)-q(i,jm2)
  174 continue
  175 continue
      do 176 i=1,mr
         t = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
         q(i,j) = t
         b2(i) = q(i,nlast)+t
         b3(i) = q(i,1)+2.*t
  176 continue
  177 continue
      k1 = kr+2*jr-1
      k2 = kr+jr
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (k2+istag-2,1,zero,fnum,tcos(k4))
      k4 = k1+k2+1
      call cosgen (jr-1,1,zero,one,tcos(k4))
      call merge (tcos,k1,k2,k1+k2,jr-1,0)
      k3 = k1+k2+lr
      call cosgen (jr,1,half,zero,tcos(k3+1))
      k4 = k3+jr+1
      call cosgen (kr,1,half,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,kr,k1)
      if (lr .eq. 0) go to 178
      call cosgen (lr,1,half,fden,tcos(k4))
      call merge (tcos,k3,jr,k3+jr,lr,k3-lr)
      call cosgen (kr,1,half,fden,tcos(k4))
  178 k3 = kr
      k4 = kr
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 179 i=1,mr
         b(i) = b(i)+b2(i)+b3(i)
  179 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 180 i=1,mr
         q(i,j) = q(i,j)+b(i)
         b(i) = q(i,1)+2.*q(i,j)
  180 continue
      call cosgen (jr,1,half,zero,tcos)
      call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
      if (jr .ne. 1) go to 182
      do 181 i=1,mr
         q(i,1) = b(i)
  181 continue
      go to 194
  182 continue
      do 183 i=1,mr
         q(i,1) = .5*q(i,1)-q(i,jm1)+b(i)
  183 continue
      go to 194
  184 continue
      if (n .ne. 2) go to 188
!
!     case  n = 2
!
      do 185 i=1,mr
         b(i) = q(i,1)
  185 continue
      tcos(1) = 0.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 186 i=1,mr
         q(i,1) = b(i)
         b(i) = 2.*(q(i,2)+b(i))*fistag
  186 continue
      tcos(1) = -fistag
      tcos(2) = 2.
      call trix (2,0,mr,a,bb,c,b,tcos,d,w)
      do 187 i=1,mr
         q(i,1) = q(i,1)+b(i)
  187 continue
      jr = 1
      i2r = 0
      go to 194
  188 continue
!
!     case of general n and nr = 2 .
!
      do 189 i=1,mr
         ii = ip+i
         b3(i) = 0.
         b(i) = q(i,1)+2.*p(ii)
         q(i,1) = .5*q(i,1)-q(i,jm1)
         b2(i) = 2.*(q(i,1)+q(i,nlast))
  189 continue
      k1 = kr+jr-1
      tcos(k1+1) = -2.
      k4 = k1+3-istag
      call cosgen (kr+istag-2,1,zero,fnum,tcos(k4))
      k4 = k1+kr+1
      call cosgen (jr-1,1,zero,one,tcos(k4))
      call merge (tcos,k1,kr,k1+kr,jr-1,0)
      call cosgen (kr,1,half,fden,tcos(k1+1))
      k2 = kr
      k4 = k1+k2+1
      call cosgen (lr,1,half,fden,tcos(k4))
      k3 = lr
      k4 = 0
      call tri3 (mr,a,bb,c,k,b,b2,b3,tcos,d,w,w2,w3)
      do 190 i=1,mr
         b(i) = b(i)+b2(i)
  190 continue
      tcos(1) = 2.
      call trix (1,0,mr,a,bb,c,b,tcos,d,w)
      do 191 i=1,mr
         q(i,1) = q(i,1)+b(i)
  191 continue
      go to 194
  192 do 193 i=1,mr
         b(i) = q(i,nlast)
  193 continue
      go to 196
  194 continue
!
!     start back substitution.
!
      j = nlast-jr
      do 195 i=1,mr
         b(i) = q(i,nlast)+q(i,j)
  195 continue
  196 jm2 = nlast-i2r
      if (jr .ne. 1) go to 198
      do 197 i=1,mr
         q(i,nlast) = 0.
  197 continue
      go to 202
  198 continue
      if (nrod .ne. 0) go to 200
      do 199 i=1,mr
         ii = ip+i
         q(i,nlast) = p(ii)
  199 continue
      ip = ip-mr
      go to 202
  200 do 201 i=1,mr
         q(i,nlast) = q(i,nlast)-q(i,jm2)
  201 continue
  202 continue
      call cosgen (kr,1,half,fden,tcos)
      call cosgen (lr,1,half,fden,tcos(kr+1))
      if (lr .ne. 0) go to 204
      do 203 i=1,mr
         b(i) = fistag*b(i)
  203 continue
  204 continue
      call trix (kr,lr,mr,a,bb,c,b,tcos,d,w)
      do 205 i=1,mr
         q(i,nlast) = q(i,nlast)+b(i)
  205 continue
      nlastp = nlast
  206 continue
      jstep = jr
      jr = i2r
      i2r = i2r/2
      if (jr .eq. 0) go to 222
      go to (207,208),mixbnd
  207 jstart = 1+jr
      go to 209
  208 jstart = jr
  209 continue
      kr = kr-jr
      if (nlast+jr .gt. n) go to 210
      kr = kr-jr
      nlast = nlast+jr
      jstop = nlast-jstep
      go to 211
  210 continue
      jstop = nlast-jr
  211 continue
      lr = kr-jr
      call cosgen (jr,1,half,zero,tcos)
      do 221 j=jstart,jstop,jstep
         jm2 = j-jr
         jp2 = j+jr
         if (j .ne. jr) go to 213
         do 212 i=1,mr
            b(i) = q(i,j)+q(i,jp2)
  212    continue
         go to 215
  213    continue
         do 214 i=1,mr
            b(i) = q(i,j)+q(i,jm2)+q(i,jp2)
  214    continue
  215    continue
         if (jr .ne. 1) go to 217
         do 216 i=1,mr
            q(i,j) = 0.
  216    continue
         go to 219
  217    continue
         jm1 = j-i2r
         jp1 = j+i2r
         do 218 i=1,mr
            q(i,j) = .5*(q(i,j)-q(i,jm1)-q(i,jp1))
  218    continue
  219    continue
         call trix (jr,0,mr,a,bb,c,b,tcos,d,w)
         do 220 i=1,mr
            q(i,j) = q(i,j)+b(i)
  220    continue
  221 continue
      nrod = 1
      if (nlast+i2r .le. n) nrod = 0
      if (nlastp .ne. nlast) go to 194
      go to 206
  222 continue
!
!     return storage requirements for p vectors.
!
      w(1) = ipstor
      return
      end

      subroutine poisp2 (m,n,a,bb,c,q,idimq,b,b2,b3,w,w2,w3,d,tcos,p)
!
!     subroutine to solve poisson equation with periodic boundary
!     conditions.
!
      integer         m, n, idimq
      BIGREAL         a(*)       ,bb(*)      ,c(*)       ,q(idimq,1) , &
     &                b(*)       ,b2(*)      ,b3(*)      ,w(*)       , &
     &                w2(*)      ,w3(*)      ,d(*)       ,tcos(*)    , &
     &                p(*)
!  local variables
      integer         i, j, mr, nr, nrm1, nrmj, nrpj, ipstor, lh
      BIGREAL         s, t

      mr = m
      nr = (n+1)/2
      nrm1 = nr-1
      if (2*nr .ne. n) go to 107
!
!     even number of unknowns
!
      do 102 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 101 i=1,mr
            s = q(i,nrmj)-q(i,nrpj)
            t = q(i,nrmj)+q(i,nrpj)
            q(i,nrmj) = s
            q(i,nrpj) = t
  101    continue
  102 continue
      do 103 i=1,mr
         q(i,nr) = 2.*q(i,nr)
         q(i,n) = 2.*q(i,n)
  103 continue
      call poisd2 (mr,nrm1,1,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr+1,1,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d, &
     &             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 105 j=1,nrm1
         nrmj = nr-j
         nrpj = nr+j
         do 104 i=1,mr
            s = .5*(q(i,nrpj)+q(i,nrmj))
            t = .5*(q(i,nrpj)-q(i,nrmj))
            q(i,nrmj) = s
            q(i,nrpj) = t
  104    continue
  105 continue
      do 106 i=1,mr
         q(i,nr) = .5*q(i,nr)
         q(i,n) = .5*q(i,n)
  106 continue
      go to 118
  107 continue
!
!     odd  number of unknowns
!
      do 109 j=1,nrm1
         nrpj = n+1-j
         do 108 i=1,mr
            s = q(i,j)-q(i,nrpj)
            t = q(i,j)+q(i,nrpj)
            q(i,j) = s
            q(i,nrpj) = t
  108    continue
  109 continue
      do 110 i=1,mr
         q(i,nr) = 2.*q(i,nr)
  110 continue
      lh = nrm1/2
      do 112 j=1,lh
         nrmj = nr-j
         do 111 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  111    continue
  112 continue
      call poisd2 (mr,nrm1,2,a,bb,c,q,idimq,b,w,d,tcos,p)
      ipstor = w(1)
      call poisn2 (mr,nr,2,1,a,bb,c,q(1,nr),idimq,b,b2,b3,w,w2,w3,d, &
     &             tcos,p)
      ipstor = max0(ipstor,int(w(1)))
      do 114 j=1,nrm1
         nrpj = nr+j
         do 113 i=1,mr
            s = .5*(q(i,nrpj)+q(i,j))
            t = .5*(q(i,nrpj)-q(i,j))
            q(i,nrpj) = t
            q(i,j) = s
  113    continue
  114 continue
      do 115 i=1,mr
         q(i,nr) = .5*q(i,nr)
  115 continue
      do 117 j=1,lh
         nrmj = nr-j
         do 116 i=1,mr
            s = q(i,j)
            q(i,j) = q(i,nrmj)
            q(i,nrmj) = s
  116    continue
  117 continue
  118 continue
!
!     return storage requirements for p vectors.
!
      w(1) = ipstor
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

