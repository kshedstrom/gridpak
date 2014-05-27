
      subroutine blktri (iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y, &
     &                   ierror,w)
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
! dimension of           an(n),bn(n),cn(n),am(m),bm(m),cm(m),y(idimy,n),
! arguments              w(see argument list)
!
! latest revision        november 1988
!
! usage                  call blktri (iflg,np,n,an,bn,cn,mp,m,am,bm,
!                                     cm,idimy,y,ierror,w)
!
! purpose                blktri solves a system of linear equations
!                        of the form
!
!                        an(j)*x(i,j-1) + am(i)*x(i-1,j) +
!                        (bn(j)+bm(i))*x(i,j) + cn(j)*x(i,j+1) +
!                        cm(i)*x(i+1,j) = y(i,j)
!
!                        for i = 1,2,...,m  and  j = 1,2,...,n.
!
!                        i+1 and i-1 are evaluated modulo m and
!                        j+1 and j-1 modulo n, i.e.,
!
!                        x(i,0) = x(i,n),  x(i,n+1) = x(i,1),
!                        x(0,j) = x(m,j),  x(m+1,j) = x(1,j).
!
!                        these equations usually result from the
!                        discretization of separable elliptic
!                        equations.  boundary conditions may be
!                        dirichlet, neumann, or periodic.
!
! arguments
!
! on input               iflg
!
!                          = 0  initialization only.
!                               certain quantities that depend on np,
!                               n, an, bn, and cn are computed and
!                               stored in the work array w.
!
!                          = 1  the quantities that were computed
!                               in the initialization are used
!                               to obtain the solution x(i,j).
!
!                               note:
!                               a call with iflg=0 takes
!                               approximately one half the time
!                               as a call with iflg = 1.
!                               however, the initialization does
!                               not have to be repeated unless np,
!                               n, an, bn, or cn change.
!
!                        np
!                          = 0  if an(1) and cn(n) are not zero,
!                               which corresponds to periodic
!                               bounary conditions.
!
!                          = 1  if an(1) and cn(n) are zero.
!
!                        n
!                          the number of unknowns in the j-direction.
!                          n must be greater than 4.
!                          the operation count is proportional to
!                          mnlog2(n), hence n should be selected
!                          less than or equal to m.
!
!                        an,bn,cn
!                          one-dimensional arrays of length n
!                          that specify the coefficients in the
!                          linear equations given above.
!
!                        mp
!                          = 0  if am(1) and cm(m) are not zero,
!                               which corresponds to periodic
!                               boundary conditions.
!
!                          = 1  if am(1) = cm(m) = 0  .
!
!                        m
!                          the number of unknowns in the i-direction.
!                           m must be greater than 4.
!
!                        am,bm,cm
!                          one-dimensional arrays of length m that
!                          specify the coefficients in the linear
!                          equations given above.
!
!                        idimy
!                          the row (or first) dimension of the
!                          two-dimensional array y as it appears
!                          in the program calling blktri.
!                          this parameter is used to specify the
!                          variable dimension of y.
!                          idimy must be at least m.
!
!                        y
!                          a two-dimensional array that specifies
!                          the values of the right side of the linear
!                          system of equations given above.
!                          y must be dimensioned at least m*n.
!
!                        w
!                          a one-dimensional array that must be
!                          provided by the user for work space.
!                          if np=1 define k=int(log2(n))+1 and
!                          set l=2**(k+1) then w must have dimension
!                          (k-2)*l+k+5+max(2n,6m)
!
!                          if np=0 define k=int(log2(n-1))+1 and
!                          set l=2**(k+1) then w must have dimension
!                          (k-2)*l+k+5+2n+max(2n,6m)
!
!                          **important**
!                          for purposes of checking, the required
!                          dimension of w is computed by blktri and
!                          stored in w(1) in floating point format.
!
! arguments
!
! on output              y
!                          contains the solution x.
!
!                        ierror
!                          an error flag that indicates invalid
!                          input parameters.  except for number zer0,
!                          a solution is not attempted.
!
!                        = 0  no error.
!                        = 1  m is less than 5
!                        = 2  n is less than 5
!                        = 3  idimy is less than m.
!                        = 4  blktri failed while computing results
!                             that depend on the coefficient arrays
!                             an, bn, cn.  check these arrays.
!                        = 5  an(j)*cn(j-1) is less than 0 for some j.
!
!                             possible reasons for this condition are
!                             1. the arrays an and cn are not correct
!                             2. too large a grid spacing was used
!                                in the discretization of the elliptic
!                                equation.
!                             3. the linear equations resulted from a
!                                partial differential equation which
!                                was not elliptic.
!
!                        w
!                           contains intermediate values that must
!                           not be destroyed if blktri will be called
!                           again with iflg=1. w(1) contains the
!                           number of locations required by w in
!                           floating point format.
!
!
! special conditions     the algorithm may fail if abs(bm(i)+bn(j))
!                        is less than abs(am(i))+abs(an(j))+
!                        abs(cm(i))+abs(cn(j))
!                        for some i and j. the algorithm will also
!                        fail if an(j)*cn(j-1) is less than zero for
!                        some j.
!                        see the description of the output parameter
!                        ierror.
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
! history                written by paul swarztrauber at ncar in the
!                        early 1970's.  rewritten and released in
!                        january, 1980.
!
! algorithm              generalized cyclic reduction
!
! portability            fortran 77.  approximate machine accuracy
!                        is computed in function epmach.
!
! references             swarztrauber,p. and r. sweet, 'efficient
!                        fortran subprograms for the solution of
!                        elliptic equations'
!                        ncar tn/ia-109, july, 1975, 138 pp.
!
!                        swarztrauber p. n.,a direct method for
!                        the discrete solution of separable
!                        elliptic equations, s.i.a.m.
!                        j. numer. anal.,11(1974) pp. 1136-1150.
!***********************************************************************
#include "griddefs.h"
      integer         iflg       ,np         ,n          ,mp         ,  &
     &                m          ,idimy      ,ierror
      BIGREAL         an(*)      ,bn(*)      ,cn(*)      ,am(*)      ,  &
     &                bm(*)      ,cm(*)      ,y(idimy,1) ,w(*)
      external        prod       ,prodp      ,cprod      ,cprodp
      integer         npp        ,k          ,nm         ,ncmplx     ,  &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,  &
     &                nm         ,ncmplx     ,ik

!  local variables
      integer         nh         ,nl         ,iwah       ,iw1        ,  &
     &                iwbh       ,iw2        ,iw3        ,iwd        ,  &
     &                iww        ,iwu
!
! test m and n for the proper form
!
      nm = n
      ierror = 0
      if (m < 5) then
        ierror = 1; return
      endif
      if (nm < 3) then
        ierror = 2; return
      endif
      if (idimy < m) then
        ierror = 3; return
      endif
      nh = n
      npp = np
      if (npp==0) nh = nh+1
      ik = 2
      k = 1
  109 ik = ik+ik
      k = k+1
      if (nh > ik) goto 109
      nl = ik
      ik = ik+ik
      nl = nl-1
      iwah = (k-2)*ik+k+6
      if (npp==0) then
!
!     divide w into working sub arrays
!
        iwbh = iwah+nm+nm
        iw1 = iwbh
        w(1) = float(iw1-1+max0(2*nm,6*m))
        nm = nm-1
      else
        iw1 = iwah
        iwbh = iw1+nm
        w(1) = float(iw1-1+max0(2*nm,6*m))
      endif
!
! subroutine comp b computes the roots of the b polynomials
!
      if (ierror /= 0) return
      iw2 = iw1+m
      iw3 = iw2+m
      iwd = iw3+m
      iww = iwd+m
      iwu = iww+m
      if (iflg==0) then
        call compb (nl,ierror,an,bn,cn,w(2),w(iwah),w(iwbh))
        return
      endif
      if (mp==0) then
!
! subroutine blktr1 solves the linear system
!
        call blktr1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2), &
     &             w(iw3),w(iwd),w(iww),w(iwu),prodp,cprodp)
      else
        call blktr1 (nl,an,bn,cn,m,am,bm,cm,idimy,y,w(2),w(iw1),w(iw2), &
     &             w(iw3),w(iwd),w(iww),w(iwu),prod,cprod)
      endif
      return
      end
      subroutine blktr1 (n,an,bn,cn,m,am,bm,cm,idimy,y,b,w1,w2,w3,wd,   &
     &                   ww,wu,prdct,cprdct)
!
! blktr1 solves the linear system
!
! b  contains the roots of all the b polynomials
! w1,w2,w3,wd,ww,wu  are all working arrays
! prdct  is either prodp or prod depending on whether the boundary
! conditions in the m direction are periodic or not
! cprdct is either cprodp or cprod which are the complex versions
! of prodp and prod. these are called in the event that some
! of the roots of the b sub p polynomial are complex
!
!
      integer         n          ,m          ,idimy
      BIGREAL         an(*)      ,bn(*)      ,cn(*)      ,am(*)      ,  &
     &                bm(*)      ,cm(*)      ,b(*)       ,w1(*)      ,  &
     &                w2(*)      ,w3(*)      ,wd(*)      ,ww(*)      ,  &
     &                wu(*)      ,y(idimy,n)
      integer         npp        ,k          ,nm         ,ncmplx     ,  &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,  &
     &                nm         ,ncmplx     ,ik
      integer         kdo        ,l          ,i2         ,ir         ,  &
     &                im2        ,nm2        ,i1         ,irm1       ,  &
     &                im3        ,nm3        ,i3         ,im1        ,  &
     &                nm1        ,if         ,i          ,i4         ,  &
     &                idxa       ,idxc       ,nc         ,na         ,  &
     &                ipi1       ,ipi2       ,ipi3       ,ip1        ,  &
     &                ip2        ,ip3        ,np1        ,np2        ,  &
     &                np3        ,j          ,iz         ,nz         ,  &
     &                izr        ,ll         ,ifd        ,ip         ,  &
     &                np         ,imi1       ,imi2
      BIGREAL         dum
      external        prdct      ,cprdct
!
! begin reduction phase
!
      kdo = k-1
      do 109 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         call indxb (i2,ir,im2,nm2)
         call indxb (i1,irm1,im3,nm3)
         call indxb (i3,irm1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3, &
     &               m,am,bm,cm,wd,ww,wu)
         if = 2**k
         do 108 i=i4,if,i4
            if (i <= nm) then
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            call indxc (i,ir,idxc,nc)
            if (i < if) then
            call indxa (i,ir,idxa,na)
            call indxb (i-i1,irm1,im1,nm1)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am,  &
     &                  bm,cm,wd,ww,wu)
            if (ipi2 <= nm) then
              call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,       &
     &                  y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
              call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am,&
     &                  bm,cm,wd,ww,wu)
            else
              do j=1,m
                w3(j) = 0.
                w2(j) = 0.
              end do
            endif
            do j=1,m
               y(j,i) = w1(j)+w2(j)+y(j,i)
            end do
           end if
           end if
  108    continue
  109 continue
      if (npp) 132,110,132
!
!     the periodic case is treated using the capacitance matrix method
!
  110 if = 2**k
      i = if/2
      i1 = i/2
      call indxb (i-i1,k-2,im1,nm1)
      call indxb (i+i1,k-2,ip1,np1)
      call indxb (i,k-1,iz,nz)
      call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,y(1,i),w1,m,am, &
     &            bm,cm,wd,ww,wu)
      izr = i
      do 111 j=1,m
         w2(j) = w1(j)
  111 continue
      do 113 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i = i2
         call indxc (i,ir,idxc,nc)
         call indxb (i,ir,iz,nz)
         call indxb (i-i1,ir-1,im1,nm1)
         call indxb (i+i1,ir-1,ip1,np1)
         call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w1,w1,m,am,bm, &
     &               cm,wd,ww,wu)
         do 112 j=1,m
            w1(j) = y(j,i)+w1(j)
  112    continue
         call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,w1,m,am, &
     &               bm,cm,wd,ww,wu)
  113 continue
      do 118 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 117 i=i2,ifd,i4
            if (i-i2-izr) 117,114,117
  114       if (i-nm) 115,115,118
  115       call indxa (i,ir,idxa,na)
            call indxb (i,ir,iz,nz)
            call indxb (i-i1,ir-1,im1,nm1)
            call indxb (i+i1,ir-1,ip1,np1)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w2,w2,m,am, &
     &                  bm,cm,wd,ww,wu)
            do 116 j=1,m
               w2(j) = y(j,i)+w2(j)
  116       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w2,w2,m, &
     &                  am,bm,cm,wd,ww,wu)
            izr = i
            if (i-nm) 117,119,117
  117    continue
  118 continue
  119 do 120 j=1,m
         y(j,nm+1) = y(j,nm+1)-cn(nm+1)*w1(j)-an(nm+1)*w2(j)
  120 continue
      call indxb (if/2,k-1,im1,nm1)
      call indxb (if,k-1,ip,np)
      if (ncmplx) 121,122,121
  121 call cprdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
     &             y(1,nm+1),m,am,bm,cm,w1,w3,ww)
      go to 123
  122 call prdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1), &
     &            y(1,nm+1),m,am,bm,cm,wd,ww,wu)
  123 do 124 j=1,m
         w1(j) = an(1)*y(j,nm+1)
         w2(j) = cn(nm)*y(j,nm+1)
         y(j,1) = y(j,1)-w1(j)
         y(j,nm) = y(j,nm)-w2(j)
  124 continue
      do 126 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         i1 = i2/2
         i = i4
         call indxa (i,ir,idxa,na)
         call indxb (i-i2,ir,im2,nm2)
         call indxb (i-i2-i1,ir-1,im3,nm3)
         call indxb (i-i1,ir-1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,w1,w1,m,am, &
     &               bm,cm,wd,ww,wu)
         call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w1,w1,m,am,bm, &
     &               cm,wd,ww,wu)
         do 125 j=1,m
            y(j,i) = y(j,i)-w1(j)
  125    continue
  126 continue
!
      izr = nm
      do 131 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         do 130 i=i4,if,i4
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            if (ipi2-izr) 127,128,127
  127       if (i-izr) 130,131,130
  128       call indxc (i,ir,idxc,nc)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,w2,w2,m, &
     &                  am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w2,w2,m,am, &
     &                  bm,cm,wd,ww,wu)
            do 129 j=1,m
               y(j,i) = y(j,i)-w2(j)
  129       continue
            izr = i
            go to 131
  130    continue
  131 continue
!
! begin back substitution phase
!
  132 do 144 ll=1,k
         l = k-ll+1
         ir = l-1
         irm1 = ir-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 143 i=i2,ifd,i4
            if (i-nm) 133,133,143
  133       imi1 = i-i1
            imi2 = i-i2
            ipi1 = i+i1
            ipi2 = i+i2
            call indxa (i,ir,idxa,na)
            call indxc (i,ir,idxc,nc)
            call indxb (i,ir,iz,nz)
            call indxb (imi1,irm1,im1,nm1)
            call indxb (ipi1,irm1,ip1,np1)
            if (i-i2) 134,134,136
  134       do 135 j=1,m
               w1(j) = 0.
  135       continue
            go to 137
  136       call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2), &
     &                  w1,m,am,bm,cm,wd,ww,wu)
  137       if (ipi2-nm) 140,140,138
  138       do 139 j=1,m
               w2(j) = 0.
  139       continue
            go to 141
  140       call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2), &
     &                  w2,m,am,bm,cm,wd,ww,wu)
  141       do 142 j=1,m
               w1(j) = y(j,i)+w1(j)+w2(j)
  142       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i), &
     &                  m,am,bm,cm,wd,ww,wu)
  143    continue
  144 continue
      return
      end
      BIGREAL function bsrh (xll,xrr,iz,c,a,bh,f,sgn)
      integer         iz
      BIGREAL         a(*)       ,c(*)       ,bh(*)      ,xll        , &
     &                xrr        ,f          ,sgn
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      BIGREAL         xl         ,xr         ,dx         ,x
      xl = xll
      xr = xrr
      dx = .5*abs(xr-xl)
  101 x = .5*(xl+xr)
      if (sgn*f(x,iz,c,a,bh)) 103,105,102
  102 xr = x
      go to 104
  103 xl = x
  104 dx = .5*dx
      if (dx-cnv) 105,105,101
  105 bsrh = .5*(xl+xr)
      return
      end
      subroutine compb (n,ierror,an,bn,cn,b,ah,bh)
!
!     compb computes the roots of the b polynomials using subroutine
!     tevls which is a modification the eispack program tqlrat.
!     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
!     less than zero for some j.  ah,bh are temporary work arrays.
!
      integer         n          ,ierror
      BIGREAL         an(*)      ,bn(*)      ,cn(*)      ,b(*)       , &
     &                ah(*)      ,bh(*)
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      BIGREAL         dum        ,epmach     ,bnorm      ,tmp        , &
     &                arg        ,d1         ,d2         ,d3
      integer         j          ,if         ,kdo        ,l          , &
     &                i2         ,ir         ,ipl        ,i4         , &
     &                i          ,ifd        ,ib         ,nb         , &
     &                ls         ,js         ,jf         ,lh         , &
     &                nmp        ,l1         ,l2         ,j1         , &
     &                j2         ,n2m2
      eps = epmach(dum)
      bnorm = abs(bn(1))
      do 102 j=2,nm
         tmp = abs(bn(j))
         bnorm = max(bnorm,tmp)
         arg = an(j)*cn(j-1)
         if (arg) 119,101,101
  101    b(j) = sign(sqrt(arg),an(j))
  102 continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         ipl = i4-1
         ifd = if-i4
         do 107 i=i4,ifd,i4
            call indxb (i,l,ib,nb)
            if (nb) 108,108,103
  103       js = i-ipl
            jf = js+nb-1
            ls = 0
            do 104 j=js,jf
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
  104       continue
            call tevls (nb,bh,ah,ierror)
            if (ierror) 118,105,118
  105       lh = ib-1
            do 106 j=1,nb
               lh = lh+1
               b(lh) = -bh(j)
  106       continue
  107    continue
  108 continue
      do 109 j=1,nm
         b(j) = -bn(j)
  109 continue
      if (npp) 117,110,117
  110 nmp = nm+1
      nb = nm+nmp
      do 112 j=1,nb
         l1 = mod(j-1,nmp)+1
         l2 = mod(j+nm-1,nmp)+1
         arg = an(l1)*cn(l2)
         if (arg) 119,111,111
  111    bh(j) = sign(sqrt(arg),-an(l1))
         ah(j) = -bn(l1)
  112 continue
      call tevls (nb,ah,bh,ierror)
      if (ierror) 118,113,118
  113 call indxb (if,k-1,j2,lh)
      call indxb (if/2,k-1,j1,lh)
      j2 = j2+1
      lh = j2
      n2m2 = j2+nm+nm-2
  114 d1 = abs(b(j1)-b(j2-1))
      d2 = abs(b(j1)-b(j2))
      d3 = abs(b(j1)-b(j2+1))
      if ((d2 .lt. d1) .and. (d2 .lt. d3)) go to 115
      b(lh) = b(j2)
      j2 = j2+1
      lh = lh+1
      if (j2-n2m2) 114,114,116
  115 j2 = j2+1
      j1 = j1+1
      if (j2-n2m2) 114,114,116
  116 b(lh) = b(n2m2+1)
      call indxb (if,k-1,j1,j2)
      j2 = j1+nmp+nmp
#if DBLEPREC
      call ppadd (nm+1,ierror,an,cn,dcmplx(b(j1:j2)),b(j1),b(j2))
#else
      call ppadd (nm+1,ierror,an,cn,cmplx(b(j1:j2)),b(j1),b(j2))
#endif  /* DBLEPREC */
  117 return
  118 ierror = 4
      return
  119 ierror = 5
      return
      end
      subroutine cprod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,w,y)
!
! prod applies a sequence of matrix operations to the vector x and
! stores the result in yy           (complex case)
! aa   array containing scalar multipliers of the vector x
! nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
! bd,bm1,bm2 are arrays containing roots of certian b polynomials
! na is the length of the array aa
! x,yy the matrix operations are applied to x and the result is yy
! a,b,c  are arrays which contain the tridiagonal matrix
! m  is the order of the matrix
! d,w,y are working arrays
! isgn  determines whether or not a change in sign is made
!
      integer         nd         ,nm1        ,nm2        ,na         , &
     &                m
      BIGCOMPLEX      y(*)       ,d(*)       ,w(*)       ,bd(*)      , &
     &                crt        ,den        ,y1         ,y2
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,x(*)       , &
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
      integer         j          ,mm         ,id         ,m1         , &
     &                m2         ,ia         ,iflg       ,k
      BIGREAL         rt         ,zero
      zero = 0.
      do 101 j=1,m
         y(j) = cmplx(x(j),zero)
  101 continue
      mm = m-1
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 109,109,103
  103 crt = bd(id)
      id = id-1
!
! begin solution to system
!
      d(m) = a(m)/(b(m)-crt)
      w(m) = y(m)/(b(m)-crt)
      do 104 j=2,mm
         k = m-j
         den = b(k+1)-crt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  104 continue
      den = b(1)-crt-c(1)*d(2)
#if DBLEPREC
      if (cdabs(den)) 105,106,105
#else
      if (cabs(den)) 105,106,105
#endif  /* DBLEPREC */
  105 y(1) = (y(1)-c(1)*w(2))/den
      go to 107
  106 y(1) = (1.,0.)
  107 do 108 j=2,m
         y(j) = w(j)-d(j)*y(j-1)
  108 continue
  109 if (m1) 110,110,112
  110 if (m2) 121,121,111
  111 rt = bm2(m2)
      m2 = m2-1
      go to 117
  112 if (m2) 113,113,114
  113 rt = bm1(m1)
      m1 = m1-1
      go to 117
  114 if (abs(bm1(m1))-abs(bm2(m2))) 116,116,115
  115 rt = bm1(m1)
      m1 = m1-1
      go to 117
  116 rt = bm2(m2)
      m2 = m2-1
  117 y1 = (b(1)-rt)*y(1)+c(1)*y(2)
      if (mm-2) 120,118,118
!
! matrix multiplication
!
  118 do 119 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  119 continue
  120 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)
      y(m-1) = y1
      iflg = 1
      go to 102
  121 if (ia) 124,124,122
  122 rt = aa(ia)
      ia = ia-1
      iflg = 1
!
! scalar multiplication
!
      do 123 j=1,m
         y(j) = rt*y(j)
  123 continue
  124 if (iflg) 125,125,102
  125 do 126 j=1,m
#if DBLEPREC
         yy(j) = dble(y(j))
#else
         yy(j) = real(y(j))
#endif  /* DBLEPREC */
  126 continue
      return
      end
      subroutine cprodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,u,y)
!
! prodp applies a sequence of matrix operations to the vector x and
! stores the result in yy       periodic boundary conditions
! and  complex  case
!
! bd,bm1,bm2 are arrays containing roots of certian b polynomials
! nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
! aa   array containing scalar multipliers of the vector x
! na is the length of the array aa
! x,yy the matrix operations are applied to x and the result is yy
! a,b,c  are arrays which contain the tridiagonal matrix
! m  is the order of the matrix
! d,u,y are working arrays
! isgn  determines whether or not a change in sign is made
!
      integer         nd         ,nm1        ,nm2        ,na         , &
     &                m
      BIGCOMPLEX      y(*)       ,d(*)       ,u(*)       ,v          , &
     &                den        ,bh         ,ym         ,am         , &
     &                y1         ,y2         ,yh         ,bd(*)      , &
     &                crt
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,x(*)       , &
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,yy(*)
      integer         j          ,mm         ,mm2        ,id         , &
     &                m1         ,m2         ,ia         ,iflg
      integer         k
      BIGREAL         rt         ,zero
      zero = 0.
      do 101 j=1,m
         y(j) = cmplx(x(j),zero)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 111,111,103
  103 crt = bd(id)
      id = id-1
      iflg = 1
!
! begin solution to system
!
      bh = b(m)-crt
      ym = y(m)
      den = b(1)-crt
      d(1) = c(1)/den
      u(1) = a(1)/den
      y(1) = y(1)/den
#if DBLEPREC
      v = dcmplx(c(m),zero)
#else
      v = cmplx(c(m),zero)
#endif  /* DBLEPREC */
      if (mm2-2) 106,104,104
  104 do 105 j=2,mm2
         den = b(j)-crt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         y(j) = (y(j)-a(j)*y(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*y(j-1)
         v = -v*d(j-1)
  105 continue
  106 den = b(m-1)-crt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*y(m-2)
      den = bh-am*d(m-1)
#if DBLEPREC
      if (cdabs(den)) 107,108,107
#else
      if (cabs(den)) 107,108,107
#endif  /* DBLEPREC */
  107 y(m) = (ym-am*y(m-1))/den
      go to 109
  108 y(m) = (1.,0.)
  109 y(m-1) = y(m-1)-d(m-1)*y(m)
      do 110 j=2,mm
         k = m-j
         y(k) = y(k)-d(k)*y(k+1)-u(k)*y(m)
  110 continue
  111 if (m1) 112,112,114
  112 if (m2) 123,123,113
  113 rt = bm2(m2)
      m2 = m2-1
      go to 119
  114 if (m2) 115,115,116
  115 rt = bm1(m1)
      m1 = m1-1
      go to 119
  116 if (abs(bm1(m1))-abs(bm2(m2))) 118,118,117
  117 rt = bm1(m1)
      m1 = m1-1
      go to 119
  118 rt = bm2(m2)
      m2 = m2-1
!
! matrix multiplication
!
  119 yh = y(1)
      y1 = (b(1)-rt)*y(1)+c(1)*y(2)+a(1)*y(m)
      if (mm-2) 122,120,120
  120 do 121 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  121 continue
  122 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)+c(m)*yh
      y(m-1) = y1
      iflg = 1
      go to 102
  123 if (ia) 126,126,124
  124 rt = aa(ia)
      ia = ia-1
      iflg = 1
!
! scalar multiplication
!
      do 125 j=1,m
         y(j) = rt*y(j)
  125 continue
  126 if (iflg) 127,127,102
  127 do 128 j=1,m
#if DBLEPREC
         yy(j) = dble(y(j))
#else
         yy(j) = real(y(j))
#endif  /* DBLEPREC */
  128 continue
      return
      end
      subroutine indxa (i,ir,idxa,na)
      integer         i          ,ir         ,idxa       ,na
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
  101 na = 0
  102 return
      end
      subroutine indxb (i,ir,idx,idp)
!
! b(idx) is the location of the first root of the b(i,ir) polynomial
!
      integer         i          ,ir         ,idx        ,idp
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      integer         izh        ,ipl        ,id
      idx = MAX(i,1)
      idp = 0
      if (ir) 107,101,103
  101 if (i-nm) 102,102,107
  102 idx = i
      idp = 1
      return
  103 izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
  104 idp = 0
      return
  105 if (i+ipl-nm) 107,107,106
  106 idp = nm+ipl-i+1
  107 return
      end
      subroutine indxc (i,ir,idxc,nc)
      integer         i          ,ir         ,idxc       ,nc
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
      subroutine ppadd (n,ierror,a,c,cbp,bp,bh)
!
!     ppadd computes the eigenvalues of the periodic tridiagonal matrix
!     with coefficients an,bn,cn
!
! n is the order of the bh and bp polynomials
!     on output bp contians the eigenvalues
! cbp is the same as bp except type complex
! bh is used to temporarily store the roots of the b hat polynomial
! which enters through bp
!
      integer         n          ,ierror
      BIGCOMPLEX      cf         ,cx         ,fsg        ,hsg        , &
     &                dd         ,f          ,fp         ,fpp        , &
     &                cdis       ,r1         ,r2         ,r3         , &
     &                cbp(*)
      BIGREAL         a(*)       ,c(*)       ,bp(*)      ,bh(*)
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         eps        ,cnv
      common /cblkt/  npp        ,k          ,eps        ,cnv        , &
     &                nm         ,ncmplx     ,ik
      external        psgf       ,ppspf      ,ppsgf
      BIGREAL         scnv       ,xl         ,db         ,psgf       , &
     &                bsrh       ,xr         ,psg        ,xm
      integer         iz         ,j          ,nt         ,modiz      , &
     &                is         ,ig         ,if         ,it         , &
     &                icv        ,i2         ,i3         ,nhalf      , &
     &                izm        ,izm2
      BIGREAL         sgn        ,ppsgf      ,zero
      zero = 0.
      scnv = sqrt(cnv)
      iz = n
      izm = iz-1
      izm2 = iz-2
      if (bp(n)-bp(1)) 101,142,103
  101 do 102 j=1,n
         nt = n-j
         bh(j) = bp(nt+1)
  102 continue
      go to 105
  103 do 104 j=1,n
         bh(j) = bp(j)
  104 continue
  105 ncmplx = 0
      modiz = mod(iz,2)
      is = 1
      if (modiz) 106,107,106
  106 if (a(1)) 110,142,107
  107 xl = bh(1)
      db = bh(3)-bh(1)
  108 xl = xl-db
      if (psgf(xl,iz,c,a,bh)) 108,108,109
  109 sgn = -1.
      cbp(1) = cmplx(bsrh(xl,bh(1),iz,c,a,bh,psgf,sgn),zero)
      is = 2
  110 if = iz-1
      if (modiz) 111,112,111
  111 if (a(1)) 112,142,115
  112 xr = bh(iz)
      db = bh(iz)-bh(iz-2)
  113 xr = xr+db
      if (psgf(xr,iz,c,a,bh)) 113,114,114
  114 sgn = 1.
      cbp(iz) = cmplx(bsrh(bh(iz),xr,iz,c,a,bh,psgf,sgn),zero)
      if = iz-2
  115 do 136 ig=is,if,2
         xl = bh(ig)
         xr = bh(ig+1)
         sgn = -1.
         xm = bsrh(xl,xr,iz,c,a,bh,ppspf,sgn)
         psg = psgf(xm,iz,c,a,bh)
         if (abs(psg)-eps) 118,118,116
  116    if (psg*ppsgf(xm,iz,c,a,bh)) 117,118,119
!
!     case of a real zero
!
  117    sgn = 1.
         cbp(ig) = cmplx(bsrh(bh(ig),xm,iz,c,a,bh,psgf,sgn),zero)
         sgn = -1.
         cbp(ig+1) = cmplx(bsrh(xm,bh(ig+1),iz,c,a,bh,psgf,sgn),zero)
         go to 136
!
!     case of a multiple zero
!
  118    cbp(ig) = cmplx(xm,zero)
         cbp(ig+1) = cmplx(xm,zero)
         go to 136
!
!     case of a complex zero
!
  119    it = 0
         icv = 0
         cx = cmplx(xm,zero)
  120    fsg = (1.,0.)
         hsg = (1.,0.)
         fp = (0.,0.)
         fpp = (0.,0.)
         do 121 j=1,iz
            dd = 1./(cx-bh(j))
            fsg = fsg*a(j)*dd
            hsg = hsg*c(j)*dd
            fp = fp+dd
            fpp = fpp-dd*dd
  121    continue
         if (modiz) 123,122,123
  122    f = (1.,0.)-fsg-hsg
         go to 124
  123    f = (1.,0.)+fsg+hsg
  124    i3 = 0
#if DBLEPREC
         if (cdabs(fp)) 126,126,125
#else
         if (cabs(fp)) 126,126,125
#endif  /* DBLEPREC */
  125    i3 = 1
         r3 = -f/fp
  126    i2 = 0
#if DBLEPREC
         if (cdabs(fpp)) 132,132,127
#else
         if (cabs(fpp)) 132,132,127
#endif  /* DBLEPREC */
  127    i2 = 1
#if DBLEPREC
         cdis = cdsqrt(fp**2-2.*f*fpp)
#else
         cdis = csqrt(fp**2-2.*f*fpp)
#endif  /* DBLEPREC */
         r1 = cdis-fp
         r2 = -fp-cdis
#if DBLEPREC
         if (cdabs(r1)-cdabs(r2)) 129,129,128
#else
         if (cabs(r1)-cabs(r2)) 129,129,128
#endif  /* DBLEPREC */
  128    r1 = r1/fpp
         go to 130
  129    r1 = r2/fpp
  130    r2 = 2.*f/fpp/r1
#if DBLEPREC
         if (cdabs(r2) .lt. cdabs(r1)) r1 = r2
#else
         if (cabs(r2) .lt. cabs(r1)) r1 = r2
#endif  /* DBLEPREC */
         if (i3) 133,133,131
#if DBLEPREC
  131    if (cdabs(r3) .lt. cdabs(r1)) r1 = r3
#else
  131    if (cabs(r3) .lt. cabs(r1)) r1 = r3
#endif  /* DBLEPREC */
         go to 133
  132    r1 = r3
  133    cx = cx+r1
         it = it+1
         if (it .gt. 50) go to 142
#if DBLEPREC
         if (cdabs(r1) .gt. scnv) go to 120
#else
         if (cabs(r1) .gt. scnv) go to 120
#endif  /* DBLEPREC */
         if (icv) 134,134,135
  134    icv = 1
         go to 120
  135    cbp(ig) = cx
         cbp(ig+1) = conjg(cx)
  136 continue
#if DBLEPREC
      if (cdabs(cbp(n))-cdabs(cbp(1))) 137,142,139
#else
      if (cabs(cbp(n))-cabs(cbp(1))) 137,142,139
#endif  /* DBLEPREC */
  137 nhalf = n/2
      do 138 j=1,nhalf
         nt = n-j
         cx = cbp(j)
         cbp(j) = cbp(nt+1)
         cbp(nt+1) = cx
  138 continue
  139 ncmplx = 1
      do 140 j=2,iz
#if DBLEPREC
         if (dimag(cbp(j))) 143,140,143
#else
         if (aimag(cbp(j))) 143,140,143
#endif  /* DBLEPREC */
  140 continue
      ncmplx = 0
      do 141 j=2,iz
#if DBLEPREC
         bp(j) = dble(cbp(j))
#else
         bp(j) = real(cbp(j))
#endif  /* DBLEPREC */
  141 continue
      go to 143
  142 ierror = 4
  143 continue
      return
      end
      subroutine prod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u)
!
! prod applies a sequence of matrix operations to the vector x and
! stores the result in y
! bd,bm1,bm2 are arrays containing roots of certian b polynomials
! nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
! aa   array containing scalar multipliers of the vector x
! na is the length of the array aa
! x,y  the matrix operations are applied to x and the result is y
! a,b,c  are arrays which contain the tridiagonal matrix
! m  is the order of the matrix
! d,w,u are working arrays
! is  determines whether or not a change in sign is made
!
      integer         nd         ,nm1        ,nm2        ,na         , &
     &                m
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,x(*)       , &
     &                y(*)       ,d(*)       ,w(*)       ,bd(*)      , &
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,u(*)
      integer         j          ,mm         ,id         ,ibr        , &
     &                m1         ,m2         ,ia         ,k
      BIGREAL         rt         ,den
      do 101 j=1,m
         w(j) = x(j)
         y(j) = w(j)
  101 continue
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
!
! scalar multiplication
!
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 125,125,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
!
! begin solution to system
!
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do 107 j=2,mm
         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  107 continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = 1.
      if (den) 108,109,108
  108 w(1) = (y(1)-c(1)*w(2))/den
  109 do 110 j=2,m
         w(j) = w(j)-d(j)*w(j-1)
  110 continue
      if (na) 113,113,102
  111 do 112 j=1,m
         y(j) = w(j)
  112 continue
      ibr = 1
      go to 102
  113 if (m1) 114,114,115
  114 if (m2) 111,111,120
  115 if (m2) 117,117,116
  116 if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
  117 if (ibr) 118,118,119
  118 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
  119 rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
  120 if (ibr) 121,121,122
  121 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
  122 rt = rt-bm2(m2)
      m2 = m2-1
  123 do 124 j=1,m
         y(j) = y(j)+rt*w(j)
  124 continue
      go to 102
  125 return
      end
      subroutine prodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,u,w)
!
! prodp applies a sequence of matrix operations to the vector x and
! stores the result in y        periodic boundary conditions
!
! bd,bm1,bm2 are arrays containing roots of certian b polynomials
! nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
! aa   array containing scalar multipliers of the vector x
! na is the length of the array aa
! x,y  the matrix operations are applied to x and the result is y
! a,b,c  are arrays which contain the tridiagonal matrix
! m  is the order of the matrix
! d,u,w are working arrays
! is  determines whether or not a change in sign is made
!
      integer         nd         ,nm1        ,nm2        ,na         , &
     &                m
      BIGREAL         a(*)       ,b(*)       ,c(*)       ,x(*)       , &
     &                y(*)       ,d(*)       ,u(*)       ,bd(*)      , &
     &                bm1(*)     ,bm2(*)     ,aa(*)      ,w(*)
      integer         j          ,mm         ,mm2        ,id         , &
     &                ibr        ,m1         ,m2         ,ia         , &
     &                k
      BIGREAL         rt         ,ym         ,den        ,v          , &
     &                bh         ,am
      do 101 j=1,m
         y(j) = x(j)
         w(j) = y(j)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 128,128,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
!
! begin solution to system
!
      bh = b(m)-rt
      ym = y(m)
      den = b(1)-rt
      d(1) = c(1)/den
      u(1) = a(1)/den
      w(1) = y(1)/den
      v = c(m)
      if (mm2-2) 109,107,107
  107 do 108 j=2,mm2
         den = b(j)-rt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         w(j) = (y(j)-a(j)*w(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*w(j-1)
         v = -v*d(j-1)
  108 continue
  109 den = b(m-1)-rt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*w(m-2)
      den = bh-am*d(m-1)
      if (den) 110,111,110
  110 w(m) = (ym-am*w(m-1))/den
      go to 112
  111 w(m) = 1.
  112 w(m-1) = w(m-1)-d(m-1)*w(m)
      do 113 j=2,mm
         k = m-j
         w(k) = w(k)-d(k)*w(k+1)-u(k)*w(m)
  113 continue
      if (na) 116,116,102
  114 do 115 j=1,m
         y(j) = w(j)
  115 continue
      ibr = 1
      go to 102
  116 if (m1) 117,117,118
  117 if (m2) 114,114,123
  118 if (m2) 120,120,119
  119 if (abs(bm1(m1))-abs(bm2(m2))) 123,123,120
  120 if (ibr) 121,121,122
  121 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 114,122,122
  122 rt = rt-bm1(m1)
      m1 = m1-1
      go to 126
  123 if (ibr) 124,124,125
  124 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 114,125,125
  125 rt = rt-bm2(m2)
      m2 = m2-1
  126 do 127 j=1,m
         y(j) = y(j)+rt*w(j)
  127 continue
      go to 102
  128 return
      end
      subroutine tevls (n,d,e2,ierr)
!
      integer         i          ,j          ,l          ,m          , &
     &                n          ,ii         ,l1         ,mml        , &
     &                ierr
      BIGREAL         d(n)       ,e2(n)
      BIGREAL         b          ,c          ,f          ,g          , &
     &                h          ,p          ,r          ,s
!
!     real sqrt,abs,sign
!
      integer         npp        ,k          ,nm         ,ncmplx     , &
     &                ik
      BIGREAL         machep     ,cnv
      common /cblkt/  npp        ,k          ,machep     ,cnv        , &
     &                nm         ,ncmplx     ,ik
      integer         nhalf      ,ntop
      BIGREAL         dhold
!
!     this subroutine is a modification of the eispack subroutine tqlrat
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input-
!
!        n is the order of the matrix,
!
!        d contains the diagonal elements of the input matrix,
!
!        e2 contains the                subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output-
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues,
!
!        e2 has been destroyed,
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!
!
!     ********** machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!
!                **********
!
      ierr = 0
      if (n .eq. 1) go to 115
!
      do 101 i=2,n
         e2(i-1) = e2(i)*e2(i)
  101 continue
!
      f = 0.0
      b = 0.0
      e2(n) = 0.0
!
      do 112 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to 102
         b = h
         c = b*b
!
!     ********** look for small squared sub-diagonal element **********
!
  102    do 103 m=l,n
            if (e2(m) .le. c) go to 104
!
!     ********** e2(n) is always zero, so there is no exit
!                through the bottom of the loop **********
!
  103    continue
!
  104    if (m .eq. l) go to 108
  105    if (j .eq. 30) go to 114
         j = j+1
!
!     ********** form shift **********
!
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
!
         do 106 i=l1,n
            d(i) = d(i)-h
  106    continue
!
         f = f+h
!
!     ********** rational ql transformation **********
!
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
!
!     ********** for i=m-1 step -1 until l do -- **********
!
         do 107 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
  107    continue
!
         e2(l) = s*g
         d(l) = h
!
!     ********** guard against underflowed h **********
!
         if (h .eq. 0.0) go to 108
         if (abs(e2(l)) .le. abs(c/h)) go to 108
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to 105
  108    p = d(l)+f
!
!     ********** order eigenvalues **********
!
         if (l .eq. 1) go to 110
!
!     ********** for i=l step -1 until 2 do -- **********
!
         do 109 ii=2,l
            i = l+2-ii
            if (p .ge. d(i-1)) go to 111
            d(i) = d(i-1)
  109    continue
!
  110    i = 1
  111    d(i) = p
  112 continue
!
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
         ntop = n-i
         dhold = d(i)
         d(i) = d(ntop+1)
         d(ntop+1) = dhold
  113 continue
      go to 115
!
!     ********** set error -- no convergence to an
!                eigenvalue after 30 iterations **********
!
  114 ierr = l
  115 return
!
!     ********** last card of tqlrat **********
!
!
! revision history---
!
! september 1973    version 1
! april     1976    version 2
! january   1978    version 3
! december  1979    version 3.1
! february  1985    documentation upgrade
! november  1988    version 3.2, fortran 77 changes
!-----------------------------------------------------------------------
      end
