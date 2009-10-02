      subroutine sepx4 (iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n, &
     &                  nbdcnd,bdc,bdd,cofx,grhs,usol,idmn,w,pertrb, &
     &                  ierror)
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
! dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
! arguments              usol(idmn,n+1),     grhs(idmn,n+1),
!                        w (see argument list)
!
! latest revision        november 1988
!
! purpose                sepx4 solves for either the second-order
!                        finite difference approximation or a
!                        fourth-order approximation to a separable
!                        elliptic equation
!
!                          af(x)*uxx+bf(x)*ux+cf(x)*u+uyy = g(x,y)
!
!                        on a rectangle (x greater than or equal to
!                        a and less than or equal to b, y greater than
!                        or equal to c and less than or equal to d).
!                        any combination of periodic or mixed boundary
!                        conditions is allowed.  if boundary
!                        conditions in the x direction are periodic
!                        (see mbdcnd=0 below) then the coefficients
!                        must satisfy
!
!                          af(x)=c1,bf(x)=0,cf(x)=c2 for all x.
!
!                        here c1,c2 are constants, c1.gt.0.
!
!                        the possible boundary conditions are:
!                        in the x-direction:
!                          (0) periodic, u(x+b-a,y)=u(x,y) for
!                              all y,x
!                          (1) u(a,y), u(b,y) are specified for all y
!                          (2) u(a,y), du(b,y)/dx+beta*u(b,y) are
!                              specified for all y
!                          (3) du(a,y)/dx+alpha*u(a,y),du(b,y)/dx+
!                              beta*u(b,y) are specified for all y
!                          (4) du(a,y)/dx+alpha*u(a,y),u(b,y) are
!                              specified for all y
!
!                        in the y-direction:
!                          (0) periodic, u(x,y+d-c)=u(x,y) for all x,y
!                          (1) u(x,c),u(x,d) are specified for all x
!                          (2) u(x,c),du(x,d)/dy are specified for
!                              all x
!                          (3) du(x,c)/dy,du(x,d)/dy are specified for
!                              all x
!                          (4) du(x,c)/dy,u(x,d) are specified for
!                              all x
!
! usage                  call sepx4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,
!                                   beta,c,d,n,nbdcnd,bdc,bdd,cofx,
!                                   grhs,usol,idmn,w,pertrb,ierror)
!
! arguments
! on input               iorder
!                          = 2 if a second-order approximation is
!                              sought
!                          = 4 if a fourth-order approximation is
!                              sought
!
!                        a,b
!                          the range of the x-independent variable,
!                          i.e., x is greater than or equal to a
!                          and less than or equal to b.  a must be
!                          less than b.
!
!                        m
!                          the number of panels into which the
!                          interval (a,b) is subdivided.  hence,
!                          there will be m+1 grid points in the x-
!                          direction given by xi=a+(i-1)*dlx
!                          for i=1,2,...,m+1 where dlx=(b-a)/m is
!                          the panel width.  m must be less than
!                          idmn and greater than 5.
!
!                        mbdcnd
!                          indicates the type of boundary condition
!                          at x=a and x=b
!                          = 0 if the solution is periodic in x, i.e.,
!                              u(x+b-a,y)=u(x,y)  for all y,x
!                          = 1 if the solution is specified at x=a
!                              and x=b, i.e., u(a,y) and u(b,y) are
!                              specified for all y
!                          = 2 if the solution is specified at x=a
!                              and the boundary condition is mixed at
!                              x=b, i.e., u(a,y) and
!                              du(b,y)/dx+beta*u(b,y) are specified
!                              for all y
!                          = 3 if the boundary conditions at x=a and
!                              x=b are mixed, i.e.,
!                              du(a,y)/dx+alpha*u(a,y) and
!                              du(b,y)/dx+beta*u(b,y) are specified
!                              for all y
!                          = 4 if the boundary condition at x=a is
!                              mixed and the solution is specified
!                              at x=b, i.e., du(a,y)/dx+alpha*u(a,y)
!                              and u(b,y) are specified for all y
!
!                        bda
!                          a one-dimensional array of length n+1 that
!                          specifies the values of
!                          du(a,y)/dx+ alpha*u(a,y) at x=a, when
!                          mbdcnd=3 or 4.
!                          bda(j) = du(a,yj)/dx+alpha*u(a,yj),
!                          j=1,2,...,n+1
!                          when mbdcnd has any other value, bda is
!                          a dummy parameter.
!
!                        alpha
!                          the scalar multiplying the solution in case
!                          of a mixed boundary condition at x=a
!                          (see argument bda).  if mbdcnd is not equal
!                          to either 3 or 4, then alpha is a dummy
!                          parameter.
!
!                        bdb
!                          a one-dimensional array of length n+1 that
!                          specifies the values of
!                          du(b,y)/dx+ beta*u(b,y) at x=b.
!                          when mbdcnd=2 or 3
!                          bdb(j) = du(b,yj)/dx+beta*u(b,yj),
!                          j=1,2,...,n+1
!                          when mbdcnd has any other value, bdb is
!                          a dummy parameter.
!
!                        beta
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at x=b
!                          (see argument bdb).  if mbdcnd is not equal
!                          to 2 or 3, then beta is a dummy parameter.
!
!                        c,d
!                          the range of the y-independent variable,
!                          i.e., y is greater than or equal to c and
!                          less than or equal to d.  c must be less
!                          than d.
!
!                        n
!                          the number of panels into which the
!                          interval (c,d) is subdivided.  hence,
!                          there will be n+1 grid points in the y-
!                          direction given by yj=c+(j-1)*dly for
!                          j=1,2,...,n+1 where dly=(d-c)/n is the
!                          panel width.  in addition, n must be
!                          greater than 4.
!
!                        nbdcnd
!                          indicates the types of boundary conditions
!                          at y=c and y=d
!                          = 0 if the solution is periodic in y,
!                              i.e., u(x,y+d-c)=u(x,y) for all x,y
!                          = 1 if the solution is specified at y=c
!                              and y = d, i.e., u(x,c)  and u(x,d)
!                              are specified for all x
!                          = 2 if the solution is specified at y=c
!                              and the boundary condition is mixed
!                              at y=d, i.e., du(x,c)/dy and u(x,d)
!                              are specified for all x
!                          = 3 if the boundary conditions are mixed
!                              at y=cand y=d i.e.,
!                              du(x,d)/dy and du(x,d)/dy are
!                              specified for all x
!                          = 4 if the boundary condition is mixed
!                              at y=c and the solution is specified
!                              at y=d, i.e. du(x,c)/dy+gama*u(x,c)
!                              and u(x,d) are specified for all x
!
!                        bdc
!                          a one-dimensional array of length m+1 that
!                          specifies the value du(x,c)/dy at y=c.
!
!                          when nbdcnd=3 or 4
!                            bdc(i) = du(xi,c)/dy i=1,2,...,m+1.
!
!                          when nbdcnd has any other value, bdc is
!                          a dummy parameter.
!
!                        bdd
!                          a one-dimensional array of length m+1 that
!                          specified the value of du(x,d)/dy at y=d.
!
!                          when nbdcnd=2 or 3
!                            bdd(i)=du(xi,d)/dy i=1,2,...,m+1.
!
!                          when nbdcnd has any other value, bdd is
!                          a dummy parameter.
!
!                        cofx
!                          a user-supplied subprogram with parameters
!                          x, afun, bfun, cfun which returns the
!                          values of the x-dependent coefficients
!                          af(x), bf(x), cf(x) in the elliptic
!                          equation at x.  if boundary conditions in
!                          the x direction are periodic then the
!                          coefficients must satisfy af(x)=c1,bf(x)=0,
!                          cf(x)=c2 for all x.  here c1.gt.0
!                          and c2 are constants.
!
!                          note that cofx must be declared external
!                          in the calling routine.
!
!                        grhs
!                          a two-dimensional array that specifies the
!                          values of the right-hand side of the
!                          elliptic equation, i.e.,grhs(i,j)=g(xi,yi),
!                          for i=2,...,m, j=2,...,n.  at the
!                          boundaries, grhs is defined by
!
!                          mbdcnd   grhs(1,j)   grhs(m+1,j)
!                          ------   ---------   -----------
!                            0      g(a,yj)     g(b,yj)
!                            1         *           *
!                            2         *        g(b,yj)  j=1,2,...,n+1
!                            3      g(a,yj)     g(b,yj)
!                            4      g(a,yj)        *
!
!                          nbdcnd   grhs(i,1)   grhs(i,n+1)
!                          ------   ---------   -----------
!                            0      g(xi,c)     g(xi,d)
!                            1         *           *
!                            2         *        g(xi,d)  i=1,2,...,m+1
!                            3      g(xi,c)     g(xi,d)
!                            4      g(xi,c)        *
!
!                          where * means these quantites are not used.
!                          grhs should be dimensioned idmn by at least
!                          n+1 in the calling routine.
!
!                        usol
!                          a two-dimensional array that specifies the
!                          values of the solution along the boundaries.
!                          at the boundaries, usol is defined by
!
!                          mbdcnd   usol(1,j)   usol(m+1,j)
!                          ------   ---------   -----------
!                            0         *           *
!                            1      u(a,yj)     u(b,yj)
!                            2      u(a,yj)        *     j=1,2,...,n+1
!                            3         *           *
!                            4         *        u(b,yj)
!
!                          nbdcnd   usol(i,1)   usol(i,n+1)
!                          ------   ---------   -----------
!                            0         *           *
!                            1      u(xi,c)     u(xi,d)
!                            2      u(xi,c)        *     i=1,2,...,m+1
!                            3         *           *
!                            4         *        u(xi,d)
!
!                          where * means the quantites are not used
!                          in the solution.
!
!                          if iorder=2, the user may equivalence grhs
!                          and usol to save space.  note that in this
!                          case the tables specifying the boundaries
!                          of the grhs and usol arrays determine the
!                          boundaries uniquely except at the corners.
!                          if the tables call for both g(x,y) and
!                          u(x,y) at a corner then the solution must
!                          be chosen.
!                          for example, if mbdcnd=2 and nbdcnd=4,
!                          then u(a,c), u(a,d),u(b,d) must be chosen
!                          at the corners in addition to g(b,c).
!
!                          if iorder=4, then the two arrays, usol and
!                          grhs, must be distinct.
!
!                          usol should be dimensioned idmn by at least
!                          n+1 in the calling routine.
!
!                        idmn
!                          the row (or first) dimension of the arrays
!                          grhs and usol as it appears in the program
!                          calling sepeli.  this parameter is used
!                          to specify the variable dimension of grhs
!                          and usol.  idmn must be at least 7 and
!                          greater than or equal to m+1.
!
!                        w
!                          a one-dimensional array that must be
!                          provided by the user for work space.
!                          10*n+(16+int(log2(n+1)))*(m+1)+11 will
!                          suffice as a length for w.  the actual
!                          length of w in the calling routine
!                          must be set in w(1) (see ierror=11).
!
!
! on output              usol
!                          contains the approximate solution to the
!                          elliptic equation. usol(i,j) is the
!                          approximation to u(xi,yj) for i=1,2...,m+1
!                          and j=1,2,...,n+1.  the approximation has
!                          error o(dlx**2+dly**2) if called with
!                          iorder=2 and o(dlx**4+dly**4) if called
!                          with iorder=4.
!
!                        w
!                          contains intermediate values that must not
!                          be destroyed if sepeli is called again
!                          with intl=1.  in addition w(1) contains
!                          the exact minimal length (in floating point)
!                          required for the work space (see ierror=11).
!
!                        pertrb
!                          if a combination of periodic or derivative
!                          boundary conditions (i.e., alpha=beta=0 if
!                          mbdcnd=3) is specified and if cf(x)=0 for
!                          all x then a solution to the discretized
!                          matrix equation may not exist
!                          (reflecting the non-uniqueness of solutions
!                          to the pde).
!                          pertrb is a constant calculated and
!                          subtracted from the right hand side of the
!                          matrix equation insuring the existence of a
!                          solution.  sepx4 computes this solution
!                          which is a weighted minimal least squares
!                          solution to the original problem.  if
!                          singularity is not detected pertrb=0.0 is
!                          returned by sepx4.
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters or failure to find a solution
!
!                          =  0 no error
!                          =  1 if a greater than b or c greater
!                               than d
!                          =  2 if mbdcnd less than 0 or mbdcnd
!                               greater than 4
!                          =  3 if nbdcnd less than 0 or nbdcnd
!                               greater than 4
!                          =  4 if attempt to find a solution fails.
!                               (the linear system generated is not
!                               diagonally dominant.)
!                          =  5 if idmn is too small (see discussion
!                               of idmn)
!                          =  6 if m is too small or too large
!                               (see discussion of m)
!                          =  7 if n is too small (see discussion of n)
!                          =  8 if iorder is not 2 or 4
!                          =  9 if intl is not 0 or 1
!                          = 10 if afun is less than or equal to zero
!                               for some interior mesh point xi some
!                               interior mesh point (xi,yj)
!                          = 11 if the work space length input in w(1)
!                               is less than the exact minimal work
!                               space length required output in w(1).
!                          = 12 if mbdcnd=0 and af(x)=cf(x)=constant
!                               or bf(x)=0 for all x is not true.
!
! special conditions     none
!
! i/o                    none
!
! required library       comf, genbun, gnbnaux, and sepaux
! files                  from fishpak
!
!
! precision              single
!
! required library       none
! files
!
! language               fortran
!
! history                sepx4 was developed at ncar by john c.
!                        adams of the scientific computing division
!                        in october 1978.  the basis of this code is
!                        ncar routine sepeli.  both packages were
!                        released on ncar's public libraries in
!                        january 1980.
!
! portability            fortran 77
!
! algorithm              sepx4 automatically discretizes the separable
!                        elliptic equation which is then solved by a
!                        generalized cyclic reduction algorithm in the
!                        subroutine pois.  the fourth order solution
!                        is obtained using the technique of defferred
!                        corrections referenced below.
!
! timing                 when possible, sepx4 should be used instead
!                        of package sepeli.  the increase in speed
!                        is at least a factor of three.
!
! references             keller, h.b., numerical methods for two-point
!                        boundary-value problems, blaisdel (1968),
!                        waltham, mass.
!
!                        swarztrauber, p., and r. sweet (1975):
!                        efficient fortran subprograms for the
!                        solution of elliptic partial differential
!                        equations.  ncar technical note
!                          ncar-tn/ia-109, pp. 135-137.
!***********************************************************************
#include "griddefs.h"
      integer         iorder, m, mbdcnd, n, nbdcnd, idmn, ierror
      BIGREAL         grhs(idmn,1)           ,usol(idmn,1)
      BIGREAL         bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     , &
     &                w(*)
      BIGREAL         a, b, alpha, beta, c, d, pertrb
      external        cofx

!  local variables
      integer         k, l, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, &
     &                i11, i12, i13, log2n, linput, loutpt, length

!
!     check input parameters
!
      call c4kprm(iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,idmn,ierror)
      if (ierror .ne. 0) return
!
!     compute minimum work space and check work space length input
!
      l = n+1
      if (nbdcnd .eq. 0) l = n
      k = m+1
      l = n+1
!     estimate log base 2 of n
      log2n=int(alog(float(n+1))/alog(2.0)+0.5)
      length=4*(n+1)+(10+log2n)*(m+1)
      ierror = 11
      linput = int(w(1)+0.5)
      loutpt = length+6*(k+l)+1
      w(1) = float(loutpt)
      if (loutpt .gt. linput) return
      ierror = 0
!
!     set work space indices
!
      i1 = length+2
      i2 = i1+l
      i3 = i2+l
      i4 = i3+l
      i5 = i4+l
      i6 = i5+l
      i7 = i6+l
      i8 = i7+k
      i9 = i8+k
      i10 = i9+k
      i11 = i10+k
      i12 = i11+k
      i13 = 2
      call s4elip(iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n, &
     &            nbdcnd,bdc,bdd,cofx,w(i1),w(i2),w(i3), &
     &            w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10),w(i11), &
     &            w(i12),grhs,usol,idmn,w(i13),pertrb,ierror)
      return
      end

      subroutine s4elip(iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n, &
     &                  nbdcnd,bdc,bdd,cofx,an,bn,cn,dn,un,zn,am,bm, &
     &                  cm,dm,um,zm,grhs,usol,idmn,w,pertrb,ierror)
!
!     s4elip sets up vectors and arrays for input to blktri
!     and computes a second order solution in usol.  a return jump to
!     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
!     solution is generated in usol.
!
      integer         iorder, m, mbdcnd, n, nbdcnd, idmn, ierror
      BIGREAL         bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     , &
     &                w(*)
      BIGREAL         grhs(idmn,1)           ,usol(idmn,1)
      BIGREAL         an(*)      ,bn(*)      ,cn(*)      ,dn(*)      , &
     &                un(*)      ,zn(*)
      BIGREAL         am(*)      ,bm(*)      ,cm(*)      ,dm(*)      , &
     &                um(*)      ,zm(*)
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
      logical         singlr
      external        cofx

!  local variables
      integer         i, j, i1, mp, np, iord, ieror
      BIGREAL         a, b, c, d, alpha, beta, pertrb, ai, bi, ci, &
     &                xi, axi, bxi, cxi, dyj, eyj, fyj, ax1, cxm, &
     &                dy1, fyn, gama, xnu, prtrb
!
!     set parameters internally
!
      kswx = mbdcnd+1
      kswy = nbdcnd+1
      k = m+1
      l = n+1
      ait = a
      bit = b
      cit = c
      dit = d
      dly=(dit-cit)/float(n)
!
!     set right hand side values from grhs in usol on the interior
!     and non-specified boundaries.
!
      do  20 i=2,m
         do  10 j=2,n
      usol(i,j)=dly**2*grhs(i,j)
   10    continue
   20 continue
      if (kswx.eq.2 .or. kswx.eq.3) go to  40
      do  30 j=2,n
      usol(1,j)=dly**2*grhs(1,j)
   30 continue
   40 continue
      if (kswx.eq.2 .or. kswx.eq.5) go to  60
      do  50 j=2,n
      usol(k,j)=dly**2*grhs(k,j)
   50 continue
   60 continue
      if (kswy.eq.2 .or. kswy.eq.3) go to  80
      do  70 i=2,m
      usol(i,1)=dly**2*grhs(i,1)
   70 continue
   80 continue
      if (kswy.eq.2 .or. kswy.eq.5) go to 100
      do  90 i=2,m
      usol(i,l)=dly**2*grhs(i,l)
   90 continue
  100 continue
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.3) &
     &usol(1,1)=dly**2*grhs(1,1)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.3) &
     &usol(k,1)=dly**2*grhs(k,1)
      if (kswx.ne.2 .and. kswx.ne.3 .and. kswy.ne.2 .and. kswy.ne.5) &
     &usol(1,l)=dly**2*grhs(1,l)
      if (kswx.ne.2 .and. kswx.ne.5 .and. kswy.ne.2 .and. kswy.ne.5) &
     &usol(k,l)=dly**2*grhs(k,l)
      i1 = 1
!
!     set switches for periodic or non-periodic boundaries
!
      mp=1
      if(kswx.eq.1) mp=0
      np=nbdcnd
!
!     set dlx,dly and size of block tri-diagonal system generated
!     in nint,mint
!
      dlx = (bit-ait)/float(m)
      mit = k-1
      if (kswx .eq. 2) mit = k-2
      if (kswx .eq. 4) mit = k
      dly = (dit-cit)/float(n)
      nit = l-1
      if (kswy .eq. 2) nit = l-2
      if (kswy .eq. 4) nit = l
      tdlx3 = 2.0*dlx**3
      dlx4 = dlx**4
      tdly3 = 2.0*dly**3
      dly4 = dly**4
!
!     set subscript limits for portion of array to input to blktri
!
      is = 1
      js = 1
      if (kswx.eq.2 .or. kswx.eq.3) is = 2
      if (kswy.eq.2 .or. kswy.eq.3) js = 2
      ns = nit+js-1
      ms = mit+is-1
!
!     set x - direction
!
      do 110 i=1,mit
         xi = ait+float(is+i-2)*dlx
         call cofx (xi,ai,bi,ci)
         axi = (ai/dlx-0.5*bi)/dlx
         bxi = -2.*ai/dlx**2+ci
         cxi = (ai/dlx+0.5*bi)/dlx
      am(i)=dly**2*axi
      bm(i)=dly**2*bxi
      cm(i)=dly**2*cxi
  110 continue
!
!     set y direction
!
      dyj=1.0
      eyj=-2.0
      fyj=1.0
      do 120 j=1,nit
         an(j) = dyj
         bn(j) = eyj
         cn(j) = fyj
  120 continue
!
!     adjust edges in x direction unless periodic
!
      ax1 = am(1)
      cxm = cm(mit)
      go to (170,130,150,160,140),kswx
!
!     dirichlet-dirichlet in x direction
!
  130 am(1) = 0.0
      cm(mit) = 0.0
      go to 170
!
!     mixed-dirichlet in x direction
!
  140 am(1) = 0.0
      bm(1) = bm(1)+2.*alpha*dlx*ax1
      cm(1) = cm(1)+ax1
      cm(mit) = 0.0
      go to 170
!
!     dirichlet-mixed in x direction
!
  150 am(1) = 0.0
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*beta*dlx*cxm
      cm(mit) = 0.0
      go to 170
!
!     mixed - mixed in x direction
!
  160 continue
      am(1) = 0.0
      bm(1) = bm(1)+2.*dlx*alpha*ax1
      cm(1) = cm(1)+ax1
      am(mit) = am(mit)+cxm
      bm(mit) = bm(mit)-2.*dlx*beta*cxm
      cm(mit) = 0.0
  170 continue
!
!     adjust in y direction unless periodic
!
      dy1 = an(1)
      fyn = cn(nit)
      gama=0.0
      xnu=0.0
      go to (220,180,200,210,190),kswy
!
!     dirichlet-dirichlet in y direction
!
  180 continue
      an(1) = 0.0
      cn(nit) = 0.0
      go to 220
!
!     mixed-dirichlet in y direction
!
  190 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      cn(nit) = 0.0
      go to 220
!
!     dirichlet-mixed in y direction
!
  200 an(1) = 0.0
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.*dly*xnu*fyn
      cn(nit) = 0.0
      go to 220
!
!     mixed - mixed direction in y direction
!
  210 continue
      an(1) = 0.0
      bn(1) = bn(1)+2.*dly*gama*dy1
      cn(1) = cn(1)+dy1
      an(nit) = an(nit)+fyn
      bn(nit) = bn(nit)-2.0*dly*xnu*fyn
      cn(nit) = 0.0
  220 if (kswx .eq. 1) go to 270
!
!     adjust usol along x edge
!
      do 260 j=js,ns
         if (kswx.ne.2 .and. kswx.ne.3) go to 230
         usol(is,j) = usol(is,j)-ax1*usol(1,j)
         go to 240
  230    usol(is,j) = usol(is,j)+2.0*dlx*ax1*bda(j)
  240    if (kswx.ne.2 .and. kswx.ne.5) go to 250
         usol(ms,j) = usol(ms,j)-cxm*usol(k,j)
         go to 260
  250    usol(ms,j) = usol(ms,j)-2.0*dlx*cxm*bdb(j)
  260 continue
  270 if (kswy .eq. 1) go to 320
!
!     adjust usol along y edge
!
      do 310 i=is,ms
         if (kswy.ne.2 .and. kswy.ne.3) go to 280
         usol(i,js) = usol(i,js)-dy1*usol(i,1)
         go to 290
  280    usol(i,js) = usol(i,js)+2.0*dly*dy1*bdc(i)
  290    if (kswy.ne.2 .and. kswy.ne.5) go to 300
         usol(i,ns) = usol(i,ns)-fyn*usol(i,l)
         go to 310
  300    usol(i,ns) = usol(i,ns)-2.0*dly*fyn*bdd(i)
  310 continue
  320 continue
!
!     save adjusted edges in grhs if iorder=4
!
      if (iorder .ne. 4) go to 350
      do 330 j=js,ns
         grhs(is,j) = usol(is,j)
         grhs(ms,j) = usol(ms,j)
  330 continue
      do 340 i=is,ms
         grhs(i,js) = usol(i,js)
         grhs(i,ns) = usol(i,ns)
  340 continue
  350 continue
      iord = iorder
      pertrb = 0.0
!
!     check if operator is singular
!
      call c4ksng(mbdcnd,nbdcnd,alpha,beta,cofx,singlr)
!
!     compute non-zero eigenvector in null space of transpose
!     if singular
!
      if (singlr) call septri (mit,am,bm,cm,dm,um,zm)
      if (singlr) call septri (nit,an,bn,cn,dn,un,zn)
!
!     adjust right hand side if necessary
!
  360 continue
      if (singlr) call seport (usol,idmn,zn,zm,pertrb)
!
!     compute solution
!
!     save adjusted right hand side in grhs
      do 444 j=js,ns
      do 444 i=is,ms
      grhs(i,j)=usol(i,j)
  444 continue
      call genbun(np,nit,mp,mit,am,bm,cm,idmn,usol(is,js),ieror,w)
!     check if error detected in pois
!     this can only correspond to ierror=12
      if(ieror.eq.0) go to 224
!     set error flag if improper coefficients input to pois
      ierror=12
      return
  224 continue
      if (ierror .ne. 0) return
!
!     set periodic boundaries if necessary
!
      if (kswx .ne. 1) go to 380
      do 370 j=1,l
         usol(k,j) = usol(1,j)
  370 continue
  380 if (kswy .ne. 1) go to 400
      do 390 i=1,k
         usol(i,l) = usol(i,1)
  390 continue
  400 continue
!
!     minimize solution with respect to weighted least squares
!     norm if operator is singular
!
      if (singlr) call sepmin (usol,idmn,zn,zm,prtrb)
!
!     return if deferred corrections and a fourth order solution are
!     not flagged
!
      if (iord .eq. 2) return
      iord = 2
!
!     compute new right hand side for fourth order solution
!
      call d4fer(cofx,idmn,usol,grhs)
      go to 360
      end

      subroutine c4kprm(iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,idmn, &
     &   ierror)
      integer     iorder, m, mbdcnd, n, nbdcnd, idmn, ierror
      BIGREAL     a, b, c, d
      external    cofx

!  local variables
      integer     i
      BIGREAL     ai, bi, ci, xi, dlx
!
!     this program checks the input parameters for errors
!
!
!
!     check definition of solution region
!
      ierror = 1
      if (a.ge.b .or. c.ge.d) return
!
!     check boundary switches
!
      ierror = 2
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) return
      ierror = 3
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) return
!
!     check first dimension in calling routine
!
      ierror = 5
      if (idmn .lt. 7) return
!
!     check m
!
      ierror = 6
      if (m.gt.(idmn-1) .or. m.lt.6) return
!
!     check n
!
      ierror = 7
      if (n .lt. 5) return
!
!     check iorder
!
      ierror = 8
      if (iorder.ne.2 .and. iorder.ne.4) return
!
!     check intl
!
!
!     check that equation is elliptic
!
      dlx = (b-a)/float(m)
      do  30 i=2,m
         xi = a+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
      if (ai.gt.0.0) go to 10
      ierror=10
      return
   10 continue
   30 continue
!
!     no error found
!
      ierror = 0
      return
      end

      subroutine c4ksng(mbdcnd,nbdcnd,alpha,beta,cofx,singlr)
      integer         mbdcnd, nbdcnd
      BIGREAL         alpha, beta
      external        cofx
      logical         singlr
!
!     this subroutine checks if the pde   sepeli
!     must solve is a singular operator
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
      integer         i
      BIGREAL         ai, bi, ci, xi

      singlr = .false.
!
!     check if the boundary conditions are
!     entirely periodic and/or mixed
!
      if ((mbdcnd.ne.0 .and. mbdcnd.ne.3) .or. &
     &    (nbdcnd.ne.0 .and. nbdcnd.ne.3)) return
!
!     check that mixed conditions are pure neuman
!
      if (mbdcnd .ne. 3) go to  10
      if (alpha.ne.0.0 .or. beta.ne.0.0) return
   10 continue
!
!     check that non-derivative coefficient functions
!     are zero
!
      do  30 i=is,ms
         xi = ait+float(i-1)*dlx
         call cofx (xi,ai,bi,ci)
         if (ci .ne. 0.0) return
   30 continue
!
!     the operator must be singular if this point is reached
!
      singlr = .true.
      return
      end

      subroutine d4fer(cofx,idmn,usol,grhs)
      integer         idmn
      BIGREAL         grhs(idmn,1)    ,usol(idmn,1)
      external        cofx
!
!     this subroutine first approximates the truncation error given by
!     trun1(x,y)=dlx**2*tx+dly**2*ty where
!     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
!     at the boundaries if periodic(here uxxx,uxxxx are the third
!     and fourth partial derivatives of u with respect to x).
!     tx is of the form afun(x)/3.0*(uxxxx/4.0+uxxx/dlx)
!     at x=a or x=b if the boundary condition there is mixed.
!     tx=0.0 along specified boundaries.  ty has symmetric form
!     in y with x,afun(x),bfun(x) replaced by y,dfun(y),efun(y).
!     the second order solution in usol is used to approximate
!     (via second order finite differencing) the truncation error
!     and the result is added to the right hand side in grhs
!     and then transferred to usol to be used as a new right
!     hand side when calling blktri for a fourth order solution.
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
      integer         i, j
      BIGREAL         ai, bi, ci, xi, uxxx, uxxxx, uyyy, uyyyy, &
     &                tx, ty
!
!
!     compute truncation error approximation over the entire mesh
!
         do  30 i=is,ms
            xi = ait+float(i-1)*dlx
            call cofx (xi,ai,bi,ci)
         do 30 j=js,ns
!
!     compute partial derivative approximations at (xi,yj)
!
            call sepdx (usol,idmn,i,j,uxxx,uxxxx)
            call sepdy (usol,idmn,i,j,uyyy,uyyyy)
            tx = ai*uxxxx/12.0+bi*uxxx/6.0
             ty=uyyyy/12.0
!
!     reset form of truncation if at boundary which is non-periodic
!
            if (kswx.eq.1 .or. (i.gt.1 .and. i.lt.k)) go to  10
            tx = ai/3.0*(uxxxx/4.0+uxxx/dlx)
   10       if (kswy.eq.1 .or. (j.gt.1 .and. j.lt.l)) go to  20
            ty = (uyyyy/4.0+uyyy/dly)/3.0
   20 grhs(i,j)=grhs(i,j)+dly**2*(dlx**2*tx+dly**2*ty)
   30    continue
!
!     reset the right hand side in usol
!
      do  60 i=is,ms
         do  50 j=js,ns
            usol(i,j) = grhs(i,j)
   50    continue
   60 continue
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

