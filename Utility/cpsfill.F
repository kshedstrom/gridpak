!  *** In gridpak version 5.4 ***** October 18, 2001 ****************
!  Kate Hedstrom  (kate@arsc.edu)
!  John Wilkin  (wilkin@imcs.rutgers.edu)
!
!  emulate our old modified conrec - with some hints from cpcnrc
!
!  modification:  John Wilkin, 16 December, 1991
!  Enabled hardwiring of contour levels for use with plothist
! *******************************************************************

      subroutine cpshift(zdat,K,M,N,finc,shift,scalel)
      logical      shift
      integer      K, M, N
#include "griddefs.h"
      real         zdat(K,*), finc, scalel

!  dimensions of work arrays
      integer       lrwk, liwk, lama
      parameter (   lrwk=5000, liwk=2000, lama=24000   )
      integer       iwrk(liwk), iama(lama)
      real          rwrk(lrwk), ciu, clev, spval

! Declare the contour-line drawing routine
      external cpdrpl

      logical          colour, tallflg, mapflag
      integer          lcflag, ncontr
      common / cpflg / colour, tallflg, lcflag, ncontr, mapflag
      integer          i, ncl, iclu, llpf

! Extra parameters if we decide to hardwire the contour levels
      real             cpmin, cpmax, clevel
      integer          nclevels
      common / cfix  / cpmin, cpmax, nclevels

! if colour is true we actually want colour fill plots so
! divert to cpsfill

      if (colour) then
        call cpsfill(zdat,K,M,N,ncontr,lcflag,tallflg,.false.)
        return
      endif

      spval = 1.e+20

!  turn on drawing of thick edge of plot
      call cpseti('PAI - array index',-1)
      call cpseti('CLU - contour level use',1)
      call cpsetr('CLL - contour line width',2.)

!  turn on drawing of thick edge of special values
      call cpseti('PAI - array index',-2)
      call cpseti('CLU - contour level use',1)
      call cpsetr('CLL - contour line width',2.)

!  arrange for selection of contour levels

      if (nclevels .eq. 0) then

!  default scheme -- cpshift follows contour level selection hints
!  from input argument list

        if (finc .lt. 0.) then
          call cpseti('CLS - contour level selector',max(1,int(-finc)))
          call cpsetr('CIS - contour interval selector',0.)
        else if (finc .eq. 0.) then
          call cpseti('CLS - contour level selector',16)
          call cpsetr('CIS - contour interval selector',0.)
        else
          call cpseti('CLS - contour level selector',1)
          call cpsetr('CIS - contour interval selector',finc)
        end if

      else

!  override selection hints
!  cpmin, cpmax and nclevels are input to plothist

        call cpseti('CLS - contour level selector',0)
        call cpseti('NCL - number of contour levels',nclevels+1)
        do i=1,nclevels+1
          clevel = cpmin + float(i-1)/float(nclevels)*(cpmax-cpmin)
          call cpseti('PAI - array index',i)
          call cpsetr('CLV - contour level',clevel)
          call cpseti('CLU - line type',3)
        enddo

      endif


!  we will call set
      call cpseti('SET',0)

!  don't label highs and lows
      call cpsetc('HLT - high/low label',' ')

!  initialize CONPACK
      call cprect(zdat,K,M,N,rwrk,lrwk,iwrk,liwk)
      call cppkcl(zdat,rwrk,iwrk)

!  go through the contours, set dash patterns, and shift if desired
      call cpgetr('CIU - contour interval used',ciu)
      call cpgeti('NCL - number of contour levels',ncl)

      do 100 i=1,ncl
        call cpseti('PAI - array index',i)
        call cpgeti('CLU - contour level use flag',iclu)
        call cpgetr('CLV - contour level',clev)
        if (shift) then
          clev = clev + .5*ciu
        end if
        call cpsetr('CLV - contour level',clev)
        if (clev .lt. 0.) then
          if (scalel .gt. 1.5) then
            call cpseti('CLD - contour line dash pattern',13107)
          else
            call cpseti('CLD - contour line dash pattern',21845)
          end if
        else
          call cpseti('CLD - contour line dash pattern',65535)
        end if
        call cpseti('CLU - contour level use flag',iclu)
  100 continue
      call cpsetr('CIU - contour interval used',ciu)
      call cpsetr('ILS - information label size',.012*scalel)

      call cpseti('LLP',1)
      call cpgeti('LLP - line label positioning flag',llpf)

! draw the contour lines, masking them if necessary.
      if (llpf.le.1) then
        call cpcldr (zdat,rwrk,iwrk)
      else
        call arinam (iama,lama)
        call cplbam (zdat,rwrk,iwrk,iama)
        call cpcldm (zdat,rwrk,iwrk,iama,cpdrpl)
      end if

!  draw labels
      call gsclip(0)
      call cplbdr(zdat,rwrk,iwrk)
      call gsclip(1)

      return
      end

! *******************************************************************

      subroutine cpsfill(zdat,K,M,N,ncontr,lcflag,tallflg,fixflag)

!  John Wilkin (wilkin@flood.ml.csiro.au)
!  Kate Hedstrom (kate@ahab.rutgers.edu)
!  Produces colour fill plots

!  abs(lcflag) determines which colour bar to use
!  If lcflag < 0 solid contour lines are overplotted on the colours

!  This version is presently set-up to request 'ncontr' contour levels
!  splitting the range into 'ncontr+1' equal bands

!  tallflg added for vertical label bars

!  PROBLEMS:
!  Since conpack can end up selecting more than ncontr contours (see the
!  documentation on conpack parameter CLS) problems can arise
!  when more levels are selected than entries defined in the colour table
!  In this event a warning is printed that cpsfill ran out of colours

!  To add more colours increase ncolor in cpsfill, defclr and colram,
!  and add more rgb colour triplets to the data statements in defclr.
!  Don't get caught by the fortran limit on 19 continuations!

!  **********************************************************************

      integer         K, M, N, ncontr, lcflag
      real            zdat(K,*)
      logical         tallflg, fixflag

!  dimensions of work arrays
      integer         lrwk, liwk, lawk, nlevb, lcra
      parameter (     lrwk=100000, liwk=100000, lawk=5000000, nlevb=30, &
     &                lcra=100000                                    )
      integer         iwrk(liwk), iama(lawk), iara(12), igra(nlevb)
      real            rwrk(lrwk), xcra(lcra), ycra(lcra)
      integer         iasf(13)

      integer         ncolor
      parameter (     ncolor=17      )
      real            cdepth(ncolor)
      real            ciu, clev, spval
      real            fl, fr, fb, ft, ul, ur, ub, ut
      integer         i, ncl, ncused, iclu, ll
      external        colram

!  labels and indicies for labelbar
      character*9     llbs(ncolor-1)
      integer         lfin(ncolor)

! Extra parameters if we decide to hardwire the contour levels
      real            cpmin, cpmax, clevel
      integer         nclevels
      common / cfix / cpmin, cpmax, nclevels

      data iasf  / 13*1 /
      data cdepth / 20, 30, 40, 50, 60, 100, 200, 400, 600, 800, &
     &              1000, 1500, 2000, 2500, 3000, 4000, 5000 /
!     data cdepth / 19, 50, 75, 100, 150, 200, 400, 600, 800,
!    &              1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500 /

      spval = 1.e+20

!  we will call set
      call cpseti('SET',0)

!  set all aspect source flags to "individual"
      call gsasf(iasf)

!  force solid fill
      call gsfais(1)

!  define colour indicies
      call defclr(lcflag)

!  define the list of indicies required by the labelbar routine
      do i=1,ncolor
        lfin(i)=i+1
      end do

!  don't label highs and lows or lines
      call cpsetc('HLT - high/low label',' ')
      call cpseti('LLP - line label position',0)

!  kludge for NCAR version 3.1.3
      call cpseti('PAI - parameter array identifier  ',-1)
      call cpseti('AIA - AREA identifier outside grid',-1)

!  turn on drawing of thick edge of plot (will be drawn if contour lines
!  requested -- lcflag < 0 )
      call cpseti('PAI - array index',-1)
      call cpseti('CLU - contour level use',1)
      call cpsetr('CLL - contour line width',2.)

!  arrange for selection of contour levels

      if (fixflag) then
        call cpseti('CLS - contour level selector',0)
        call cpseti('NCL - number of contour levels', ncontr)

        do i=1,ncontr
          call cpseti('PAI - array index',i)
          call cpsetr('CLV - contour level',cdepth(i))
          call cpseti('CLU - line type',3)
        enddo
      else
        if (nclevels .eq. 0) then

!  default scheme -- request conpack find ncontr "nice" contour levels
          if (ncontr .ge. ncolor) then
            print *,' cpsfill: ncontr > ncolors ... trouble ahead'
          endif
          call cpseti('CLS - contour level selector',ncontr)
          call cpsetr('CIS - contour interval selector',0.)

        else

!  override selection hints
!  cpmin, cpmax and nclevels are input to plothist

          call cpseti('CLS - contour level selector',0)
          call cpseti('NCL - number of contour levels',nclevels+1)
          do i=1,nclevels+1
            clevel = cpmin + float(i-1)/float(nclevels)*(cpmax-cpmin)
            call cpseti('PAI - array index',i)
            call cpsetr('CLV - contour level',clevel)
            call cpseti('CLU - line type',3)
          enddo

        endif
      endif

!  initialize CONPACK
      call cprect(zdat,K,M,N,rwrk,lrwk,iwrk,liwk)

!  draw contour plot
      call arinam(iama,lawk)
      call cpclam(zdat,rwrk,iwrk,iama)
      call arscam(iama,xcra,ycra,lcra,iara,igra,nlevb,colram)

!  find contour levels and enter these in labels for labelbar
      call cpgetr('CIU - contour interval used',ciu)
      call cpgeti('NCL - number of contour levels',ncl)
      if (ncl .gt. ncolor-1) then
        print *,' cpsfill: warning -- ran out of colours'
      endif
      ncused = min(ncolor-1,ncl)

      do 110 i=1,ncused
        call cpseti('PAI - array index',i)
        call cpgeti('CLU - contour level use flag',iclu)
        call cpgetr('CLV - contour level',clev)
        write(llbs(i),100)clev
  100   format(1pg9.2)
  110 continue

!  create labelbar

      call sfseti('TYPE of FILL',0)

      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      if (tallflg) then
        call lblbar(1,fr+0.01,fr+0.09,fb,ft,ncused+1,.15,1.,lfin,0, &
     &            llbs,ncused,1)
      else
        call lblbar(0,fl,fr,fb-0.05,fb-0.01,ncused+1,1.,.333,lfin,0, &
     &            llbs,ncused,1)
      end if

      call gsplci(1)

!  put black contour lines over the colored map
      if (lcflag .lt. 0) call cpcldr(zdat,rwrk,iwrk)

!  draw labels
!     call gsclip(0)
!     call cplbdr(zdat,rwrk,iwrk)
!     call gsclip(1)

!  reset min, max
      call cpsetr('CMN - contour min',1.)
      call cpsetr('CMX - contour max',0.)

      return
      end

!  ********************************************************************

      subroutine defclr(lcflag)
      integer    lcflag

!  defines the colour indicies

      integer         ncolor, krgb, i
      parameter (     ncolor=17      )
      real            rgb(3,ncolor,6)
      real            rgb1(3,ncolor),rgb2(3,ncolor),rgb3(3,ncolor), &
     &                rgb4(3,ncolor),rgb5(3,ncolor),rgb6(3,ncolor)
      equivalence     (rgb(1,1,1),rgb1(1,1)),(rgb(1,1,2),rgb2(1,1))
      equivalence     (rgb(1,1,3),rgb3(1,1)),(rgb(1,1,4),rgb4(1,1))
      equivalence     (rgb(1,1,5),rgb5(1,1)),(rgb(1,1,6),rgb6(1,1))

!  Kate's colors (red to blue)
      data rgb1  /    0.00, 0.25, 1.00, &
     &                0.00, 0.50, 1.00, &
     &                0.00, 0.75, 1.00, &
     &                0.00, 1.00, 1.00, &
     &                0.25, 1.00, 1.00, &
     &                0.50, 1.00, 1.00, &
     &                0.90, 0.90, 0.90, &
     &                1.00, 0.63, 0.63, &
     &                1.00, 0.50, 0.50, &
     &                1.00, 0.38, 0.38, &
     &                1.00, 0.25, 0.25, &
     &                1.00, 0.13, 0.13, &
     &                1.00, 0.00, 0.00, &
     &                1.00, 0.00, 0.00, &
     &                1.00, 0.00, 0.00, &
     &                1.00, 0.00, 0.00, &
     &                1.00, 0.00, 0.00 /
!  John's colours (grey to blue like GEBCO bathymetry charts)
!  second attempt -- not so dark on the Tek printer
      data rgb2    /  0.65, 0.60, 0.60, &
     &                0.65, 0.75, 0.75, &
     &                0.65, 0.90, 0.90, &
     &                0.60, 1.00, 1.00, &
     &                0.50, 1.00, 1.00, &
     &                0.40, 1.00, 1.00, &
     &                0.25, 1.00, 1.00, &
     &                0.10, 1.00, 1.00, &
     &                0.05, 0.95, 1.00, &
     &                0.00, 0.85, 1.00, &
     &                0.00, 0.75, 1.00, &
     &                0.00, 0.65, 1.00, &
     &                0.00, 0.50, 1.00, &
     &                0.00, 0.40, 0.85, &
     &                0.00, 0.30, 0.75, &
     &                0.00, 0.20, 0.60, &
     &                0.00, 0.00, 0.50 /
!  John's colours (modified for Dale's printer)
      data rgb3    /  0.65, 0.60, 0.60, &
     &                0.65, 0.75, 0.75, &
     &                0.65, 0.90, 0.90, &
     &                0.60, 1.00, 1.00, &
     &                0.50, 1.00, 1.00, &
     &                0.40, 1.00, 1.00, &
     &                0.20, 1.00, 1.00, &
     &                0.00, 1.00, 1.00, &
     &                0.05, 0.95, 1.00, &
     &                0.00, 0.80, 1.00, &
     &                0.00, 0.60, 1.00, &
     &                0.00, 0.40, 1.00, &
     &                0.00, 0.20, 1.00, &
     &                0.00, 0.00, 1.00, &
     &                0.00, 0.00, 0.80, &
     &                0.00, 0.00, 0.60, &
     &                0.00, 0.00, 0.40 /
!  Kate's colours (for Arctic psi)
      data rgb4    /  0.10, 1.00, 0.10, &
     &                0.32, 1.00, 0.30,       &! -8 to -6
     &                0.47, 1.00, 0.42, &
     &                0.62, 1.00, 0.54,       &! -4 to -2
     &                0.77, 1.00, 0.66, &
     &                1.00, 0.97, 0.70,       &! 0 to 2
     &                1.00, 0.84, 0.50, &
     &                1.00, 0.72, 0.30,       &! 4 to 6
     &                0.90, 0.60, 0.30, &
     &                0.86, 0.54, 0.40,       &! 8 to 10
     &                0.77, 0.47, 0.50, &
     &                0.69, 0.45, 0.64,       &! 12 to 14
     &                0.62, 0.41, 0.74, &
     &                0.58, 0.32, 0.80,       &! 16 to 18
     &                0.50, 0.28, 0.84, &
     &                0.44, 0.20, 0.88,       &! not used
     &                0.00, 0.00, 0.60 /
!  shades of grey - light to dark
      data rgb5    /  0.94, 0.94, 0.94, &
     &                0.88, 0.88, 0.88, &
     &                0.82, 0.82, 0.82, &
     &                0.76, 0.76, 0.76, &
     &                0.70, 0.70, 0.70, &
     &                0.64, 0.64, 0.64, &
     &                0.58, 0.58, 0.58, &
     &                0.52, 0.52, 0.52, &
     &                0.46, 0.46, 0.46, &
     &                0.40, 0.40, 0.40, &
     &                0.34, 0.34, 0.34, &
     &                0.28, 0.28, 0.28, &
     &                0.22, 0.22, 0.22, &
     &                0.16, 0.16, 0.16, &
     &                0.10, 0.10, 0.10, &
     &                0.05, 0.05, 0.05, &
     &                0.00, 0.00, 0.00 /
!  shades of grey - dark to light
      data rgb6    /  0.16, 0.16, 0.16, &
     &                0.22, 0.22, 0.22, &
     &                0.28, 0.28, 0.28, &
     &                0.34, 0.34, 0.34, &
     &                0.40, 0.40, 0.40, &
     &                0.46, 0.46, 0.46, &
     &                0.52, 0.52, 0.52, &
     &                0.58, 0.58, 0.58, &
     &                0.64, 0.64, 0.64, &
     &                0.70, 0.70, 0.70, &
     &                0.76, 0.76, 0.76, &
     &                0.82, 0.82, 0.82, &
     &                0.88, 0.88, 0.88, &
     &                0.94, 0.94, 0.94, &
     &                0.96, 0.96, 0.96, &
     &                0.98, 0.98, 0.98, &
     &                1.00, 1.00, 1.00 /

      krgb = abs(lcflag)
      if (krgb .gt. 6) call crash(' cpsfill: invalid lcflag',lcflag)

      do 100 i=1,ncolor
        call gscr(1,i+1,rgb(1,i,krgb),rgb(2,i,krgb),rgb(3,i,krgb))
  100 continue

      return
      end


! ********************************************************************

      subroutine colram(xcra,ycra,ncra,iaia,igia,naia)
      real            xcra(*), ycra(*)
      integer         iaia(*), igia(*)
      integer         ncolor, ncra, naia, ifll, i
      parameter (     ncolor=17      )

!  The arrays xcra and ycra, for indices 1 to ncra, contain the x
!  and y coordinates of points defining a polygon.  The area
!  identifiers in the array iaia, each with an associated group
!  identifier in the array igia, tell us whether the polygon
!  is to be color-filled or not.

      ifll = 1

!  if any of the area identifiers is negative, don't fill the
!  polygon.
      do 100 i=1,naia
        if (iaia(i) .lt. 0) ifll = 0
  100 continue

!  otherwise, fill the polygon in the color implied by its area
!  identifier relative to edge group 3 (the contour-line group)

      if (ifll .ne. 0) then
        ifll = 0
        do 110 i=1,naia
          if (igia(i) .eq. 3) ifll = iaia(i)
  110   continue
        if (ifll .gt. 0 .and. ifll .le. ncolor) then
          call gsfaci(ifll+1)
          call gfa(ncra-1, xcra, ycra)
        end if
      end if
!     print *,'in colram'
!     do i=1,naia
!       print *,'iaia ', i, iaia(i), igia(i)
!     enddo
!     print *,'ifll', ifll

      return
      end
