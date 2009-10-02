#include "griddefs.h"
      subroutine opencdf (ncname)
!
!=======================================================================
!                                                                    ===
!  This routine opens an existing NetCDF file and inquires about it  ===
!  contents, and checks for consistency with model dimensions.       ===
!                                                                    ===
!  On Input:                                                         ===
!                                                                    ===
!     ncname   Input NetCDF file name.                               ===
!                                                                    ===
!  Calls:   (NetCDF library)                                         ===
!           crash                                                    ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global variables.
!-----------------------------------------------------------------------
!
      use netcdf

#include "bathy.h"
#include "ncgrid.h"
!
!-----------------------------------------------------------------------
!  Define local variables.
!-----------------------------------------------------------------------
!
      integer attype, dimid, dimsiz, i, ncid, ndims, &
     &        ngatts, nvatts, recdim
      character*20  dimnam
      character*(*) ncname
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!  Open input NetCDF file.
!
      rcode=nf90_open(TRIM(ncname),nf90_nowrite,ncid)
      if ((rcode.ne.0).or.(ncid.eq.-1)) then
        write(stdout,10) TRIM(ncname)
        call crash ('OPENCDF',1)
      endif
!
!  Inquire and get global "type" attribute.
!
      rcode = nf90_inquire_attribute(ncid,nf90_global,'type',attype)
      if (rcode.eq.0) then
        rcode = nf90_get_att(ncid,nf90_global,'type',type)
        if(rcode.ne.0) then
          write(stdout,20) 'type (global)', TRIM(ncname)
          call crash ('OPENCDF',1)
        endif
      else
        write(stdout,30) 'type (global)',TRIM(ncname)
      endif
!
!  Inquire about the contents of input NetCDF file:  Inquire about the
!  dimensions and variables.
!
      rcode = nf90_inquire(ncid,ndims,nvars,ngatts,recdim)
      if (nvars.gt.maxvar) then
        write(stdout,40) maxvar, nvars
        call crash ('OPENCDF',1)
      endif
      if (rcode.eq.0) then
!
!  Inquire about dimensions.  Check dimensions for consistency.
!
        do i=1,ndims
          dimid=i
          rcode = nf90_inquire_dimension(ncid,dimid,dimnam,dimsiz)
          if (rcode.ne.0) then
            write(stdout,50) dimid, TRIM(ncname)
            call crash ('OPENCDF',1)
          endif
          if ((TRIM(dimnam) .eq. 'xi_rho') .or. &
     &        (TRIM(dimnam) .eq. 'xi_v')) then
            if (dimsiz.ne.Lp) then
              write(stdout,60) TRIM(dimnam), dimsiz, Lp
              call crash ('OPENCDF',1)
            endif
          elseif ((TRIM(dimnam) .eq. 'xi_u') .or. &
     &            (TRIM(dimnam) .eq. 'xi_psi')) then
            if (dimsiz.ne.L) then
              write(stdout,60) TRIM(dimnam), dimsiz, L
              call crash ('OPENCDF',1)
            endif
          elseif ((TRIM(dimnam) .eq. 'eta_rho') .or. &
     &            (TRIM(dimnam) .eq. 'eta_u')) then
            if (dimsiz.ne.Mp) then
              write(stdout,60) TRIM(dimnam), dimsiz, Mp
              call crash ('OPENCDF',1)
            endif
          elseif ((TRIM(dimnam) .eq. 'eta_v') .or. &
     &            (TRIM(dimnam) .eq. 'eta_psi')) then
            if (dimsiz.ne.M) then
              write(stdout,60) TRIM(dimnam), dimsiz, M
              call crash ('OPENCDF',1)
            endif
          endif
        enddo
!
!  Inquire about variables.
!
        do i=1,nvars
          varid=i
          rcode = nf90_inquire_variable(ncid,varid,varnam(i),vartyp, &
     &                 nvdims(i),vdims(:,i),nvatts)
          if (rcode.ne.0) then
            write(stdout,70) varid, TRIM(ncname)
            call crash ('OPENCDF',1)
          endif
        enddo
      else
        write(stdout,80) TRIM(ncname)
        call crash ('OPENCDF',1)
      endif
!
!-----------------------------------------------------------------------
!  Inquire size of unlimited time record dimension.
!-----------------------------------------------------------------------
!
      bathsize=0
      if (recdim.ne.-1) then
        rcode = nf90_inquire_dimension(ncid,recdim,len=bathsize)
        if (rcode.ne.0) then
          write(stdout,90) 'bath', TRIM(ncname)
          call crash ('OPENCDF',1)
        endif
      endif
!
!  Close input NetCDF file.
!
      rcode = nf90_close(ncid)
!
  10  format(/' OPENCDF - unable to open input NetCDF file: ',a)
  20  format(/' OPENCDF - error while reading attribute: ',a,2x, &
     &        ' in input NetCDF file: ',a)
  30  format(/' OPENCDF - cannot find attribute: ',a,2x, &
     &        ' in input NetCDF file: ',a)
  40  format(/' OPENCDF - too small dimension parameter, maxvar = ',2i5, &
     &       /,11x,'change file  ncscrum.h  and recompile.')
  50  format(/' OPENCDF - error while reading dimension ID: ',i3,2x, &
     &        ' in input NetCDF file: ',a)
  60  format(/' OPENCDF - inconsistent size of dimension: ',a,2x, &
     &       2i5)
  70  format(/' OPENCDF - error while inquiring information for ', &
     &        ' variable ID: ',i3,2x,' in input NetCDF file: ',a)
  80  format(/' OPENCDF - unable to inquire about contents of', &
     &        ' input NetCDF file: ',a)
  90  format(/' OPENCDF - error inquiring dimension: ',a,2x, &
     &        ' in input NetCDF file: ',a)
      return
      end
