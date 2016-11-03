!
!=======================================================================
!  Include file "ncscrum.h"
!=======================================================================
!
!  maxvar      Maximum number of variables for input NetCDF files.
!  ncgridid    NetCDF ID for grid file.
!  nvars       Number of variables defined in current input NetCDF file.
!  nvdims      Number of dimensions for each variables in current input
!              NetCDF file.
!  rcode       Error code returned by NetCDF library (0 for no errors).
!  bathsize    Size of unlimited time record dimension in current input
!              NetCDF file.
!  varid       Generic ID for arbitrary variables in NetCDF files.
!  varnam      Names of all variables in current input NetCDF file.
!  vdims       Dimension IDs for each of the variables in current input
!              NetCDF file.
!
!=======================================================================
!
      integer maxvar
      parameter (maxvar=100)
!
      integer           bathindex, nvdims(maxvar), &
     &                  vdims(5,maxvar), nvars, ncgridid, rcode, &
     &                  varid, vartyp, bathsize
      common /incgrid/  bathindex, nvdims, vdims, nvars, &
     &                  ncgridid, rcode, varid, vartyp, &
     &                  bathsize
!
      character*5       version
      integer           patchlevel, stdout
      parameter   (     version='5.3  '  )
      parameter   (     patchlevel=0 )
      parameter   (     stdout=6 )

      character*44      date_str
      character*120     history
      character*15      varnam(maxvar)
      common /cncgrid/  date_str, history, varnam
      character*1024    CPPoptions
      character*80      gridfile, grid1_file, gridid, type
      common /strings/  CPPoptions, gridfile, grid1_file, gridid, type
