#include "griddefs.h"
      subroutine def_grid
!
!=======================================================================
!                                                                    ===
!  This routine creates a Gridpak NetCDF file, it defines its        ===
!  dimensions, attributes, and variables.                            ===
!                                                                    ===
!                                                                    ===
!  Calls:  crash, opencdf, plus NetCDF library                       ===
!                                                                    ===
!  WARNING:  The character argument to the NetCDF routines are       ===
!            (upper/lower) case sensitive.                           ===
!                                                                    ===
!=======================================================================
!
!-----------------------------------------------------------------------
!  Define global variables.
!-----------------------------------------------------------------------
!
      use netcdf

#include "grid.h"
#include "ncgrid.h"
!
!-----------------------------------------------------------------------
!  Define local variables.
!-----------------------------------------------------------------------
!
      integer bathdim, epdim, erdim, eudim, evdim, &
     &        onedim, twodim, xpdim, xrdim, xudim, xvdim
      integer h2dgrd(3), p2dgrd(2), t2dgrd(2), &
     &        u2dgrd(2), v2dgrd(2)
      BIGREAL c0, c1
      parameter  ( c0=0.d0, c1=1.d0 )
!
!=======================================================================
!  Begin executable code.
!=======================================================================
!
!=======================================================================
!  Create a new grid NetCDF file.
!=======================================================================
!
      rcode = nf90_create(TRIM(gridfile), nf90_clobber, ncgridid)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
      if ((rcode.ne.0).or.(ncgridid.eq.-1)) then
        write(stdout,10) TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      rcode = nf90_def_dim(ncgridid,'xi_psi',L,xpdim)
      rcode = nf90_def_dim(ncgridid,'xi_rho',Lp,xrdim)
      rcode = nf90_def_dim(ncgridid,'xi_u',L,xudim)
      rcode = nf90_def_dim(ncgridid,'xi_v',Lp,xvdim)
      rcode = nf90_def_dim(ncgridid,'eta_psi',M,epdim)
      rcode = nf90_def_dim(ncgridid,'eta_rho',Mp,erdim)
      rcode = nf90_def_dim(ncgridid,'eta_u',Mp,eudim)
      rcode = nf90_def_dim(ncgridid,'eta_v',M,evdim)
      rcode = nf90_def_dim(ncgridid,'one',1,onedim)
      rcode = nf90_def_dim(ncgridid,'two',2,twodim)
      rcode = nf90_def_dim(ncgridid,'bath',nf90_unlimited,bathdim)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Define dimension vectors for staggered 2D psi type
!  variables.
!
      p2dgrd(1)=xpdim
      p2dgrd(2)=epdim
!
!  Define dimension vectors for staggered 2D tracer type
!  variables.
!
      t2dgrd(1)=xrdim
      t2dgrd(2)=erdim
      h2dgrd(1)=xrdim
      h2dgrd(2)=erdim
      h2dgrd(3)=bathdim
!
!  Define dimension vectors for staggered 2D u-momemtum type
!  variables.
!
      u2dgrd(1)=xudim
      u2dgrd(2)=eudim
!
!  Define dimension vectors for staggered 2D v-momemtum type
!  variables.
!
      v2dgrd(1)=xvdim
      v2dgrd(2)=evdim
!
!-----------------------------------------------------------------------
!  Define type of floating-point variables: single or double precision.
!-----------------------------------------------------------------------
!
#if DBLEPREC || defined cray
      vartyp=nf90_double
#else
      vartyp=nf90_float
#endif /* DBLEPREC */
!
!-----------------------------------------------------------------------
!  Create history attribute.
!-----------------------------------------------------------------------
!
      call get_date (date_str)
      if (LEN_TRIM(date_str).gt.0) then
        history='Gridpak, Version '//version//', '//TRIM(date_str)
      else
        history='Gridpak, Version '//version
      endif
      rcode = nf90_put_att(ncgridid,nf90_global,'type',&
     &              'Gridpak file')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Put global attributes to NetCDF file.
!
      rcode = nf90_put_att(ncgridid,nf90_global,'gridid', gridid)
      rcode = nf90_put_att(ncgridid,nf90_global,'history', &
     &               history)
      rcode = nf90_put_att(ncgridid,nf90_global,'CPP-options', &
     &               CPPoptions)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Domain Length.
!
        rcode = nf90_def_var(ncgridid,'xl',vartyp,varid=varid)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
        rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'domain length in the XI-direction')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
        rcode = nf90_put_att(ncgridid,varid,'units', 'meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
        rcode = nf90_def_var(ncgridid,'el',vartyp,varid=varid)
        rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'domain length in the ETA-direction')
        rcode = nf90_put_att(ncgridid,varid,'units', 'meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Projection parameters.
!
      rcode = nf90_def_var(ncgridid,'JPRJ',nf90_char,(/twodim/),varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Map projection type')
      rcode = nf90_put_att(ncgridid,varid,'option_ME', 'Mercator')
      rcode = nf90_put_att(ncgridid,varid,'option_ST', 'Stereographic')
      rcode = nf90_put_att(ncgridid,varid,'option_LC', &
     &              'Lambert conformal conic')
      rcode=nf90_def_var(ncgridid,'PLAT',nf90_double,(/twodim/),varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Reference latitude(s) for map projection')
      rcode = nf90_put_att(ncgridid,varid,'units', 'degree_north')
      rcode=nf90_def_var(ncgridid,'PLONG',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Reference longitude for map projection')
      rcode = nf90_put_att(ncgridid,varid,'units', 'degree_east')
      rcode=nf90_def_var(ncgridid,'ROTA',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Rotation angle for map projection')
      rcode = nf90_put_att(ncgridid,varid,'units', 'degree')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif

      rcode=nf90_def_var(ncgridid,'JLTS',nf90_char,(/twodim/),varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'How limits of map are chosen')
      rcode = nf90_put_att(ncgridid,varid,'option_CO', &
     &              'P1, .. P4 define two opposite corners of domain')
      rcode = nf90_put_att(ncgridid,varid,'option_MA', &
     &              'Maximum (whole world)')
      rcode = nf90_put_att(ncgridid,varid,'option_AN', &
     &              'Angles - P1..P4 define angles to edge of domain')
      rcode = nf90_put_att(ncgridid,varid,'option_LI', &
     &              'Limits - P1..P4 define limits in u,v space')
      rcode=nf90_def_var(ncgridid,'P1',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Map limit parameter number 1')
      rcode=nf90_def_var(ncgridid,'P2',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Map limit parameter number 2')
      rcode=nf90_def_var(ncgridid,'P3',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Map limit parameter number 3')
      rcode=nf90_def_var(ncgridid,'P4',nf90_double,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Map limit parameter number 4')
!     rcode=nf90_def_var(ncgridid,'JGRD',nf90_double,varid=varid)
!     rcode = nf90_put_att(ncgridid,varid,'long_name',
!    &              'Density of lat,lon lines to draw on plots')
!     rcode = nf90_put_att(ncgridid,varid,'units', 'degree')
      rcode=nf90_def_var(ncgridid,'XOFF',nf90_float,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Offset in x direction')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode = nf90_def_var(ncgridid,'YOFF',nf90_float,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Offset in y direction')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  depthmin and depthmax.
!
      rcode = nf90_def_var(ncgridid,'depthmin',vartyp,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Shallow bathymetry clipping depth')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode = nf90_def_var(ncgridid,'depthmax',vartyp,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Deep bathymetry clipping depth')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Spherical geometry logical flag and beta-plane parameters.
!
      rcode=nf90_def_var(ncgridid,'spherical',nf90_int,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Grid type logical switch')
      rcode = nf90_put_att(ncgridid,varid,'flag_values', '0, 1')
      rcode = nf90_put_att(ncgridid,varid,'flag_meanings',  &
     &              'Cartesian, spherical')
      rcode=nf90_def_var(ncgridid,'f0',vartyp,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Coriolis parameter central value on a beta-plane')
      rcode = nf90_put_att(ncgridid,varid,'_FillValue',c0)
      rcode=nf90_def_var(ncgridid,'dfdy',vartyp,varid=varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Coriolis parameter gradient on a beta-plane')
      rcode = nf90_put_att(ncgridid,varid,'_FillValue',c0)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Bathymetry.
!
      rcode=nf90_def_var(ncgridid,'hraw',vartyp,h2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Working bathymetry at RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'h',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Final bathymetry at RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Coriolis Parameter.
!
      rcode=nf90_def_var(ncgridid,'f',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Coriolis parameter at RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units', 'second-1')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Curvilinear coordinates metrics.
!
      rcode=nf90_def_var(ncgridid,'pm',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'curvilinear coordinate metric in XI')
      rcode = nf90_put_att(ncgridid,varid,'units', 'meter-1')
      rcode=nf90_def_var(ncgridid,'pn',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'curvilinear coordinate metric in ETA')
      rcode = nf90_put_att(ncgridid,varid,'units', 'meter-1')
      rcode=nf90_def_var(ncgridid,'dndx',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'XI-derivative of inverse metric factor pn')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'dmde',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'ETA-derivative of inverse metric factor pm')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  x,y location on the grids.
!
      rcode=nf90_def_var(ncgridid,'x_rho',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'X-location of RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'y_rho',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Y-location of RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'x_psi',vartyp,p2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'X-location of PSI-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'y_psi',vartyp,p2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Y-location of PSI-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'x_u',vartyp,u2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'X-location of U-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'y_u',vartyp,u2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Y-location of U-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'x_v',vartyp,v2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'X-location of V-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      rcode=nf90_def_var(ncgridid,'y_v',vartyp,v2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'Y-location of V-points')
      rcode = nf90_put_att(ncgridid,varid,'units','meter')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Latitude, longitude on the grids.
!
      rcode=nf90_def_var(ncgridid,'lat_rho',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'latitude of RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_north')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'latitude')
      rcode=nf90_def_var(ncgridid,'lon_rho',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'longitude of RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_east')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'longitude')
      rcode=nf90_def_var(ncgridid,'lat_psi',vartyp,p2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'latitude of PSI-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_north')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'latitude')
      rcode=nf90_def_var(ncgridid,'lon_psi',vartyp,p2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'longitude of PSI-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_east')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'longitude')
      rcode=nf90_def_var(ncgridid,'lat_u',vartyp,u2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'latitude of U-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_north')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'latitude')
      rcode=nf90_def_var(ncgridid,'lon_u',vartyp,u2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'longitude of U-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_east')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'longitude')
      rcode=nf90_def_var(ncgridid,'lat_v',vartyp,v2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'latitude of V-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_north')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'latitude')
      rcode=nf90_def_var(ncgridid,'lon_v',vartyp,v2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'longitude of V-points')
      rcode = nf90_put_att(ncgridid,varid,'units', &
     &              'degree_east')
      rcode = nf90_put_att(ncgridid,varid,'standard_name', &
     &              'longitude')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  0/1 masks.
!
      rcode=nf90_def_var(ncgridid,'mask_rho',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'mask on RHO-points')
      rcode = nf90_put_att(ncgridid,varid,'flag_values', '0., 1.')
      rcode = nf90_put_att(ncgridid,varid,'flag_meanings', 'land, water')
      rcode=nf90_def_var(ncgridid,'mask_u',vartyp,u2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'mask on U-points')
      rcode = nf90_put_att(ncgridid,varid,'flag_values', '0., 1.')
      rcode = nf90_put_att(ncgridid,varid,'flag_meanings', 'land, water')
      rcode=nf90_def_var(ncgridid,'mask_v',vartyp,v2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'mask on V-points')
      rcode = nf90_put_att(ncgridid,varid,'flag_values', '0., 1.')
      rcode = nf90_put_att(ncgridid,varid,'flag_meanings', 'land, water')
      rcode=nf90_def_var(ncgridid,'mask_psi',vartyp,p2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'mask on PSI-points')
      rcode = nf90_put_att(ncgridid,varid,'flag_values', '0., 1.')
      rcode = nf90_put_att(ncgridid,varid,'flag_meanings', 'land, water')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!  Angle.
!
      rcode=nf90_def_var(ncgridid,'angle',vartyp,t2dgrd,varid)
      rcode = nf90_put_att(ncgridid,varid,'long_name', &
     &              'angle between XI-axis and EAST')
      rcode = nf90_put_att(ncgridid,varid,'units','radian')
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
      rcode = nf90_enddef(ncgridid)
      if (rcode.ne.0) then
        print *, trim(nf90_strerror(rcode))
      endif
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Curvilinear transformation metrics.
!
      spherical = 0
      rcode=nf90_inq_varid(ncgridid,'spherical',varid)
      rcode = nf90_put_var(ncgridid,varid,spherical)
      if (rcode.ne.0) then
        write(stdout,20) 'spherical', TRIM(gridfile)
        print *, trim(nf90_strerror(rcode))
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'pm',varid)
      rcode = nf90_put_var(ncgridid,varid,pm)
      if (rcode.ne.0) then
        write(stdout,20) 'pm', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'pn',varid)
      rcode = nf90_put_var(ncgridid,varid,pn)
      if (rcode.ne.0) then
        write(stdout,20) 'pn', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'dndx',varid)
      rcode = nf90_put_var(ncgridid,varid,dndx)
      if (rcode.ne.0) then
        write(stdout,20) 'dndx', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'dmde',varid)
      rcode = nf90_put_var(ncgridid,varid,dmde)
      if (rcode.ne.0) then
        write(stdout,20) 'dmde', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
!
!  Domain Length.
!
      rcode=nf90_inq_varid(ncgridid,'xl',varid)
      rcode = nf90_put_var(ncgridid,varid,xl)
      if (rcode.ne.0) then
        write(stdout,20) 'xl', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'el',varid)
      rcode = nf90_put_var(ncgridid,varid,el)
      if (rcode.ne.0) then
        write(stdout,20) 'el', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
!
!  x,y grids.
!
      rcode=nf90_inq_varid(ncgridid,'x_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,xr)
      if (rcode.ne.0) then
        write(stdout,20) 'x_rho', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_rho',varid)
      rcode = nf90_put_var(ncgridid,varid,yr)
      if (rcode.ne.0) then
        write(stdout,20) 'y_rho', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,xp)
      if (rcode.ne.0) then
        write(stdout,20) 'x_psi', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_psi',varid)
      rcode = nf90_put_var(ncgridid,varid,yp)
      if (rcode.ne.0) then
        write(stdout,20) 'y_psi', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_u',varid)
      rcode = nf90_put_var(ncgridid,varid,xu)
      if (rcode.ne.0) then
        write(stdout,20) 'x_u', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_u',varid)
      rcode = nf90_put_var(ncgridid,varid,yu)
      if (rcode.ne.0) then
        write(stdout,20) 'y_u', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'x_v',varid)
      rcode = nf90_put_var(ncgridid,varid,xv)
      if (rcode.ne.0) then
        write(stdout,20) 'x_v', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
      rcode=nf90_inq_varid(ncgridid,'y_v',varid)
      rcode = nf90_put_var(ncgridid,varid,yv)
      if (rcode.ne.0) then
        write(stdout,20) 'y_v', TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
!
!  Close the file.
!
      rcode = nf90_close(ncgridid)
      if (rcode.ne.0) then
        write(stdout,30) TRIM(gridfile)
        call crash ('DEF_GRID',1)
      endif
!
!=======================================================================
!
  10  format(/' DEF_GRID - unable to create grid NetCDF file: ',a)
  20  format(/' DEF_GRID - error while writing variable: ',a,/,11x, &
     &        'into grid NetCDF file: ',a)
  30  format(/' DEF_GRID - unable to close grid NetCDF file: ',a)
      return
      end
