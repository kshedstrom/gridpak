! define one of these as 1 for reading the respective bathymetry.
! Otherwise, it defaults to using the BATHY_FILE environment variable,
! here with the dimensions of the Smith and Sandwell subset for USwest.
#undef ETOPO5
#undef ETOPO2
#define ETOPO1 1
#undef GEBCO
#undef IBCAO
#undef ARDEM

! for 64-bit output
#define DBLEPREC      1

! to draw coastlines on some plots
#define DRAW_COASTS   1

! to keep ellipsoidal terms in Earth's shape
#define ELLIPSOID     1

! for averaging bathymetry in gridbox (for EW/NS grids only)
#undef IMG_AVG

#undef KEEP_SHALLOW

! for NCAR graphics (3.2 or better) */
#define PLOTS	     1
! for X windows rather than metafile */
#undef X_WIN

#undef SYS_POTS       /* unimplimented system calls */
#undef XPOTS1	      /* read ipot1 file */
#undef XPOTS2	      /* read ipot2 file */
#undef XPOTS3	      /* read ipot3 file */
#undef XPOTS4	      /* read ipot4 file */

#ifdef cray
#undef DCOMPLEX
#define DBLEPREC      1	/* for 64-bit output */
#define BIGREAL real
#define SMALLREAL real
#define BIGCOMPLEX complex
#define FLoaT float
#else
#if DBLEPREC
#define DCOMPLEX      1    /* for compilers which support complex*16 */
#define SMALLREAL real
#define BIGREAL real*8
#define BIGCOMPLEX complex*16
#define FLoaT dfloat
#else
#undef DCOMPLEX
#define BIGREAL real
#define SMALLREAL real
#define BIGCOMPLEX complex
#define FLoaT float
#endif  /* DBLEPREC */
#endif  /* cray */
