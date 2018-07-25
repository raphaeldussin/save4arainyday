PROGRAM detrend_precip

!------------------------------------------------------
! A program which reads precip input, applies weights
! and corrections from Gilles Garric and outputs it
! as a netcdf file
!
! Raphael Dussin, fev 2011
!
!------------------------------------------------------

USE NETCDF

   IMPLICIT NONE

   INTEGER                                   :: narg, iargc                     ! command line arguments
   INTEGER                                   :: ncid_precip, ncid_trend, ncid_out  ! netcdf id's
   INTEGER                                   :: ncid_mask
   INTEGER                                   :: irec                            ! records loop index
   INTEGER                                   :: lonin_varid, latin_varid
   INTEGER                                   :: lonout_varid, latout_varid, tout_varid
   INTEGER                                   :: prein_varid, trend_varid, preout_varid
   INTEGER                                   :: mask_varid

   INTEGER                                   :: npx, npy, npt, npxw, npyw, nptw ! dimensions in precip and trend
   INTEGER                                   :: xid, yid, tid, xwid, ywid, twid ! netcdf dimensions id's
   INTEGER                                   :: xoutid, youtid, toutid

   INTEGER                                   :: start(3), count(3)

   INTEGER                                   :: jj, ji

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zprein, ztrend, zpreout
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: zlon, zlat, ztim
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zmask
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ztmp

   REAL(KIND=8)                              :: thisyear, startyear, tnow

   CHARACTER(LEN=256)                        :: cfilein, ctrendin
   CHARACTER(LEN=256), PARAMETER             :: cfileout="precip_detrended.nc"
   CHARACTER(LEN=64)                         :: cvar_rad="precip", cvar_wgt="trend", cvar_out="precip"
   CHARACTER(LEN=256)                        :: cmask="mask.nc", cthisyear, cstartyear

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()
   IF ( narg /= 4 ) THEN
      PRINT *,' Usage : detrend_precip precip_file trend_file current_year start_year '
      PRINT *,' mask.nc must be in your directory ' 
      STOP
   ENDIF

   CALL getarg (1, cfilein)
   CALL getarg (2, ctrendin )
   CALL getarg (3, cthisyear )
   CALL getarg (4, cstartyear )
   READ(cthisyear,*) thisyear
   READ(cstartyear,*) startyear

   PRINT *,'------------------------------------------------------------------'
   PRINT *,' Detrending file ', trim(cfilein), ' with coef ', trim(ctrendin)
   PRINT *,' This year is ', trim(cthisyear)
   PRINT *,' Detrending from ', trim(cstartyear)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Verify that input files have the same size

   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )

   CALL check( NF90_INQ_DIMID(ncid_precip, 'lon', xid) )
   CALL check( NF90_INQ_DIMID(ncid_precip, 'lat', yid) )
   CALL check( NF90_INQ_DIMID(ncid_precip, 'time', tid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, xid,len=npx) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, yid,len=npy) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, tid,len=npt) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in radsw file : npx = ', npx
   PRINT *, 'in radsw file : npy = ', npy
   PRINT *, 'in radsw file : npt = ', npt

   CALL check( NF90_CLOSE(ncid_precip) )

   !!--------------------------------------------------------------
   CALL check( NF90_OPEN(trim(ctrendin), nf90_nowrite, ncid_trend) )

   CALL check( NF90_INQ_DIMID(ncid_trend, 'lon', xwid) )
   CALL check( NF90_INQ_DIMID(ncid_trend, 'lat', ywid) )
   CALL check( NF90_INQ_DIMID(ncid_trend, 'time', twid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_trend, xwid,len=npxw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_trend, ywid,len=npyw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_trend, twid,len=nptw) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in weight file : npx = ', npxw
   PRINT *, 'in weight file : npy = ', npyw
   PRINT *, 'in weight file : npt = ', nptw
   PRINT *, '-------------------------------------------------------------'

   CALL check( NF90_CLOSE(ncid_trend) )

   IF ( npx /= npxw .OR. npy /= npyw ) THEN
      PRINT *, 'Problem with dimensions in input files'
      STOP
   ENDIF

   !--------------------------------------------------------------------
   ALLOCATE( zprein(npx,npy) , ztrend(npx,npy), zpreout(npx,npy) )
   ALLOCATE( zlon(npx) , zlat(npy) )
   ALLOCATE( zmask(npx,npy) )
   ALLOCATE( ztmp(npx,npy), ztim(npt) )

   count = (/ npx, npy, 1 /)
   start = (/ 1, 1, 1 /)

   !--------------------------------------------------------------------
   !! Get lon and lat
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )
   CALL check( NF90_INQ_VARID(ncid_precip, 'lon', lonin_varid) )
   CALL check( NF90_INQ_VARID(ncid_precip, 'lat', latin_varid) )
   CALL check( NF90_GET_VAR(ncid_precip, lonin_varid, zlon ))
   CALL check( NF90_GET_VAR(ncid_precip, latin_varid, zlat ))
   CALL check( NF90_CLOSE(ncid_precip) )
   
!   PRINT *, 'Max of lon is ', MAXVAL(zlon)
!   PRINT *, 'Max of lat is ', MAXVAL(zlat)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Get the mask and gphi from mask.nc and mesh_hgr.nc

   CALL check( NF90_OPEN(trim(cmask), nf90_nowrite, ncid_mask) )
   CALL check( NF90_INQ_VARID(ncid_mask, 'lsm', mask_varid) )
   CALL check( NF90_GET_VAR(ncid_mask, mask_varid, zmask ) )
   CALL check( NF90_CLOSE(ncid_mask) )


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  create the netcdf output file

   CALL check( nf90_create(trim(cfileout), nf90_clobber, ncid_out) )

   CALL check( nf90_def_dim(ncid_out, 'lon0', npx, xoutid) )
   CALL check( nf90_def_dim(ncid_out, 'lat0', npy, youtid) )
   CALL check( nf90_def_dim(ncid_out, 'time', NF90_UNLIMITED, toutid) )

   CALL check( nf90_def_var(ncid_out, 'lon0', NF90_REAL, (/ xoutid /), lonout_varid) )
   CALL check( nf90_def_var(ncid_out, 'lat0', NF90_REAL, (/ youtid /), latout_varid) )
   CALL check( nf90_def_var(ncid_out, 'time', NF90_REAL, (/ toutid /), tout_varid) )
   CALL check( nf90_def_var(ncid_out, trim(cvar_out), NF90_REAL, (/ xoutid, youtid, toutid /), preout_varid) )

   ! Assign units attributes to coordinate variables.
   CALL check( nf90_put_att(ncid_out, latout_varid, "units", "degrees_north") )
   CALL check( nf90_put_att(ncid_out, lonout_varid, "units", "degrees_east")  )
   CALL check( nf90_put_att(ncid_out, preout_varid, "units", "kg.m-2.s-1") )
   CALL check( nf90_put_att(ncid_out, tout_varid,   "units", "unknown") )


   ! End define mode.
   CALL check( nf90_enddef(ncid_out) )

   CALL check( nf90_put_var(ncid_out, latout_varid, zlat) )
   CALL check( nf90_put_var(ncid_out, lonout_varid, zlon) )

   DO irec=1,npt
      ztim(irec) = irec
   ENDDO

   CALL check( nf90_put_var(ncid_out, tout_varid, ztim) )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read the netcdf files

   ! radsw
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )
   CALL check( NF90_INQ_VARID(ncid_precip, trim(cvar_rad), prein_varid) )
   ! weight
   CALL check( NF90_OPEN(trim(ctrendin), nf90_nowrite, ncid_trend) )
   CALL check( NF90_INQ_VARID(ncid_trend, trim(cvar_wgt), trend_varid) )

   ! Read the trend once
   CALL check( NF90_GET_VAR(ncid_trend, trend_varid, ztrend, start = start, &
               &                 count = count) )

   !ztrend = ztrend / 86400 ! set in kg/m2/s

   DO irec=1,npt

      tnow = thisyear + (( (1. * irec) -1. ) / npt )

      PRINT *, 'Working on timestep : ', irec , ' ; time = ', tnow

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_precip, prein_varid, zprein, start = start, &
                  &                 count = count) )


      ztmp = zprein - ( ztrend * ( tnow - startyear ) )

      !! Avoid negative precip
      WHERE( ztmp < 0. ) ztmp = 0.

      !! this is my final array
      zpreout(:,:) = ztmp(:,:) * zmask(:,:)

      CALL check( NF90_PUT_VAR(ncid_out, preout_varid, REAL(zpreout), start = start, &
                  &            count = count) )


   ENDDO

   CALL check( NF90_CLOSE(ncid_precip) )
   CALL check( NF90_CLOSE(ncid_trend) )
   CALL check( nf90_close(ncid_out) )


CONTAINS

      SUBROUTINE check(status)
        INTEGER, INTENT(in) :: status
    
        IF(status /= nf90_noerr) THEN
          PRINT *, trim(nf90_strerror(status))
          STOP 2
        ENDIF
      END SUBROUTINE check

! shapiro modified to work in my framework

      SUBROUTINE Shapiro_1D(ptabin, kiter, cd_overlap, ptabout) 


      !!=====================================================================
      !!
      !! Description: This function applies a 1D Shapiro filter
      !!              (3 points filter) horizontally to a 2D field
      !!              in regular grid
      !! Arguments :
      !!            ptabin      : Input variable to filter
      !!            zla_mask    : Input mask variable
      !!            kiter       : Number of Shapiro filter iterations
      !!            cd_overlap  : Logical argument for periodical condition
      !!                          (global ocean case)
      !!            ptabout     : Output filtered variable
      !!
      !! History : 08/2009  S. CAILLEAU : from 1st version of N. FERRY
      !!           09/2009  C. REGNIER  : Corrections
      !!
      !!=====================================================================

      INTEGER,                  INTENT(IN)  :: kiter
      REAL(8), DIMENSION(:,:),  INTENT(IN)  :: ptabin   
      LOGICAL,                  INTENT(IN)  :: cd_overlap  
      REAL(8), DIMENSION(:,:),  INTENT(OUT) :: ptabout 

      ! * Local variable
      REAL(8), DIMENSION(:,:), ALLOCATABLE  :: zvarout
      REAL(8), PARAMETER                    :: rp_aniso_diff_XY=2.25 !  anisotrope case
      REAL(8)                               :: zalphax,zalphay, znum, zden,test
      INTEGER                               :: ji, jj, jn

! rp_aniso_diff_XY=2.25 : valeur trouvée empiriquement pour 140 itération pour
! le filtre de Shapiro et
! pour un rapport d'anisotopie de 1.5 : on filtre de plus rapidement en x qu'en
! y.
!------------------------------------------------------------------------------

    ALLOCATE( zvarout(npx,npy) )

!   Global ocean case : only this case is treated anyway
    IF ( cd_overlap ) THEN 

       ptabout(:,:) = ptabin(:,:)
       zvarout(:,:) = ptabin(:,:)       

       zalphax=1./2.
       zalphay=1./2.

!  Dx/Dy=rp_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rp_aniso_diff_XY correspondant est:

       IF ( rp_aniso_diff_XY >=  1. ) zalphay=zalphay/rp_aniso_diff_XY
       IF ( rp_aniso_diff_XY <   1. ) zalphax=zalphax*rp_aniso_diff_XY

       DO jn = 1,kiter
           DO jj = 2,npy-1
              DO ji = 2,npx-1
                 ! We crop on the coast
                  znum = zvarout(ji,jj)   &
                         + 0.25*zalphax*(zvarout(ji-1,jj  )-zvarout(ji,jj))*zmask(ji-1,jj  )  &
                         + 0.25*zalphax*(zvarout(ji+1,jj  )-zvarout(ji,jj))*zmask(ji+1,jj  )  &
                         + 0.25*zalphay*(zvarout(ji  ,jj-1)-zvarout(ji,jj))*zmask(ji  ,jj-1)  &
                         + 0.25*zalphay*(zvarout(ji  ,jj+1)-zvarout(ji,jj))*zmask(ji  ,jj+1)
                  ptabout(ji,jj)=znum*zmask(ji,jj)+ptabin(ji,jj)*(1.-zmask(ji,jj))
               ENDDO  ! end loop ji
           ENDDO  ! end loop jj

           ! we deal with east/west overlap of the ORCA grid
           ptabout(1,:)   = ptabout(npx-1,:)
           ptabout(npx,:) = ptabout(2,:)
           ! we leave the jpj=1 and jpj=npjglo points unfiltered
           ! jpj=1 is land anyway, jpj=npjglo could be fixed

           zvarout(:,:) = ptabout(:,:)

         ENDDO  ! end loop jn
      ENDIF


END SUBROUTINE Shapiro_1D     

END PROGRAM
