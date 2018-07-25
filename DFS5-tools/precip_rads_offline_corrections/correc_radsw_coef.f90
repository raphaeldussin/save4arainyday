PROGRAM correc_radsw_coef

!------------------------------------------------------
! A program which reads radiative input, applies coef
! and outputs it as a netcdf file
! ( DFS 5.1 stuff )
!
! Raphael Dussin, may 2012
!
!------------------------------------------------------

USE NETCDF

   IMPLICIT NONE

   INTEGER                                   :: narg, iargc                     ! command line arguments
   INTEGER                                   :: ncid_radsw, ncid_wgt, ncid_out  ! netcdf id's
   INTEGER                                   :: ncid_mesh, ncid_mask
   INTEGER                                   :: irec                            ! records loop index
   INTEGER                                   :: lonin_varid, latin_varid
   INTEGER                                   :: lonout_varid, latout_varid
   INTEGER                                   :: radin_varid, radout_varid
   INTEGER                                   :: wgt_varid, mask_varid, phi_varid

   INTEGER                                   :: npx, npy, npt, npxw, npyw, nptw ! dimensions in rad and wgt
   INTEGER                                   :: xid, yid, tid, xwid, ywid, twid ! netcdf dimensions id's
   INTEGER                                   :: xoutid, youtid, toutid
   INTEGER                                   :: tout_varid

   INTEGER                                   :: start(3), count(3)

   INTEGER                                   :: jj, ji

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zradin, zwgt, zradout
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zmask
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ztmp
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: zlon0, zlat0, ztim

   CHARACTER(LEN=256)                        :: cfilein, cwgtin
   CHARACTER(LEN=256), PARAMETER             :: cfileout="radsw_corrected.nc"
   CHARACTER(LEN=64)                         :: cvar_rad="radsw", cvar_wgt="radsw", cvar_out="radsw"
   CHARACTER(LEN=256)                        :: cmask="lsm_erainterim.nc"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()
   IF ( narg /= 2 ) THEN
      PRINT *,' Usage : correc_radsw_coef radsw_file coef_file '
      PRINT *,' lsm_erainterim.nc must be in your directory ' 
      STOP
   ENDIF

   CALL getarg (1, cfilein)
   CALL getarg (2, cwgtin )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Verify that input files have the same size

   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_radsw) )

   CALL check( NF90_INQ_DIMID(ncid_radsw, 'lon0', xid) )
   CALL check( NF90_INQ_DIMID(ncid_radsw, 'lat0', yid) )
   CALL check( NF90_INQ_DIMID(ncid_radsw, 'time', tid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_radsw, xid,len=npx) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_radsw, yid,len=npy) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_radsw, tid,len=npt) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in radsw file : npx = ', npx
   PRINT *, 'in radsw file : npy = ', npy
   PRINT *, 'in radsw file : npt = ', npt

   CALL check( NF90_CLOSE(ncid_radsw) )

   !!--------------------------------------------------------------
   CALL check( NF90_OPEN(trim(cwgtin), nf90_nowrite, ncid_wgt) )

   CALL check( NF90_INQ_DIMID(ncid_wgt, 'lon', xwid) )
   CALL check( NF90_INQ_DIMID(ncid_wgt, 'lat', ywid) )
   CALL check( NF90_INQ_DIMID(ncid_wgt, 'time', twid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, xwid,len=npxw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, ywid,len=npyw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, twid,len=nptw) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in weight file : npx = ', npxw
   PRINT *, 'in weight file : npy = ', npyw
   PRINT *, 'in weight file : npt = ', nptw
   PRINT *, '-------------------------------------------------------------'

   CALL check( NF90_CLOSE(ncid_wgt) )

   IF ( npx /= npxw .OR. npy /= npyw ) THEN
      PRINT *, 'Problem with dimensions in input files'
      STOP
   ENDIF

   !--------------------------------------------------------------------
   ALLOCATE( zradin(npx,npy) , zwgt(npx,npy), zradout(npx,npy) )
   ALLOCATE( zlon0(npx) , zlat0(npy) )
   ALLOCATE( zmask(npx,npy) )
   ALLOCATE( ztmp(npx,npy), ztim(npt) )

   count = (/ npx, npy, 1 /)
   start = (/ 1, 1, 1 /)

   !--------------------------------------------------------------------
   !! Get lon and lat
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_radsw) )
   CALL check( NF90_INQ_VARID(ncid_radsw, 'lon0', lonin_varid) )
   CALL check( NF90_INQ_VARID(ncid_radsw, 'lat0', latin_varid) )
   CALL check( NF90_GET_VAR(ncid_radsw, lonin_varid, zlon0 ) )
   CALL check( NF90_GET_VAR(ncid_radsw, latin_varid, zlat0 ) )
   CALL check( NF90_CLOSE(ncid_radsw) )
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Get the mask 

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
   CALL check( nf90_def_var(ncid_out, trim(cvar_out), NF90_REAL, (/ xoutid, youtid, toutid /), radout_varid) )

   ! Assign units attributes to coordinate variables.
   CALL check( nf90_put_att(ncid_out, latout_varid, "units", "degrees_north") )
   CALL check( nf90_put_att(ncid_out, lonout_varid, "units", "degrees_east")  )
   CALL check( nf90_put_att(ncid_out, radout_varid, "units", "W/m2") )
   CALL check( nf90_put_att(ncid_out, radout_varid, "long_name", "Shortwave radiation") )
   CALL check( nf90_put_att(ncid_out, tout_varid,   "units", "unknown") )


   ! End define mode.
   CALL check( nf90_enddef(ncid_out) )

   CALL check( nf90_put_var(ncid_out, latout_varid, zlat0) )
   CALL check( nf90_put_var(ncid_out, lonout_varid, zlon0) )

   DO irec=1,npt
      ztim(irec) = irec
   ENDDO

   CALL check( nf90_put_var(ncid_out, tout_varid, ztim) )


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read the netcdf files

   ! radsw
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_radsw) )
   CALL check( NF90_INQ_VARID(ncid_radsw, trim(cvar_rad), radin_varid) )
   ! weight
   CALL check( NF90_OPEN(trim(cwgtin), nf90_nowrite, ncid_wgt) )
   CALL check( NF90_INQ_VARID(ncid_wgt, trim(cvar_wgt), wgt_varid) )

   CALL check( NF90_GET_VAR(ncid_wgt, wgt_varid, zwgt ) )
   CALL check( NF90_CLOSE(ncid_wgt) )

   DO irec=1,npt

      PRINT *, 'Working on timestep : ', irec

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_radsw, radin_varid, zradin, start = start, &
                  &                 count = count) )


      !! this is my final array
      zradout(:,:) = zradin(:,:) * zwgt(:,:) * zmask(:,:)

      CALL check( NF90_PUT_VAR(ncid_out, radout_varid, REAL(zradout), start = start, &
                  &            count = count) )


   ENDDO

   CALL check( NF90_CLOSE(ncid_radsw) )
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
