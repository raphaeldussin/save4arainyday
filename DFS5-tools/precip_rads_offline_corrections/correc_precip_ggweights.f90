PROGRAM correc_precip_ggweights

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
   INTEGER                                   :: ncid_precip, ncid_wgt, ncid_out  ! netcdf id's
   INTEGER                                   :: ncid_mesh, ncid_mask
   INTEGER                                   :: irec                            ! records loop index
   INTEGER                                   :: lonin_varid, latin_varid
   INTEGER                                   :: lonout_varid, latout_varid
   INTEGER                                   :: prein_varid, preout_varid
   INTEGER                                   :: wgt_varid, mask_varid, phi_varid

   INTEGER                                   :: npx, npy, npt, npxw, npyw, nptw ! dimensions in rad and wgt
   INTEGER                                   :: xid, yid, tid, xwid, ywid, twid ! netcdf dimensions id's
   INTEGER                                   :: xoutid, youtid, toutid

   INTEGER                                   :: start(3), count(3)

   INTEGER                                   :: iter_shapiro=250
   INTEGER                                   :: jj, ji

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zprein, zwgt, zpreout
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlon, zlat
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zmask, zphi, zxyt
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zprein_lr, zprein_hr, ztmp, ztmp2

   REAL(KIND=8)                              :: zzlat, zzlat1, zzlat2

   CHARACTER(LEN=256)                        :: cfilein, cwgtin
   CHARACTER(LEN=256), PARAMETER             :: cfileout="precip_corrected.nc"
   CHARACTER(LEN=64)                         :: cvar_rad="precip", cvar_wgt="socoprec", cvar_out="precip"
   CHARACTER(LEN=256)                        :: cmask="mask.nc", cmesh="mesh_hgr.nc"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()
   IF ( narg /= 2 ) THEN
      PRINT *,' Usage : correc_precip_ggweights precip_file weight_file '
      PRINT *,' mesh_hgr.nc and mask.nc must be in your directory ' 
      STOP
   ENDIF

   CALL getarg (1, cfilein)
   CALL getarg (2, cwgtin )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Verify that input files have the same size

   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )

   CALL check( NF90_INQ_DIMID(ncid_precip, 'x', xid) )
   CALL check( NF90_INQ_DIMID(ncid_precip, 'y', yid) )
   CALL check( NF90_INQ_DIMID(ncid_precip, 'time_counter', tid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, xid,len=npx) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, yid,len=npy) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_precip, tid,len=npt) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in radsw file : npx = ', npx
   PRINT *, 'in radsw file : npy = ', npy
   PRINT *, 'in radsw file : npt = ', npt

   CALL check( NF90_CLOSE(ncid_precip) )

   !!--------------------------------------------------------------
   CALL check( NF90_OPEN(trim(cwgtin), nf90_nowrite, ncid_wgt) )

   CALL check( NF90_INQ_DIMID(ncid_wgt, 'x', xwid) )
   CALL check( NF90_INQ_DIMID(ncid_wgt, 'y', ywid) )
   CALL check( NF90_INQ_DIMID(ncid_wgt, 'time_counter', twid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, xwid,len=npxw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, ywid,len=npyw) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_wgt, twid,len=nptw) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in weight file : npx = ', npxw
   PRINT *, 'in weight file : npy = ', npyw
   PRINT *, 'in weight file : npt = ', nptw
   PRINT *, '-------------------------------------------------------------'

   CALL check( NF90_CLOSE(ncid_wgt) )

   IF ( npx /= npxw .OR. npy /= npyw .OR. npt /= nptw ) THEN
      PRINT *, 'Problem with dimensions in input files'
      STOP
   ENDIF

   !--------------------------------------------------------------------
   ALLOCATE( zprein(npx,npy) , zwgt(npx,npy), zpreout(npx,npy) )
   ALLOCATE( zlon(npx,npy) , zlat(npx,npy) )
   ALLOCATE( zmask(npx,npy), zphi(npx,npy), zxyt(npx,npy) )
   ALLOCATE( zprein_lr(npx,npy), zprein_hr(npx,npy) )
   ALLOCATE( ztmp(npx,npy), ztmp2(npx,npy) )

   count = (/ npx, npy, 1 /)
   start = (/ 1, 1, 1 /)

   !--------------------------------------------------------------------
   !! Get lon and lat
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )
   CALL check( NF90_INQ_VARID(ncid_precip, 'nav_lon', lonin_varid) )
   CALL check( NF90_INQ_VARID(ncid_precip, 'nav_lat', latin_varid) )
   CALL check( NF90_GET_VAR(ncid_precip, lonin_varid, zlon, start = start, &
              &                 count = count) )
   CALL check( NF90_GET_VAR(ncid_precip, latin_varid, zlat, start = start, &
              &                 count = count) )
   CALL check( NF90_CLOSE(ncid_precip) )
   
!   PRINT *, 'Max of lon is ', MAXVAL(zlon)
!   PRINT *, 'Max of lat is ', MAXVAL(zlat)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Get the mask and gphi from mask.nc and mesh_hgr.nc

   CALL check( NF90_OPEN(trim(cmask), nf90_nowrite, ncid_mask) )
   CALL check( NF90_INQ_VARID(ncid_mask, 'tmask', mask_varid) )
   CALL check( NF90_GET_VAR(ncid_mask, mask_varid, zmask, start = start, &
              &                 count = count) )
   CALL check( NF90_CLOSE(ncid_mask) )

   CALL check( NF90_OPEN(trim(cmesh), nf90_nowrite, ncid_mesh) )
   CALL check( NF90_INQ_VARID(ncid_mesh, 'gphit', phi_varid) )
   CALL check( NF90_GET_VAR(ncid_mesh, phi_varid, zphi, start = start, &
              &                 count = count) )
   CALL check( NF90_CLOSE(ncid_mesh) )

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  create the netcdf output file

   CALL check( nf90_create(trim(cfileout), nf90_clobber, ncid_out) )

   CALL check( nf90_def_dim(ncid_out, 'x', npx, xoutid) )
   CALL check( nf90_def_dim(ncid_out, 'y', npy, youtid) )
   CALL check( nf90_def_dim(ncid_out, 'time_counter', NF90_UNLIMITED, toutid) )

   CALL check( nf90_def_var(ncid_out, 'nav_lon', NF90_REAL, (/ xoutid, youtid /), lonout_varid) )
   CALL check( nf90_def_var(ncid_out, 'nav_lat', NF90_REAL, (/ xoutid, youtid /), latout_varid) )
   CALL check( nf90_def_var(ncid_out, trim(cvar_out), NF90_REAL, (/ xoutid, youtid, toutid /), preout_varid) )

   ! Assign units attributes to coordinate variables.
   CALL check( nf90_put_att(ncid_out, latout_varid, "units", "degrees_north") )
   CALL check( nf90_put_att(ncid_out, lonout_varid, "units", "degrees_east")  )
   CALL check( nf90_put_att(ncid_out, preout_varid, "units", "kg.m-2.s-1") )


   ! End define mode.
   CALL check( nf90_enddef(ncid_out) )

   CALL check( nf90_put_var(ncid_out, latout_varid, zlat) )
   CALL check( nf90_put_var(ncid_out, lonout_varid, zlon) )


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read the netcdf files

   ! radsw
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_precip) )
   CALL check( NF90_INQ_VARID(ncid_precip, trim(cvar_rad), prein_varid) )
   ! weight
   CALL check( NF90_OPEN(trim(cwgtin), nf90_nowrite, ncid_wgt) )
   CALL check( NF90_INQ_VARID(ncid_wgt, trim(cvar_wgt), wgt_varid) )

   DO irec=1,npt

      PRINT *, 'Working on timestep : ', irec

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_precip, prein_varid, zprein, start = start, &
                  &                 count = count) )

      CALL check( NF90_GET_VAR(ncid_wgt, wgt_varid, zwgt, start = start, &
                  &                 count = count) )

!      PRINT *, '-----------------------------------------------------------------------'
!      PRINT *, 'At record ', irec, 'Max of ', trim(cvar_rad), ' is ', MAXVAL(zprein)
!      PRINT *, 'At record ', irec, 'Min of ', trim(cvar_rad), ' is ', MINVAL(zprein)
!      PRINT *, 'At record ', irec, 'Max of ', trim(cvar_wgt), ' is ', MAXVAL(zwgt)
!      PRINT *, 'At record ', irec, 'Min of ', trim(cvar_wgt), ' is ', MINVAL(zwgt)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Split into large and small scales

      CALL Shapiro_1D(zprein(:,:),iter_shapiro, .TRUE., zprein_lr)
      zprein_hr(:,:) = zprein(:,:) - zprein_lr(:,:)
      ztmp(:,:)      = zprein_lr(:,:) * zwgt(:,:) + zprein_hr(:,:) 

      !! again some thresholding for the pleasure of wasting CPU time
      DO jj=1,npy
        DO ji=1,npx
           ztmp(ji,jj)=max(ztmp(ji,jj),0.0)
        END DO
      END DO

      !! this is my final array
      zpreout(:,:) = ztmp(:,:) * zmask(:,:)

      CALL check( NF90_PUT_VAR(ncid_out, preout_varid, REAL(zpreout), start = start, &
                  &            count = count) )


   ENDDO

   CALL check( NF90_CLOSE(ncid_precip) )
   CALL check( NF90_CLOSE(ncid_wgt) )
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
