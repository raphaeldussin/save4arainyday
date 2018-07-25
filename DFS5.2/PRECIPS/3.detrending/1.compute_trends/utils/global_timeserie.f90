PROGRAM global_timeserie

!!------------------------------------------------------!!
!!------------------------------------------------------!!
!!   Compute timeserie of variable from FOTO outputs    !!
!!           and outputs it in netcdf files             !!
!!                                                      !!
!!                  Part of FARC tools                  !!
!!                                                      !!
!!                Raphael Dussin, fev 2012              !!
!!                                                      !!
!!------------------------------------------------------!!
!!------------------------------------------------------!!

USE NETCDF

   IMPLICIT NONE

   INTEGER                                   :: narg, iargc                      ! command line arguments
   INTEGER                                   :: ncid_in, ncid_mask               ! netcdf id's
   INTEGER                                   :: ncid_ts, ncid_ym, ncid_gm        ! netcdf id's
   INTEGER                                   :: irec                             ! records loop index
   INTEGER                                   :: lonin_varid,  latin_varid
   INTEGER                                   :: lonin_dimid,  latin_dimid
   INTEGER                                   :: lonout_varid, latout_varid
   INTEGER                                   :: lonout_dimid, latout_dimid
   INTEGER                                   :: in_varid,     mask_varid
   INTEGER                                   :: timein_varid, timeout_varid
   INTEGER                                   :: timein_dimid, timeout_dimid
   INTEGER                                   :: ts_varid, ym_varid, gm_varid

   INTEGER                                   :: npx, npy, npt                     ! dimensions in rad and wgt
   INTEGER                                   :: xid, yid, tid                     ! netcdf dimensions id's

   INTEGER                                   :: start(3), count(3)

   INTEGER                                   :: jj, ji, ijarg, nfil, jt, jfr
   INTEGER                                   :: nftot=0
   INTEGER                                   :: ierr1, ierr2

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zvarin, zmask
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: zlon,   zlat 
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zlon2,  zlat2
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e1,     e2

   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: zts, zym, zgm, timeout, ztim

   REAL(KIND=8)                              :: dlon,dlat,zpi,rearth
   REAL(KIND=8)                              :: zfyear=0., zlyear=0. , zratio=1.

   REAL(KIND=8)                              :: rsf, rao ! scale factor / add offset

   CHARACTER(LEN=256)                            :: cfilein, cvar_in, cvar_out
   CHARACTER(LEN=256)                            :: cfileout_ts
   CHARACTER(LEN=256)                            :: cfileout_ym
   CHARACTER(LEN=256)                            :: cfileout_gm
   CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_list                      ! list of input files
   CHARACTER(LEN=256)                            :: cldum                        ! dummy string argument

   CHARACTER(LEN=256)                            :: cmask="mask.nc"  , cvar_mask="lsm"
   CHARACTER(LEN=256)                            :: clon="lon"       , clat="lat"
   CHARACTER(LEN=256)                            :: clonalt="lon0"   , clatalt="lat0"
   CHARACTER(LEN=256)                            :: cfreq, cfyear="0", clyear="0", carea, cdataset, cratio
   CHARACTER(LEN=256)                            :: cdiroutput


  !! Show usage if no args
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : global_timeserie list_of_files -var variable_name '
     PRINT *,'                                        -fyear first_year  '
     PRINT *,'                                        -lyear last_year   '
     PRINT *,'                                        -area area_of_interest '
     PRINT *,'                                        -dataset dataset_name '
     PRINT *,'                                        -diroutput output_directory '
     PRINT *,'                                        -mask path_to_mask '
     PRINT *,'                                        -ratio multiplicative ratio '
     PRINT *,'      '
     STOP
  ENDIF

  !! Set list of files
  ALLOCATE ( cf_list(narg) )
  ijarg = 1
  nfil = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( "-dataset" )
        CALL getarg (ijarg, cldum)
        cdataset = TRIM(cldum) ; ijarg = ijarg + 1
     CASE ( "-area" )
        CALL getarg (ijarg, cldum)
        carea = TRIM(cldum) ; ijarg = ijarg + 1
     CASE ( "-fyear" )
        CALL getarg (ijarg, cldum)
        cfyear = TRIM(cldum) ; READ(cfyear,*) zfyear ; ijarg = ijarg + 1
     CASE ( "-lyear" )
        CALL getarg (ijarg, cldum)
        clyear = TRIM(cldum) ; READ(clyear,*) zlyear ; ijarg = ijarg + 1
     CASE ( "-var" )
        CALL getarg (ijarg, cldum)
        cvar_in = TRIM(cldum) ; ijarg = ijarg + 1
     CASE ( "-diroutput" )
        CALL getarg (ijarg, cldum)
        cdiroutput = TRIM(cldum) ; ijarg = ijarg + 1
     CASE ( "-mask" )
        CALL getarg (ijarg, cldum)
        cmask = TRIM(cldum) ; ijarg = ijarg + 1
     CASE ( "-ratio" )
        CALL getarg (ijarg, cldum)
        cratio = TRIM(cldum) ; READ(cratio,*) zratio ; ijarg = ijarg + 1
     CASE DEFAULT         ! then the argument is a file
        nfil          = nfil + 1
        cf_list(nfil) = TRIM(cldum)
     END SELECT
  END DO

  PRINT *, '-----------------------------------------------------------'
  PRINT *, 'variable : ', TRIM(cvar_in) , ' , starting in ', TRIM(cfyear) , ' and ending in ', TRIM(clyear)
  PRINT *, 'Area = ', TRIM(carea)
  PRINT *, nfil, ' files to proceed'
  PRINT *, '-----------------------------------------------------------'

  !! We need to known the total number of frames
  DO ijarg=1,nfil

     CALL check( NF90_OPEN(trim(cf_list(ijarg)), nf90_nowrite, ncid_in) )
     CALL check( NF90_INQ_DIMID(ncid_in, 'time', tid) )
     CALL check( NF90_INQUIRE_DIMENSION(ncid_in, tid,len=npt) )
     CALL check( NF90_CLOSE(ncid_in) )

     nftot = nftot + npt

  ENDDO

  WRITE(*,'(a,i7.0,a)') , ' allocate for ', nftot , ' frames'

  ALLOCATE ( zts(nftot), zym(nftot), zgm(nftot), timeout(nftot) )
  ALLOCATE ( ztim(1) )

  !! Read Lon/Lat and Mask
  CALL check( NF90_OPEN(trim(cmask), nf90_nowrite, ncid_mask) )

  ! first we try with regular lon/lat
  ierr1 = NF90_INQ_DIMID(ncid_mask, clon, lonin_dimid)
  ierr2 = NF90_INQ_DIMID(ncid_mask, clat, latin_dimid)
  
  ! if it fails, try the alternative lon0/lat0
  IF ( (ierr1 /= NF90_NOERR).OR.(ierr2 /= NF90_NOERR) ) THEN
     CALL check( NF90_INQ_DIMID(ncid_mask, clonalt, lonin_dimid) )
     CALL check( NF90_INQ_DIMID(ncid_mask, clatalt, latin_dimid) )
  ENDIF

  CALL check( NF90_INQUIRE_DIMENSION(ncid_mask, lonin_dimid,len=npx) )
  CALL check( NF90_INQUIRE_DIMENSION(ncid_mask, latin_dimid,len=npy) )

  ALLOCATE ( zlon(npx) , zlat(npy) )
  ALLOCATE ( zlon2(npx,npy) , zlat2(npx,npy) )
  ALLOCATE ( zmask(npx,npy) )

  start = (/ 1, 1, 1 /)
  count = (/ npx, npy, 1 /)

  ! first we try with regular lon/lat
  ierr1 = NF90_INQ_VARID(ncid_mask, clon, lonin_varid)
  ierr2 = NF90_INQ_VARID(ncid_mask, clat, latin_varid)
  
  ! if it fails, try the alternative lon0/lat0
  IF ( (ierr1 /= NF90_NOERR).OR.(ierr2 /= NF90_NOERR) ) THEN
     CALL check( NF90_INQ_VARID(ncid_mask, clonalt, lonin_varid) )
     CALL check( NF90_INQ_VARID(ncid_mask, clatalt, latin_varid) )
  ENDIF

  CALL check( NF90_INQ_VARID(ncid_mask, cvar_mask, mask_varid) )

  CALL check( NF90_GET_VAR(ncid_mask, lonin_varid, zlon, start = (/1/), &
             &                 count = (/npx/) ) )
  CALL check( NF90_GET_VAR(ncid_mask, latin_varid, zlat, start = (/1/), &
             &                 count = (/npy/) ) )
  CALL check( NF90_GET_VAR(ncid_mask, mask_varid, zmask, start = start, &
             &                 count = count) )
  CALL check( NF90_CLOSE(ncid_mask) )

  !! Create regular metrics
  zpi=ACOS(-1.)
  rearth=6371229.

  ALLOCATE( e1(npx,npy) , e2(npx,npy) )
  ALLOCATE( zvarin(npx,npy) )

  e1 = -1.
  e2 = -1.

  DO jj=2,npy-1
     DO ji=2,npx+1
     
        dlon = 0.5 * ABS( zlon(ji+1) - zlon(ji-1) ) * zpi / 180.
        dlat = 0.5 * ABS( zlat(jj+1) - zlat(jj-1) ) * zpi / 180.
        e1(ji,jj) = rearth * dlon * COS( zpi * zlat(jj) / 180. )
        e2(ji,jj) = rearth * dlat

     ENDDO
  ENDDO

  !! boundaries
  e1(:,1)   = e1(:,2)
  e1(:,npy) = e1(:,npy-1)
  e1(1,:)   = e1(2,:)
  e1(npx,:) = e1(npx-1,:)

  e2(:,1)   = e2(:,2)
  e2(:,npy) = e2(:,npy-1)
  e2(1,:)   = e2(2,:)
  e2(npx,:) = e2(npx-1,:)


  !! Create 2d array for lon and lat
  DO ji=1,npx
     zlat2(ji,:) = zlat
  ENDDO

  DO jj=1,npy
     zlon2(:,jj) = zlon
  ENDDO

  !! We set e1 to zero where we want to mask
  SELECT CASE ( carea )

     CASE ('Polar_South')
          WHERE( zlat2 > -70. ) e1(:,:) = 0.

     CASE ('Subpolar_South')
          WHERE( zlat2 < -70. ) e1(:,:) = 0.
          WHERE( zlat2 > -45. ) e1(:,:) = 0.

     CASE ('Subtropical_South')
          WHERE( zlat2 < -45. ) e1(:,:) = 0.
          WHERE( zlat2 > -25. ) e1(:,:) = 0.

     CASE ('Tropical_South')
          WHERE( zlat2 < -25. ) e1(:,:) = 0.
          WHERE( zlat2 > -10. ) e1(:,:) = 0.

     CASE ('Equatorial_Band')
          WHERE( zlat2 < -10. ) e1(:,:) = 0.
          WHERE( zlat2 >  10. ) e1(:,:) = 0.

     CASE ('Tropical_North')
          WHERE( zlat2 < 10. ) e1(:,:) = 0.
          WHERE( zlat2 > 25. ) e1(:,:) = 0.
     
     CASE ('Subtropical_North')
          WHERE( zlat2 < 25. ) e1(:,:) = 0.
          WHERE( zlat2 > 45. ) e1(:,:) = 0.

     CASE ('Subpolar_North')
          WHERE( zlat2 < 45. ) e1(:,:) = 0.
          WHERE( zlat2 > 70. ) e1(:,:) = 0.

     CASE ('Polar_North')
          WHERE( zlat2 < 70. ) e1(:,:) = 0.

     CASE ('global')
          ! Nothing to do
     CASE DEFAULT
          PRINT *, 'Uncorrect area'; STOP
  END SELECT

  !! We start the computation
  jfr = 0

  DO ijarg=1,nfil

     PRINT *, 'working on '//TRIM(cf_list(ijarg))
     ! get number of frames
     CALL check( NF90_OPEN(trim(cf_list(ijarg)), nf90_nowrite, ncid_in) )
     CALL check( NF90_INQ_DIMID(ncid_in, 'time', timein_dimid) )
     CALL check( NF90_INQUIRE_DIMENSION(ncid_in, timein_dimid,len=npt) )

     CALL check( NF90_INQ_VARID(ncid_in, TRIM(cvar_in), in_varid) )  

     ierr1 = NF90_GET_ATT(ncid_in, in_varid, 'scale_factor', rsf)
     ierr2 = NF90_GET_ATT(ncid_in, in_varid, 'add_offset',   rao)
     !!
     IF ( (ierr1 /= NF90_NOERR).OR.(ierr2 /= NF90_NOERR) ) THEN
        rsf = 1.      ;   rao = 0.
     ENDIF

     DO jt=1,npt

        jfr = jfr + 1

        start = (/ 1, 1, jt /)
        count = (/ npx, npy, 1 /)

        CALL check( NF90_GET_VAR(ncid_in, in_varid, zvarin, start = start, &
                &                 count = count ) )

        ! scale factor / offset
        zvarin = (rsf * zvarin) + rao
        zvarin = zratio * zvarin

        zts(jfr) = sum( zvarin * e1 * e2 * zmask ) / sum( e1 * e2 * zmask )
        timeout(jfr) = zfyear + (ijarg-1) + (1. * jt / npt )

     ENDDO

     CALL check( NF90_CLOSE(ncid_in) )

     zym(jfr-npt+1:jfr) = sum( zts(jfr-npt+1:jfr) ) / npt

  ENDDO

  zgm(:) = sum( zts(:) ) / nftot


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  create the netcdf output file

cvar_out = TRIM(cvar_in)//"_"//TRIM(cdataset)
cfileout_ts=TRIM(carea)//"_"//TRIM(cvar_in)//"_Dataset_"//TRIM(cdataset)//"_"//TRIM(cfyear)//"-"//TRIM(clyear)//"_full_timeserie.nc"
cfileout_ym=TRIM(carea)//"_"//TRIM(cvar_in)//"_Dataset_"//TRIM(cdataset)//"_"//TRIM(cfyear)//"-"//TRIM(clyear)//"_yearly_timeserie.nc"
cfileout_gm=TRIM(carea)//"_"//TRIM(cvar_in)//"_Dataset_"//TRIM(cdataset)//"_"//TRIM(cfyear)//"-"//TRIM(clyear)//"_mean_of_timeserie.nc"

cfileout_ts=TRIM(cdiroutput)//"/"//TRIM(cfileout_ts)
cfileout_ym=TRIM(cdiroutput)//"/"//TRIM(cfileout_ym)
cfileout_gm=TRIM(cdiroutput)//"/"//TRIM(cfileout_gm)

  !! timeserie
  CALL check( nf90_create(TRIM(cfileout_ts), nf90_clobber, ncid_ts) )
  CALL check( nf90_def_dim(ncid_ts, 'time', NF90_UNLIMITED, timeout_dimid) )
  CALL check( nf90_def_var(ncid_ts, 'time', NF90_REAL, (/ timeout_dimid /), timeout_varid) )
  CALL check( nf90_def_var(ncid_ts, trim(cvar_out), NF90_DOUBLE, (/ timeout_dimid /), ts_varid) )
  CALL check( nf90_enddef(ncid_ts) )

  CALL check( nf90_put_var(ncid_ts, timeout_varid, timeout) )
  CALL check( nf90_put_var(ncid_ts, ts_varid, zts) )
  CALL check( nf90_close(ncid_ts) )

  !! yearly means
  CALL check( nf90_create(TRIM(cfileout_ym), nf90_clobber, ncid_ym) )
  CALL check( nf90_def_dim(ncid_ym, 'time', NF90_UNLIMITED, timeout_dimid) )
  CALL check( nf90_def_var(ncid_ym, 'time', NF90_REAL, (/ timeout_dimid /), timeout_varid) )
  CALL check( nf90_def_var(ncid_ym, trim(cvar_out), NF90_DOUBLE, (/ timeout_dimid /), ym_varid) )
  CALL check( nf90_enddef(ncid_ym) )

  CALL check( nf90_put_var(ncid_ym, timeout_varid, timeout) )
  CALL check( nf90_put_var(ncid_ym, ym_varid, zym) )
  CALL check( nf90_close(ncid_ym) )

  !! global mean
  CALL check( nf90_create(TRIM(cfileout_gm), nf90_clobber, ncid_gm) )
  CALL check( nf90_def_dim(ncid_gm, 'time', NF90_UNLIMITED, timeout_dimid) )
  CALL check( nf90_def_var(ncid_gm, 'time', NF90_REAL, (/ timeout_dimid /), timeout_varid) )
  CALL check( nf90_def_var(ncid_gm, trim(cvar_out), NF90_DOUBLE, (/ timeout_dimid /), gm_varid) )
  CALL check( nf90_enddef(ncid_gm) )

  CALL check( nf90_put_var(ncid_gm, timeout_varid, timeout) )
  CALL check( nf90_put_var(ncid_gm, gm_varid, zgm) )
  CALL check( nf90_close(ncid_gm) )


  PRINT *, 'done ! '

CONTAINS

      SUBROUTINE check(status)
        INTEGER, INTENT(in) :: status

        IF(status /= nf90_noerr) THEN
          PRINT *, trim(nf90_strerror(status))
          STOP 2
        ENDIF
      END SUBROUTINE check


END PROGRAM global_timeserie
