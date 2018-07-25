PROGRAM replace_meanfield

!------------------------------------------------------
! Replaces the time-averaged mean field by an external
! field given in input
! 
!
! Raphael Dussin, mar 2011
!
!------------------------------------------------------

USE NETCDF

   IMPLICIT NONE

   INTEGER                                   :: narg, iargc                     ! command line arguments
   INTEGER                                   :: ncid_forcin, ncid_mean          ! netcdf id's
   INTEGER                                   :: ncid_forcout                    ! netcdf id's
   INTEGER                                   :: irec                            ! records loop index
   INTEGER                                   :: lonin_varid,  latin_varid       ! variables id
   INTEGER                                   :: lonout_varid, latout_varid      ! variables id
   INTEGER                                   :: forcin_varid, forcout_varid     ! variables id
   INTEGER                                   :: mean_varid                      ! variables id

   INTEGER                                   :: npx, npy, npt, npxm, npym, nptm ! dimensions in input
   INTEGER                                   :: xid, yid, tid, xmid, ymid, tmid ! netcdf dimensions id's
   INTEGER                                   :: xoutid, youtid, toutid

   INTEGER                                   :: start(3), count(3)
   INTEGER                                   :: startlon(2), countlon(2)
   INTEGER                                   :: startlat(2), countlat(2)

   INTEGER                                   :: jj, ji

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zforcin, zmean, zforcout
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: zlon, zlat
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: znewmean

   CHARACTER(LEN=256)                        :: cfilein, cmeanin
   CHARACTER(LEN=256), PARAMETER             :: cfileout="output.nc"
   CHARACTER(LEN=64)                         :: cvar_in, cvar_out

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()
   IF ( narg /= 3 ) THEN
      PRINT *,' Usage : replace_meanfield variable forcing_file mean_file '
      PRINT *,' output on output.nc                                       '
      PRINT *,' --------------------------------------------------------- ' 
      STOP
   ENDIF

   CALL getarg (1, cvar_in)
   CALL getarg (2, cfilein)
   CALL getarg (3, cmeanin)

   cvar_out = cvar_in

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Verify that input files have the same size

   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_forcin) )

   CALL check( NF90_INQ_DIMID(ncid_forcin, 'lon', xid) )
   CALL check( NF90_INQ_DIMID(ncid_forcin, 'lat', yid) )
   CALL check( NF90_INQ_DIMID(ncid_forcin, 'time', tid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_forcin, xid,len=npx) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_forcin, yid,len=npy) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_forcin, tid,len=npt) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in forcing file : npx = ', npx
   PRINT *, 'in forcing file : npy = ', npy
   PRINT *, 'in forcing file : npt = ', npt

   CALL check( NF90_CLOSE(ncid_forcin) )

   !!--------------------------------------------------------------
   CALL check( NF90_OPEN(trim(cmeanin), nf90_nowrite, ncid_mean) )

   CALL check( NF90_INQ_DIMID(ncid_mean, 'lon', xmid) )
   CALL check( NF90_INQ_DIMID(ncid_mean, 'lat', ymid) )
   CALL check( NF90_INQ_DIMID(ncid_mean, 'time', tmid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_mean, xmid,len=npxm) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_mean, ymid,len=npym) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_mean, tmid,len=nptm) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in mean file : npx = ', npxm
   PRINT *, 'in mean file : npy = ', npym
   PRINT *, 'in mean file : npt = ', nptm
   PRINT *, '-------------------------------------------------------------'

   CALL check( NF90_CLOSE(ncid_mean) )

   IF ( npx /= npxm .OR. npy /= npym ) THEN
      PRINT *, 'Problem with dimensions in input files'
      STOP
   ENDIF

   IF ( nptm /= 1 ) THEN
      PRINT *, 'Your mean file must contain only one record'
      STOP
   ENDIF

   IF ( npt /= 12 .AND. npt /= 365 .AND. npt /= 366 ) THEN
      PRINT *, 'Your forcing file time dimension looks odd'
      STOP
   ENDIF

   !--------------------------------------------------------------------
   ALLOCATE( zforcin(npx,npy) , zmean(npx,npy), zforcout(npx,npy) )
   ALLOCATE( zlon(npx) , zlat(npy) )
   ALLOCATE( znewmean(npx,npy) )

   count = (/ npx, npy, 1 /)
   start = (/ 1, 1, 1 /)
   
   countlon = (/ npx, 1 /)
   startlon = (/ 1, 1 /)

   countlat = (/ npy, 1 /)
   startlat = (/ 1, 1 /)

   !--------------------------------------------------------------------
   !! Get lon and lat
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_forcin) )
   CALL check( NF90_INQ_VARID(ncid_forcin, 'lon', lonin_varid) )
   CALL check( NF90_INQ_VARID(ncid_forcin, 'lat', latin_varid) )
   CALL check( NF90_GET_VAR(ncid_forcin, lonin_varid, zlon, start = startlon, &
              &                 count = countlon) )
   CALL check( NF90_GET_VAR(ncid_forcin, latin_varid, zlat, start = startlat, &
              &                 count = countlat) )

   CALL check( NF90_CLOSE(ncid_forcin) )
   
   PRINT *, 'Max of lon is ', MAXVAL(zlon)
   PRINT *, 'Max of lat is ', MAXVAL(zlat)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  create the netcdf output file

   CALL check( nf90_create(trim(cfileout), nf90_clobber, ncid_forcout) )

   CALL check( nf90_def_dim(ncid_forcout, 'lon', npx, xoutid) )
   CALL check( nf90_def_dim(ncid_forcout, 'lat', npy, youtid) )
   CALL check( nf90_def_dim(ncid_forcout, 'time', NF90_UNLIMITED, toutid) )

   CALL check( nf90_def_var(ncid_forcout, 'lon', NF90_REAL, (/ xoutid /), lonout_varid) )
   CALL check( nf90_def_var(ncid_forcout, 'lat', NF90_REAL, (/ youtid /), latout_varid) )
   CALL check( nf90_def_var(ncid_forcout, trim(cvar_out), NF90_REAL, (/ xoutid, youtid, toutid /), forcout_varid) )

   ! Assign units attributes to coordinate variables.
!   CALL check( nf90_put_att(ncid_forcout, latout_varid, "units", "degrees_north") )
!   CALL check( nf90_put_att(ncid_forcout, lonout_varid, "units", "degrees_east")  )
!   CALL check( nf90_put_att(ncid_forcout, forcout_varid, "units", "kg.m-2.s-1") )


   ! End define mode.
   CALL check( nf90_enddef(ncid_forcout) )

   CALL check( nf90_put_var(ncid_forcout, latout_varid, zlat) )
   CALL check( nf90_put_var(ncid_forcout, lonout_varid, zlon) )


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read the netcdf files

   ! forcing input
   CALL check( NF90_OPEN(trim(cfilein), nf90_nowrite, ncid_forcin) )
   CALL check( NF90_INQ_VARID(ncid_forcin, trim(cvar_in), forcin_varid) )
   ! new mean
   CALL check( NF90_OPEN(trim(cmeanin), nf90_nowrite, ncid_mean) )
   CALL check( NF90_INQ_VARID(ncid_mean, trim(cvar_in), mean_varid) )

   !!----------------------------------------------------------------------------
   !! Get the new mean for the forcing field
   CALL check( NF90_GET_VAR(ncid_mean, mean_varid, znewmean, start = start, &
              &                 count = count) )

   !!----------------------------------------------------------------------------
   !! computing the old mean
   DO irec=1,npt

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_forcin, forcin_varid, zforcin, start = start, &
                  &                 count = count) )

      zmean(:,:) = zmean(:,:) + zforcin(:,:)


   ENDDO

      zmean(:,:) = zmean(:,:) / npt

   !!----------------------------------------------------------------------------
   !! replace in each frame
   DO irec=1,npt

      PRINT *, 'Working on timestep : ', irec

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_forcin, forcin_varid, zforcin, start = start, &
                  &                 count = count) )

      zforcout(:,:) = zforcin(:,:) + znewmean(:,:) - zmean(:,:)

      CALL check( NF90_PUT_VAR(ncid_forcout, forcout_varid, REAL(zforcout), start = start, &
                  &            count = count) )


   ENDDO

   CALL check( NF90_CLOSE(ncid_forcin) )
   CALL check( NF90_CLOSE(ncid_mean) )
   CALL check( nf90_close(ncid_forcout) )


CONTAINS

      SUBROUTINE check(status)
        INTEGER, INTENT(in) :: status
    
        IF(status /= nf90_noerr) THEN
          PRINT *, trim(nf90_strerror(status))
          STOP 2
        ENDIF
      END SUBROUTINE check


END PROGRAM
