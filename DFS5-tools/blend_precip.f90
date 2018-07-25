PROGRAM replace_pre2field

!------------------------------------------------------
! Replaces the time-averaged pre2 field by an external
! field given in input
! 
!
! Raphael Dussin, mar 2011
!
!------------------------------------------------------

USE NETCDF

   IMPLICIT NONE

   INTEGER                                   :: narg, iargc                     ! command line arguments
   INTEGER                                   :: ncid_pre1, ncid_pre2          ! netcdf id's
   INTEGER                                   :: ncid_preout                    ! netcdf id's
   INTEGER                                   :: irec                            ! records loop index
   INTEGER                                   :: lonin_varid,  latin_varid       ! variables id
   INTEGER                                   :: lonout_varid, latout_varid      ! variables id
   INTEGER                                   :: pre1_varid, preout_varid     ! variables id
   INTEGER                                   :: pre2_varid                      ! variables id

   INTEGER                                   :: npx, npy, npt, npxm, npym, nptm ! dimensions in input
   INTEGER                                   :: xid, yid, tid, xmid, ymid, tmid ! netcdf dimensions id's
   INTEGER                                   :: xoutid, youtid, toutid

   INTEGER                                   :: start(3), count(3)
   INTEGER                                   :: startlon(2), countlon(2)
   INTEGER                                   :: startlat(2), countlat(2)

   INTEGER                                   :: jj, ji

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zpre1, zpre2, zpreout
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ztmp, ztmp2
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: zlon, zlat

   INTEGER                                   :: zlat30s, zlat25s, zlat30n, zlat25n
   REAL(KIND=4)                              :: zc1, zc2

   CHARACTER(LEN=256)                        :: cfilepre1, cfilepre2
   CHARACTER(LEN=256), PARAMETER             :: cfileout="output.nc"
   CHARACTER(LEN=64)                         :: cvar_in="precip", cvar_out="precip"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read command line
   narg= iargc()
   IF ( narg /= 2 ) THEN
      PRINT *,' Usage : blend_precip lowlat-file highlat-file             '
      PRINT *,' output on output.nc                                       '
      PRINT *,' --------------------------------------------------------- ' 
      STOP
   ENDIF

   CALL getarg (1, cfilepre1)
   CALL getarg (2, cfilepre2)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Verify that input files have the same size

   CALL check( NF90_OPEN(trim(cfilepre1), nf90_nowrite, ncid_pre1) )

   CALL check( NF90_INQ_DIMID(ncid_pre1, 'lon', xid) )
   CALL check( NF90_INQ_DIMID(ncid_pre1, 'lat', yid) )
   CALL check( NF90_INQ_DIMID(ncid_pre1, 'time', tid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre1, xid,len=npx) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre1, yid,len=npy) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre1, tid,len=npt) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in pre1g file : npx = ', npx
   PRINT *, 'in pre1g file : npy = ', npy
   PRINT *, 'in pre1g file : npt = ', npt

   CALL check( NF90_CLOSE(ncid_pre1) )

   !!--------------------------------------------------------------
   CALL check( NF90_OPEN(trim(cfilepre2), nf90_nowrite, ncid_pre2) )

   CALL check( NF90_INQ_DIMID(ncid_pre2, 'lon0', xmid) )
   CALL check( NF90_INQ_DIMID(ncid_pre2, 'lat0', ymid) )
   CALL check( NF90_INQ_DIMID(ncid_pre2, 'time', tmid) )

   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre2, xmid,len=npxm) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre2, ymid,len=npym) )
   CALL check( NF90_INQUIRE_DIMENSION(ncid_pre2, tmid,len=nptm) )

   PRINT *, '-------------------------------------------------------------'
   PRINT *, 'in pre2 file : npx = ', npxm
   PRINT *, 'in pre2 file : npy = ', npym
   PRINT *, 'in pre2 file : npt = ', nptm
   PRINT *, '-------------------------------------------------------------'

   CALL check( NF90_CLOSE(ncid_pre2) )

   IF ( npx /= npxm .OR. npy /= npym ) THEN
      PRINT *, 'Problem with dimensions in input files'
      STOP
   ENDIF

   IF ( npt /= 12 .AND. npt /= 365 .AND. npt /= 366 ) THEN
      PRINT *, 'Your pre1g file time dimension looks odd'
      !STOP
   ENDIF

   !--------------------------------------------------------------------
   ALLOCATE( zpre1(npx,npy) , zpre2(npx,npy), zpreout(npx,npy) )
   ALLOCATE( zlon(npx) , zlat(npy) )
   ALLOCATE( ztmp(npx,npy), ztmp2(npx,npy) )

   count = (/ npx, npy, 1 /)
   start = (/ 1, 1, 1 /)
   
   countlon = (/ npx, 1 /)
   startlon = (/ 1, 1 /)

   countlat = (/ npy, 1 /)
   startlat = (/ 1, 1 /)

   !--------------------------------------------------------------------
   !! Get lon and lat
   CALL check( NF90_OPEN(trim(cfilepre1), nf90_nowrite, ncid_pre1) )
   CALL check( NF90_INQ_VARID(ncid_pre1, 'lon', lonin_varid) )
   CALL check( NF90_INQ_VARID(ncid_pre1, 'lat', latin_varid) )
   CALL check( NF90_GET_VAR(ncid_pre1, lonin_varid, zlon, start = startlon, &
              &                 count = countlon) )
   CALL check( NF90_GET_VAR(ncid_pre1, latin_varid, zlat, start = startlat, &
              &                 count = countlat) )

   CALL check( NF90_CLOSE(ncid_pre1) )
   
   PRINT *, 'Max of lon is ', MAXVAL(zlon)
   PRINT *, 'Max of lat is ', MAXVAL(zlat)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  create the netcdf output file

   CALL check( nf90_create(trim(cfileout), nf90_clobber, ncid_preout) )

   CALL check( nf90_def_dim(ncid_preout, 'lon', npx, xoutid) )
   CALL check( nf90_def_dim(ncid_preout, 'lat', npy, youtid) )
   CALL check( nf90_def_dim(ncid_preout, 'time', NF90_UNLIMITED, toutid) )

   CALL check( nf90_def_var(ncid_preout, 'lon', NF90_REAL, (/ xoutid /), lonout_varid) )
   CALL check( nf90_def_var(ncid_preout, 'lat', NF90_REAL, (/ youtid /), latout_varid) )
   CALL check( nf90_def_var(ncid_preout, trim(cvar_out), NF90_REAL, (/ xoutid, youtid, toutid /), preout_varid) )

   ! Assign units attributes to coordinate variables.
!   CALL check( nf90_put_att(ncid_preout, latout_varid, "units", "degrees_north") )
!   CALL check( nf90_put_att(ncid_preout, lonout_varid, "units", "degrees_east")  )
!   CALL check( nf90_put_att(ncid_preout, preout_varid, "units", "kg.m-2.s-1") )


   ! End define mode.
   CALL check( nf90_enddef(ncid_preout) )

   CALL check( nf90_put_var(ncid_preout, latout_varid, zlat) )
   CALL check( nf90_put_var(ncid_preout, lonout_varid, zlon) )


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!  Read the netcdf files

   ! lowlat file
   CALL check( NF90_OPEN(trim(cfilepre1), nf90_nowrite, ncid_pre1) )
   CALL check( NF90_INQ_VARID(ncid_pre1, trim(cvar_in), pre1_varid) )
   ! highlat file
   CALL check( NF90_OPEN(trim(cfilepre2), nf90_nowrite, ncid_pre2) )
   CALL check( NF90_INQ_VARID(ncid_pre2, trim(cvar_in), pre2_varid) )

   !!----------------------------------------------------------------------------
   !! computing the transition lat

   zlat30s = MINLOC(ABS(zlat+30),1) ; zlat30n = MINLOC(ABS(zlat-30),1)
   zlat25s = MINLOC(ABS(zlat+25),1) ; zlat25n = MINLOC(ABS(zlat-25),1)

   PRINT *, '30S is point ', zlat30s ,' = ', zlat(zlat30s)
   PRINT *, '25S is point ', zlat25s ,' = ', zlat(zlat25s)
   PRINT *, '25N is point ', zlat25n ,' = ', zlat(zlat25n)
   PRINT *, '30N is point ', zlat30n ,' = ', zlat(zlat30n)

   DO irec=1,npt

      zpreout(:,:) = 0.

      start(3) = irec
      CALL check( NF90_GET_VAR(ncid_pre1, pre1_varid, zpre1, start = start, &
                  &                 count = count) )

      CALL check( NF90_GET_VAR(ncid_pre2, pre2_varid, zpre2, start = start, &
                  &                 count = count) )

      ! fill low latitude
      zpreout(:,zlat25n:zlat25s) = zpre1(:,zlat25n:zlat25s)
      ! fill high latitude
      zpreout(:,1:zlat30n)   = zpre2(:,1:zlat30n)
      zpreout(:,zlat30s:npy) = zpre2(:,zlat30s:npy)
      ! fill northen transition zone
      DO jj=zlat30n,zlat25n
         zc1 = FLOAT((jj - zlat30n)) / FLOAT((zlat25n - zlat30n))
         zc2 = 1 - zc1
         zpreout(:,jj) = zc1 * zpre1(:,jj) + zc2 * zpre2(:,jj)
      ENDDO
      ! fill southern transition zone
      DO jj=zlat25s,zlat30s
         zc2 = FLOAT(( jj - zlat25s)) / FLOAT(( zlat30s - zlat25s ))
         zc1 = 1 - zc2
         zpreout(:,jj) = zc1 * zpre1(:,jj) + zc2 * zpre2(:,jj)
      ENDDO

      CALL check( NF90_PUT_VAR(ncid_preout, preout_varid, REAL(zpreout), start = start, &
                  &            count = count) )


   ENDDO

   CALL check( NF90_CLOSE(ncid_pre1) )
   CALL check( NF90_CLOSE(ncid_pre2) )
   CALL check( nf90_close(ncid_preout) )


CONTAINS

      SUBROUTINE check(status)
        INTEGER, INTENT(in) :: status
    
        IF(status /= nf90_noerr) THEN
          PRINT *, trim(nf90_strerror(status))
          STOP 2
        ENDIF
      END SUBROUTINE check


END PROGRAM
