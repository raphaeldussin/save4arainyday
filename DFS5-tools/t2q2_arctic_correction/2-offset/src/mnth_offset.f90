PROGRAM MNTH_OFFSET
  !!
  USE io_ezcdf
  USE drwn
  !!
  IMPLICIT none
  !!
  LOGICAL, PARAMETER :: lp_clndr = .FALSE.
  !!
  !!
  INTEGER, PARAMETER ::   &
       !!
       &   jj_max = 42, &  ! lat = 60.7
       !!
       !! RD : first try with t2 on 1989-2009 period
       !&   jy1   = 1989,  &
       !&   jy2   = 2009,  &
       !&   jytr1 = 2010,  &
       !&   jytr2 = 2010,  &
       &   jy1   = 1989,  &
       &   jy2   = 1998,  &
       &   jytr1 = 1999,  &
       &   jytr2 = 1999,  &
       &   n_smth = 4      ! number of iterations for smoothing fields
  !!
  REAL, PARAMETER :: &
       &  lat_min = -60.,  &
       &  lat_max =  65.,  &
       &  vmax    =  1.2      ! cut extreme values
  !!
  !!################## ORCA2 : #################################
  CHARACTER(len=256), PARAMETER :: &
  !& cgrid='/home/users/dussin/STORAGE/ERA-INTERIM/MASK/lsm_interim_512x256_lonlat.nc', &
  !& cmask='/home/users/dussin/STORAGE/ERA-INTERIM/MASK/lsm_interim_512x256_lonlat.nc', &
  !& cmask_orig='/home/users/dussin/STORAGE/ERA-INTERIM/MASK/lsm_interim_512x256_lonlat.nc', &
  & cgrid='/srv/share/data81/dussin/DFS5.1/tmp/lsm_interim_512x256_lonlat.nc', &
  & cmask='/srv/share/data81/dussin/DFS5.1/tmp/lsm_interim_512x256_lonlat.nc', &
  & cmask_orig='/srv/share/data81/dussin/DFS5.1/tmp/lsm_interim_512x256_lonlat.nc', &
  & cd_out = '.'
  !!
  INTEGER, PARAMETER ::   &
       ni     = 512,    &
       nj     = 256
  !!############################################################
  !!
  INTEGER :: &
       &  ifa, iva, ifb, ivb, ifc, ivc, ifd, ivd, ife, ive, iff, ivf,  &
       &  if1, iv1, if2, iv2, if3, iv3, if4, iv4, &
       &  ji, jj, jt, jl, jj1, jj2, jm,  &
       &  nby, nbmnth,nbsnp, nds,      &
       &  nbper, jper, jnbs, jday, jsd
  !!
  CHARACTER(len=20), PARAMETER :: &
       & cv_lon = 'lon',      &
       & cv_lat = 'lat', &
       & cv_t   = 'time'
  !!
  CHARACTER(len=256), PARAMETER ::   &
       !!
       !! File with with monthly POLES temperature on ERA-40 grid :
       !& cf_pt ='/home/users/dussin/WORK/PREPA_GRD100/forcages/CORRECTION_ARCTIC/t2_POLES-ERAi_1979-1998.nc',  &
       & cf_pt ='/home/users/dussin/TOOLS/t2q2_arctic_correction/0-data/t2_POLES-ERAi_1979-1998.nc',  &
       !!
       !! Mean ERA40 on period 0 to 2 :
       !& cf_p0 = '/home/users/dussin/WORK/PREPA_GRD100/forcages/CORRECTION_ARCTIC/t2/t2_ERAinterim-512x256_1989-1998_monthly.nc'
       & cf_p0 = '/srv/share/data81/dussin/DFS5.1/tmp/drowned_t2_ERAinterim_y1979-1988_monthly.nc'
      !! RD first try
      ! & cf_p0 = '/home/users/dussin/WORK/PREPA_GRD100/forcages/CORRECTION_ARCTIC/t2/t2_ERAinterim-512x256_1989-2009_monthly.nc'
  !!
  CHARACTER(len=200) :: clnm
  !!
  CHARACTER(len=100) ::   &
       & cf_out,   &
       & cf_misc,  &
       & cdirin,   &
       & cdirout
  !!
  REAL(4), PARAMETER :: &
       &   rmv     = -9999.0,  &
       &   rmv_rst = -200.0
  !!
  CHARACTER(len=20) ::   &
       & cv_in = 't2'
  !!
  CHARACTER(len=20) ::   &
       & cv_msc
  !!
  REAL(4), DIMENSION(:), ALLOCATABLE :: vyear, vper, vm_band, vam_band
  !!
  REAL(4), DIMENSION(:,:),   ALLOCATABLE :: R_zon
  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: R_zon_m
  !!
  REAL(4), DIMENSION(:,:,:),   ALLOCATABLE :: X_smt, XmD4
  !!
  REAL(8),    DIMENSION(ni)    :: vlon
  REAL(8),    DIMENSION(nj)    :: vlat  
  REAL(8),    DIMENSION(ni,nj) :: xlon
  REAL(8),    DIMENSION(ni,nj) :: xlat  

  REAL(8),    DIMENSION(12)    :: vtime
  !!
  REAL(4),    DIMENSION(ni,nj)    :: X2D
  REAL(4),    DIMENSION(ni,nj,12) :: Xpt, Xp0, Xp1, Xp2
  REAL(4),    DIMENSION(ni,nj)    :: Xvar_m, Xm_p0, Xm_p1, Xm_p2
  REAL(4),    DIMENSION(ni,nj,12) ::  Xoffst0, Xoffst1, Xoffst2, X0, X1, X2
  !!
  INTEGER, DIMENSION(ni,nj) :: &
       & mask_orig, mask, mask0, mask1
  !!
  REAL(4) :: val1, val2, rrm, rA, rB, slp0, slp1, slp2, dlat
  !!
  INTEGER :: &
       & narg, iarg,                &
       & jlat, jy, jycpt,       &
       & id_f1, id_v1, id_f2, id_v2,&
       & ip1, jp1, ip2, jp2
  !!
  !!
  CALL GETVAR_1D(cgrid, cv_lon, ni, vlon)
  CALL GETVAR_1D(cgrid, cv_lat, nj, vlat)

  DO jj=1,nj
     xlon(:,jj) = vlon(:)
  ENDDO

  DO ji=1,ni
     xlat(ji,:) = vlat(:)
  ENDDO


  cv_msc = 'lsm'; CALL GETMASK_2D(cmask_orig, cv_msc, ni, nj, mask_orig)
  cv_msc = 'lsm'; CALL GETMASK_2D(cmask, cv_msc, ni, nj, mask)
  !!
  !!
  X2D = 0.
  !!
  DO jm = 1, 12
     !!
     vtime(jm) = REAL(jm)
     !! First we need to create offsta maps:
     CALL GETVAR_2D(if1, iv1, cf_pt, cv_in, ni, nj, 12, 0, jm, Xpt(:,:,jm))
     CALL GETVAR_2D(if2, iv2, cf_p0, cv_in, ni, nj, 12, 0, jm, Xp0(:,:,jm))
     !CALL GETVAR_2D(if3, iv3, cf_p1, cv_in, ni, nj, 12, 0, jm, Xp1(:,:,jm))
     !CALL GETVAR_2D(if4, iv4, cf_p2, cv_in, ni, nj, 12, 0, jm, Xp2(:,:,jm))
     !!
  END DO
  !!
  !!
  Xoffst0 = Xpt - Xp0 
  !Xoffst1 = Xpt - Xp1 
  !Xoffst2 = Xpt - Xp2 
  Xoffst1 = Xpt - Xp0 
  Xoffst2 = Xpt - Xp0 

  !!
!  WHERE ( Xpt == -200. )
!RD
  WHERE ( Xpt < 200. )
     Xoffst0 = 0. ;  Xoffst1 = 0. ;  Xoffst2 = 0.
  END WHERE
  !!
  X0 = Xoffst0 ; X1 = Xoffst1 ; X2 = Xoffst2
  !!
  DO jm = 1, 12
     !CALL DROWN(ni, nj, Xoffst0(:,:,jm), mask, 30)
     !CALL DROWN(ni, nj, Xoffst1(:,:,jm), mask, 30)
     !CALL DROWN(ni, nj, Xoffst2(:,:,jm), mask, 30)
     CALL DROWN(0, Xoffst0(:,:,jm), mask)
     !CALL DROWN(0, Xoffst1(:,:,jm), mask)
     !CALL DROWN(0, Xoffst2(:,:,jm), mask)
  END DO
  !!
  !!
  DO jm = 1, 12
     !!
     !DO jj = jj_max, jj_max -13, -1
     DO jj = jj_max, jj_max + 21, 1  ! RD : regle de 3 pour ERAinterim
        Xoffst0(:,jj-1,jm) = Xoffst0(:,jj,jm)/1.6
        !Xoffst1(:,jj-1,jm) = Xoffst1(:,jj,jm)/1.6
        !Xoffst2(:,jj-1,jm) = Xoffst2(:,jj,jm)/1.6
     END DO
     !!
!RD test
     DO jl = 1, n_smth
        CALL SMOOTHER(ni,nj, mask, Xoffst0(:,:,jm))
     END DO
     !!
     WHERE ( mask_orig == 0 )
        Xoffst0(:,:,jm) = rmv ;  Xoffst1(:,:,jm) = rmv ;  Xoffst2(:,:,jm) = rmv
     END WHERE
     !!
  END DO
  !!
!111 CONTINUE
  !!
  !!
!  DO jm = 1, 12
     !!
!  END DO


!  CALL exp_msk(ni, nj, mask, mask0) ; mask1 = mask0
!  CALL exp_msk(ni, nj, mask1, mask0)
  !!
!  CALL DROWN(ni, nj, Xalf0, mask0, 30)
!  CALL DROWN(ni, nj, Xalf1, mask0, 30)
!  CALL DROWN(ni, nj, Xalf2, mask0, 30)
  !!
  !!
!  WHERE ( Xalf0 > vmax ) Xalf0 = vmax
!  WHERE ( Xalf1 > vmax ) Xalf1 = vmax
!  WHERE ( Xalf2 > vmax ) Xalf2 = vmax
  !!
  !!
  !! Blending from value to 1 :
!  CALL BLENDER(ni, nj, lat_min, lat_max, vlat, Xalf0)
!  CALL BLENDER(ni, nj, lat_min, lat_max, vlat, Xalf1)
  ! CALL BLENDER_PERSIST(ni, nj, lat_min, lat_max, vlat, Xalf1)
!  CALL BLENDER(ni, nj, lat_min, lat_max, vlat, Xalf2)
  !!
  !!
  !! Smoothing
!  DO jl = 1, n_smth
!     CALL SMOOTHER(ni,nj, mask, Xalf0)
!     CALL SMOOTHER(ni,nj, mask, Xalf1)
!     CALL SMOOTHER(ni,nj, mask, Xalf2)
!  END DO
  !!
  !!
!  CALL DROWN(ni, nj, Xalf0, mask, 30)
!  CALL DROWN(ni, nj, Xalf1, mask, 30)
!  CALL DROWN(ni, nj, Xalf2, mask, 30)
  !!
  !!
!  WHERE (mask_orig == 0)
!     Xoffst0 = rmv ;      Xoffst1 = rmv ;      Xoffst2 = rmv
!  END WHERE
  !!
  cv_msc = 'offst'
  !!
  DO jt = 1, 12
  !!
     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
          &   jy1, jytr1-1, vmax
     WRITE(cf_out, '(a,"/offset_poles-eraint_",i4,"-",i4,".nc")') trim(cd_out), jy1, jytr1-1
!     CALL P2D_T(ifa, iva, ni, nj, 12, jt, vlon, vlat, vtime, Xoffst0(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)

     CALL P2D_T(ifa, iva, 12, jt, REAL(xlon), REAL(xlat), vtime, Xoffst0(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
          &       cv_msc, cv_msc, clnm, rmv)
     !!
!     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
!          &   jytr1, jytr2-1, vmax
!     WRITE(cf_out, '(a,"/offset_poles-era40_",i4,"-",i4,".nc")') trim(cd_out), jytr1, jytr2-1
!     CALL P2D_T(ifb, ivb, ni, nj, 12, jt, vlon, vlat, vtime, Xoffst1(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)
!     !!
!     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
!          &   jytr2, jy2, vmax
!     WRITE(cf_out, '(a,"/offset_poles-era40_",i4,"-",i4,".nc")') trim(cd_out), jytr2, jy2
!     CALL P2D_T(ifc, ivc, ni, nj, 12, jt, vlon, vlat, vtime, Xoffst2(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)
     !!
     !!
     !!
!     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
!          &   jy1, jytr1-1, vmax
!     WRITE(cf_out, '(a,"/B_offset_poles-era40_",i4,"-",i4,".nc")') trim(cd_out), jy1, jytr1-1
!     CALL P2D_T(ifd, ivd, ni, nj, 12, jt, vlon, vlat, vtime, X0(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)
     !!
!     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
!          &   jytr1, jytr2-1, vmax
!     WRITE(cf_out, '(a,"/B_offset_poles-era40_",i4,"-",i4,".nc")') trim(cd_out), jytr1, jytr2-1
!     CALL P2D_T(ife, ive, ni, nj, 12, jt, vlon, vlat, vtime, X1(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)
     !!
!     WRITE(clnm,'("Correction coefficient for ERA40 winds, ",i4,"-",i4,", max = ", f4.2)') &
!          &   jytr2, jy2, vmax
!     WRITE(cf_out, '(a,"/B_offset_poles-era40_",i4,"-",i4,".nc")') trim(cd_out), jytr2, jy2
!     CALL P2D_T(iff, ivf, ni, nj, 12, jt, vlon, vlat, vtime, X2(:,:,jt), cf_out, cv_lon, cv_lat, cv_t, &
!          &       cv_msc, cv_msc, clnm, rmv)
     !!
  END DO
  !!
  !!
END PROGRAM MNTH_OFFSET
!!
!!
SUBROUTINE SMOOTHER(lx,ly, msk, X)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  INTEGER, DIMENSION(lx,ly), INTENT(in) :: msk
  REAL(4), DIMENSION(lx,ly), INTENT(inout) :: X
  INTEGER i, j
  REAL(4) :: w1, w2, w3, w4, sumw
  !!
  DO i=2, lx-1
     DO j=2, ly-1
        !!
        w1 = 1. ;   w2 = 1. ;  w3 = 1. ;  w4 = 1.
        !!
        IF ( msk(i,j) == 1 ) THEN
           !!
           IF ( msk(i+1,j) == 0 ) w1 = 0.
           IF ( msk(i,j+1) == 0 ) w2 = 0.
           IF ( msk(i-1,j) == 0 ) w3 = 0.
           IF ( msk(i,j-1) == 0 ) w4 = 0.
           !!
           sumw = w1 + w2 + w3 + w4
           !!
           IF ( sumw /= 0. ) THEN
              X(i,j) = 0.5*X(i,j)  &
                   & + 0.5/sumw*(w1*X(i+1,j) + w2*X(i,j+1) + w3*X(i-1,j) + w4*X(i,j-1))
           END IF
           !!
        END IF
     END DO
  END DO
  !!
END SUBROUTINE SMOOTHER
!!
!!
SUBROUTINE BLENDER(lx, ly, l_min, l_max, vl, X)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  REAL(4), INTENT(in) :: l_min, l_max
  !!
  REAL(4), DIMENSION(ly),    INTENT(in) :: vl
  REAL(4), DIMENSION(lx,ly), INTENT(inout) :: X
  !!
  INTEGER i, j, j1, j2
  REAL(4) :: slp
  INTEGER, PARAMETER :: lbld = 10
  !!
  !!
  DO j=1, ly - 1
     IF ( (vl(j)    < l_min).and.(vl(j+1) >= l_min) ) j1 = j+1
     IF ( (vl(j+1) >= l_max).and.(vl(j)   <  l_max) ) j2 = j
  END DO
  !!
  !!
  DO i = 1, lx
     slp = (X(i,j1-lbld) - X(i,j1))/(vl(j1-lbld) - vl(j1))
     DO j=j1-1, j1-lbld, -1
        X(i,j) = X(i,j1) + slp*(vl(j) - vl(j1))
     END DO
     !!
     slp = (X(i,j2+lbld) - X(i,j2))/(vl(j2+lbld) - vl(j2))
     DO j=j2+1, j2+lbld
        X(i,j) = X(i,j2) + slp*(vl(j) - vl(j2))
     END DO
  END DO
  !!
END SUBROUTINE BLENDER
!!
!!
SUBROUTINE BLENDER_PERSIST(lx, ly, l_min, l_max, vl, X)
  !!
  IMPLICIT none
  !!
  INTEGER, INTENT(in) :: lx, ly
  REAL(4), INTENT(in) :: l_min, l_max
  !!
  REAL(4), DIMENSION(ly),    INTENT(in) :: vl
  REAL(4), DIMENSION(lx,ly), INTENT(inout) :: X
  !!
  INTEGER i, j, j1, j2
  REAL(4) :: slp
  INTEGER, PARAMETER :: lbld = 10
  !!
  !!
  DO j=1, ly - 1
     IF ( (vl(j)    < l_min).and.(vl(j+1) >= l_min) ) j1 = j+1
     IF ( (vl(j+1) >= l_max).and.(vl(j)   <  l_max) ) j2 = j
  END DO
  !!
  !!
  DO i = 1, lx
     DO j=j1-1, j1-lbld, -1
        X(i,j) = X(i,j1)
     END DO
     !!
     DO j=j2+1, j2+lbld
        X(i,j) = X(i,j2)
     END DO
  END DO
  !!
END SUBROUTINE BLENDER_PERSIST
!!
!!
!!
SUBROUTINE exp_msk(lx, ly, xmsk_in, xmsk_out)
  !!
  !!
  INTEGER, INTENT(in)  :: lx, ly
  !!
  INTEGER, DIMENSION(lx,ly), INTENT(in)   :: xmsk_in
  INTEGER, DIMENSION(lx,ly), INTENT(out)  :: xmsk_out
  !!
  !!
  INTEGER :: ji, jj
  !!
  xmsk_out = xmsk_in
  !!
  DO ji=2, lx-1
     DO jj=2, ly-1
        !!
        !!
        !! Western coasts :
        IF ( (xmsk_in(ji+1,jj) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji-1,jj) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !! Eastern coasts :
        IF ( (xmsk_in(ji-1,jj) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji+1,jj) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !!
        !! Southern coasts :
        IF ( (xmsk_in(ji,jj+1) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji,jj-1) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
        !! Northern coasts :
        IF ( (xmsk_in(ji,jj-1) == 0)      &
             &  .and.(xmsk_in(ji,jj) /= 0) ) THEN !.and.(xmsk_in(ji,jj+1) /= 0) ) THEN
           xmsk_out(ji,jj) = 0
        END IF
        !!
     END DO
  END DO
  !!
  !!
END SUBROUTINE exp_msk
!!
!!
