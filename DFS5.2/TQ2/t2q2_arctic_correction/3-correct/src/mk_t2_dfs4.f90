
PROGRAM MK_T2_DFS4
  !!
  USE io_ezcdf
  !!
  IMPLICIT none
  !!
  !! Leap years or not
  LOGICAL, PARAMETER :: lp_clndr = .TRUE.
  !!
  !!
  INTEGER, PARAMETER ::   & 
       &   jy1   = 1979,  &
       &   jyf   = 2012,  &
       &   j70   = 29,    &   !! jj = 70 degN
       &   ndtr  = 13

  REAL, PARAMETER :: val_oce = -0.7
  !!
  !!
  CHARACTER(len=512), PARAMETER :: &
       &  cgrid = '/fsnet/data/meom/workdir/dussin/MESHMASKS/lsm_erainterim.nc', &
       &  cmask_orig = '/fsnet/data/meom/workdir/dussin/MESHMASKS/lsm_erainterim.nc', &       
       &  cd_str  = '/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ANTARCTIC',  &
       !&  cf_ice_mm ='/home/users/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/0-data/ifrac_SSMI-ERAi_1979-1998.nc', &
       &  cf_ice_mm ='/fsnet/data/meom/home/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/0-data/ifrac_SSMI-ERAi_1979-1998.nc', &
       &  cd_out  = '/fsnet/data/meom/workdir/dussin/TMPDIR_TQ2_CORR_ARCTIC'
  !!
  INTEGER, PARAMETER ::   &
       ni     = 512,      &
       nj     = 256

  INTEGER :: nspy
  INTEGER :: nday
  INTEGER :: freq = 8
  LOGICAL :: isleap

  !!############################################################
  !!
  INTEGER :: ji, jj, jt, js, jjs, jjn, jm, jdm, jd, jsd, jcs
  !!
  INTEGER, DIMENSION(12) ::       &
       &   nb_m = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  !!
  CHARACTER(len=20), PARAMETER :: &
       & cv_lon = 'lon',   &
       & cv_lat = 'lat',   &
       & cv_t   = 'time',  &
       & cv_im  = 'ifrac'
!       & cv_im  = 'ice_mask'
  !!
  CHARACTER(len=200) :: &
       & clnm = 'ERA-interim arctic correction (Brodeau et al. 2007), temperature at 2m'
  !!
  CHARACTER(len=256) ::   &
       & cf_in,   &
       & cf_out,  &
       & cf_corr_out,  &
       & cf_m, cf_m1, cf_m2
  !!
  REAL(4), PARAMETER :: &
       &   rmv = -9999.0
  !!
  CHARACTER(len=20) ::   &
       & cv_in = 't2', &
       & cu    = 'K'
  !!
  CHARACTER(len=20) ::   &
       & cv_m
  !!
  REAL(8),    DIMENSION(:), ALLOCATABLE      :: vt
  REAL(8),    DIMENSION(:), ALLOCATABLE      :: vtime
  REAL(4),    DIMENSION(:,:,:), ALLOCATABLE  :: Xc_xmd, Xc_xid, Xc_corr

  REAL(8),    DIMENSION(ni)       :: vlon
  REAL(8),    DIMENSION(nj)       :: vlat, valfa
  REAL(8),    DIMENSION(ni,nj)    :: xlon, xlat

  REAL(4),    DIMENSION(ni,nj)    :: Xr4, X2D, Xmd, Xm, Xm_b, Xm_n, Xcorr
  REAL(4),    DIMENSION(ni,nj)    :: Xid, Xi, Xi_b, Xi_n
  REAL(4),    DIMENSION(ni,nj,12) :: Xoff, ice_mask
  !!
  INTEGER, DIMENSION(ni,nj) :: mask
  !!
  INTEGER :: &
       & jlat, jy, jycpt,              &
       & ifx, ivx, &
       & ifa, iva, ifb, ivb, ifc, ivc, &
       & id_f1, id_v1, id_f2, id_v2,   &
       & id_f2c, id_v2c,  &
       & jmc
  !!
  !!
  !! Getting monthly climatic sea-ice mask from SSMI :
  DO jm = 1, 12
     X2D = 0.
     CALL GETVAR_2D(ifx, ivx, TRIM(cf_ice_mm), cv_im, ni, nj, 12, 0, jm, ice_mask(:,:,jm))
  END DO
  !!
  ifx= 0 ;  ivx = 0
  jmc = 0        ! current month along the ny years...
  !!
  !!
  CALL GETVAR_1D(cgrid, cv_lon, ni, vlon)
  CALL GETVAR_1D(cgrid, cv_lat, nj, vlat)
  !!
  cv_m = 'lsm' ; CALL GETMASK_2D(cmask_orig, cv_m, ni, nj, mask)
  !!
  DO jj=1,nj
     xlon(:,jj) = vlon(:)
  ENDDO

  DO ji=1,ni
     xlat(ji,:) = vlat(:)
  ENDDO


  !!
  !!
  cv_m = 'offst'
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! LOOP ON YEARS !!
  DO jy = jy1, jyf
     !!
     IF ( MOD(jy - 1976,4) == 0 ) THEN
        isleap = .TRUE. ; nspy = freq * 366 ; nday = 366
        nb_m = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
     ELSE
        isleap = .FALSE. ; nspy = freq * 365 ; nday = 365
     ENDIF

     ALLOCATE( vt(nday) , vtime(nspy) )
     ALLOCATE( Xc_xmd(ni,nj,nday) , Xc_xid(ni,nj,nday) , Xc_corr(ni,nj,nday) )

     IF ( jy == jy1 ) THEN
        !cf_m = '/data2/dussin/WORK/PREPA_GRD100/forcages/CORRECTION_ARCTIC/offset_poles/offset_poles-eraint_1989-1998.nc'
        cf_m = '/fsnet/data/meom/home/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/2-offset/offset_poles-eraint_1979-1998.nc'
        DO jm = 1, 12
           CALL GETVAR_2D(ifx, ivx, cf_m, cv_m, ni, nj, 12, 0, jm, Xoff(:,:,jm))
        END DO
        ifx = 0 ; ivx = 0
     END IF
     !!
     !IF ( jy == jytr1 ) THEN
     !   cf_m = '/stock0/stock/brodeau/POLES/OFFSETS/offset_poles-era40_1979-1999.nc'
     !   DO jm = 1, 12
     !      CALL GETVAR_2D(ifx, ivx, cf_m, cv_m, ni, nj, 12, 0, jm, Xoff(:,:,jm))
     !   END DO
     !   ifx = 0 ; ivx = 0
     !END IF
     !!
     !IF ( jy == jytr2 ) THEN
     !   cf_m = '/stock0/stock/brodeau/POLES/OFFSETS/offset_poles-era40_2000-2007.nc'
     !   DO jm = 1, 12
     !      CALL GETVAR_2D(ifx, ivx, cf_m, cv_m, ni, nj, 12, 0, jm, Xoff(:,:,jm))
     !   END DO
     !   ifx = 0 ; ivx = 0
     !END IF
     !!
     !!
     PRINT *, 'Year =', jy
     !!
     !! On ne corrige que de 1958 a 1978 a l'equateur
     !IF ( jy == jy1   ) CALL MK_ALFA_VECT(nj, vlat, jjs, jjn, wdth_int, alfa1, valfa)
     !IF ( jy >= jytr1 ) valfa(:) = 1.
     valfa(:) = 1.
     !CALL MK_ALFA_VECT(nj, vlat, jjs, jjn, wdth_int, alfa2, valfa)
     !!
     !WRITE(cf_in, '(a,"/",a,"_INTERIM-512x256_y",i4,".nc")') &
     WRITE(cf_in, '(a,"/drowned_",a,"_ERAinterim_corr_antarctic_y",i4,".nc")') &
          &         trim(cd_str), trim(cv_in), jy
     PRINT *, 'Treated file = ', trim(cf_in)
     !!
     !WRITE(cf_out, '(a,"/",a,"_DFS5.0_y",i4,".nc")') &
     WRITE(cf_out, '(a,"/",a,"_DFS5.2_y",i4,".nc")') &
          &        trim(cd_out), trim(cv_in), jy
     PRINT *, 'Created file = ', trim(cf_out)
     !!
     !WRITE(cf_corr_out, '("corr_",a,"_DFS5.0_y",i4,".nc")') trim(cv_in), jy
     WRITE(cf_corr_out, '(a,"/","corr_",a,"_DFS5.2_y",i4,".nc")') trim(cd_out), trim(cv_in), jy
     !!
     DO jd = 1, nday
        vt(jd) = REAL(jd)
     END DO
     !!
     DO js = 1, nspy
        vtime(js) = (1. / freq) *js
     END DO
     !!
     !!
     jm  = 1 ; jmc = jmc + 1
     jdm = 1
     jd  = 1
     jcs = 0
     !!
     jsd = 1
     !!
     DO js = 1, nspy
        !!
        IF ( jdm == nb_m(jm)+1 ) THEN
           jm  = jm  + 1 ; jmc = jmc + 1 ; jdm = 1
        END IF
        !!
        !        PRINT *, 'year, month, day, day of month, snap:'
!        PRINT *, jy,    jm,    jd,  jdm,          js 
        !!
        !!
        !!
        IF ( jsd == 1 ) THEN
           !! ===========================================================
           !! What temperature offset for current day?
           !! ===========================================================
           !! Transition on the first and last ndtr days of the current month
           !!
           !! Valeur du moi courant :
           Xm = Xoff(:,:,jm)
           Xi = ice_mask(:,:,jm)
           !!
           !! Valeur a appliquer au jour courant :
           Xmd = Xoff(:,:,jm) ! pour le reste des cas...
           Xid = ice_mask(:,:,jm)
           !!
           IF ( jdm <= ndtr ) THEN
              !!
              !! Valeur du mois d'avant :
              IF ( jm == 1 ) THEN
                 Xm_b = Xoff(:,:,12)
                 Xi_b = ice_mask(:,:,12)
              ELSE
                 Xm_b = Xoff(:,:,jm-1)
                 Xi_b = ice_mask(:,:,jm-1)
              END IF
              !!
              Xmd = (Xm_b + Xm)/2 + jdm*(Xm - (Xm_b+Xm)/2)/ndtr
              Xid = (Xi_b + Xi)/2 + jdm*(Xi - (Xi_b+Xi)/2)/ndtr
              !!
           END IF
           !!
           IF ( jdm >= nb_m(jm) - ndtr + 1 ) THEN
              !!
              !! Valeur du prochain mois :
              IF ( jm == 12 ) THEN
                 Xm_n = Xoff(:,:,1)
                 Xi_n = ice_mask(:,:,1)
              ELSE
                 Xm_n = Xoff(:,:,jm+1)
                 Xi_n = ice_mask(:,:,jm+1)
              END IF
              !!
              Xmd = (Xm_n + Xm)/2 + (nb_m(jm) - jdm + 1)*(Xm - (Xm_n+Xm)/2)/ndtr
              Xid = (Xi_n + Xi)/2 + (nb_m(jm) - jdm + 1)*(Xi - (Xi_n+Xi)/2)/ndtr
              !!
           END IF
           !!
           !! Keeping tracks :
           Xc_xmd(:,:,jd) = Xmd
           Xc_xid(:,:,jd) = Xid
           !!
           !! ===========================================================
           !!
        END IF
        !!
        !!
        CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, ni, nj, nspy, 0, js, Xr4)
        !!
        !! Correction :
        !DO jj = 1, nj
        !   Xr4(:,jj) = valfa(jj)*Xr4(:,jj)
        !END DO
        !!
        !! Northern daily correction for current snapshot with Poles
        !!
        !! 1) It should only be north of 70N:
        Xcorr = 0.
        !DO jj = j70, nj  !! a corriger
        DO jj = 1, j70
           Xcorr(:,jj) = Xmd(:,jj)*Xid(:,jj) + val_oce*(1. - Xid(:,jj))
        END DO
        !!
        DO ji = 1, ni
           !DO jj = j70-1, j70-5, -1 !! a corriger
           DO jj = j70+1, j70+7  ! ou 8 ?
              IF (mask(ji,jj+1) == 1) Xcorr(ji,jj) = Xcorr(ji,jj+1)/1.3
           END DO
        END DO
        !!
        !DO ji = 1, 4
        DO ji = 1, 8 !! RD increase smooting
           CALL SMOOTHER(ni,nj,mask,Xcorr)
        END DO
        !!
        !!
        !!
        Xr4 = Xr4 + Xcorr
        !!
        !!
        WHERE ( mask == 0 )
           Xr4 = rmv
           Xcorr = rmv
        END WHERE
        !!
        IF ( jsd == 1 ) Xc_corr(:,:,jd) = Xcorr
        !!
        CALL P2D_T(id_f2, id_v2, nspy, js, REAL(xlon), REAL(xlat), vtime, &
             &     Xr4, cf_out, cv_lon, cv_lat, cv_t, cv_in, cu, clnm, rmv)
        !!
        jsd = jsd + 1        
        IF ( jsd == freq + 1 ) THEN
           jsd = 1
           jd  = jd  + 1
           jdm = jdm + 1
        END IF
        !!

        !!
     END DO
     !!
     !! Annee courante terminee :
     !!
     !!SUBROUTINE P3D_T(id_fil, id_var, lx, ly, lz, lt, lct, vlon, vlat, vdpth, vtime, &
     !!  &           x3d, cfil, cvarlon, cvarlat, cvardpth, cvartime, cvar, cunit, &
     !!  &           cln, vflag)
     !!
     !SUBROUTINE P2D_T(id_fil, id_var, lx, ly, lt, lct, vlon, vlat, vtime, &
     !       &         x2d, cfil, cvarlon, cvarlat, cvartime, cvar, cunit, cln, vflag)
     !!
     WRITE(cf_m1, '(a,"/","xmd_offset_",i4,".nc")'), trim(cd_out),jy
     WRITE(cf_m2, '(a,"/","xid_offset_",i4,".nc")'), trim(cd_out),jy
     !!
     DO jt = 1, nday
        !!
        CALL P2D_T(ifa, iva, nday, jt, REAL(xlon), REAL(xlat), vt, &
             &   Xc_corr(:,:,jt), cf_corr_out, cv_lon, cv_lat, cv_t, cv_m, cu, clnm, rmv)
        !!

        CALL P2D_T(ifb, ivb, nday, jt, REAL(xlon), REAL(xlat), vt, &
             &   Xc_xmd(:,:,jt), cf_m1, cv_lon, cv_lat, cv_t, cv_m, cu, clnm, rmv)
        !!
        CALL P2D_T(ifc, ivc, nday, jt, REAL(xlon), REAL(xlat), vt, &
             &    Xc_xid(:,:,jt), cf_m2, cv_lon, cv_lat, cv_t, cv_m, cu, clnm, rmv)
        !!
     END DO
     !!
     PRINT *, ''; PRINT *, ''; PRINT *, ''
     !!
     DEALLOCATE( vt, vtime )
     DEALLOCATE( Xc_xmd , Xc_xid , Xc_corr )
     !!
  END DO
  !!
  !!
END PROGRAM MK_T2_DFS4
!!
!!
SUBROUTINE MK_ALFA_VECT(ly, vy, js, jn, jdst, alf, va)
  !!
  IMPLICIT none
  !!
  INTEGER,                INTENT(in)  :: ly
  REAL(4), DIMENSION(ly), INTENT(in)  :: vy
  INTEGER,                INTENT(in)  :: js, jn, jdst
  REAL(4),                INTENT(in)  :: alf
  REAL(4), DIMENSION(ly), INTENT(out) :: va
  !!
  REAL(4) :: aa
  !!
  INTEGER :: j
  !!
  !!
  !! Creating alfa vector (interpolation on jdst grid points):
  va(:) = 1.
  DO j = js - jdst + 1, js
     aa = (vy(j) - vy(js - jdst + 1))/( vy(js) - vy(js - jdst + 1) )
     va(j)  = 1.*(1. - aa) + aa*alf
  END DO
  DO j = js + 1, jn - 1
     va(j) = alf
  END DO
  DO j = jn, jn + jdst - 1
     aa = (vy(j) - vy(jn + jdst - 1))/( vy(jn) - vy(jn + jdst - 1) )
     va(j) = 1.*(1. - aa) + aa*alf
  END DO
  !!
  OPEN(unit=11, file = 'vect_alf.out', status='unknown')
  DO j=1, ly
     WRITE(11,'(f,"  ",f)') vy(j), va(j)
  ENDDO
  CLOSE(11)
  !!
  !!
END SUBROUTINE MK_ALFA_VECT
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
