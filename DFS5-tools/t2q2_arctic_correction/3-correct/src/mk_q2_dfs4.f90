PROGRAM MK_Q2_DFS4
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
       &   jy1   = 1989,  &
       &   jy2   = 2009,  &
       &   jstart = 56    ! a la louche
      ! &   jstart = 128    ! a corriger a partir d'ou aplliquer la correction due aux t2 polaires corrigees
  !!
  !!
  CHARACTER(len=256), PARAMETER :: &
       &  cgrid   ='/home/users/dussin/STORAGE/ERA-INTERIM/MASK/lsm_interim_512x256_lonlat.nc', &
       &  cf_msk ='/home/users/dussin/STORAGE/ERA-INTERIM/MASK/lsm_interim_512x256_lonlat.nc', &       
       &  cd_str  = '/data2/dussin/STORAGE/FORCINGS/ERAinterim/foto_ready',  &
       &  cd_dfs4  = '.'
  !!
  INTEGER, PARAMETER ::   &
       ni     = 512,    &
       nj     = 256

  INTEGER :: nspy , freq = 8
  !!############################################################
  !!
  INTEGER :: ji, jj, js, jjs, jjn
  
  !!
  CHARACTER(len=20), PARAMETER :: &
       & cv_lon = 'lon',   &
       & cv_lat = 'lat',   &
       & cv_tm  = 'time'
  !!
  CHARACTER(len=200) :: &
       & clnm = 'ERA-interim corrected for DFS5 (Brodeau et al. 2007), specific humidity at 2m'
  
  !!
  CHARACTER(len=256) ::   &
       & cf_in,  &
       & cf_tn, cf_to,   & ! new and old temperature 
       & cf_p,   &
       & cf_out, &
       & cf_misc
  !!
  REAL(4), PARAMETER :: &
       &   rmv = -9999.0
  !!
  CHARACTER(len=20) ::   &
       & cv_in = 'q2', &
       & cv_p  = 'MSL', &
       & cv_t  = 't2', &
       & cu    = 'kg/kg'
  !!
  CHARACTER(len=20) ::   &
       & cv_msc
  !!
  REAL(8),    DIMENSION(ni)       :: vlon
  REAL(8),    DIMENSION(nj)       :: vlat
  REAL(8),    DIMENSION(ni,nj)    :: xlon, xlat
  REAL(8),    DIMENSION(:), ALLOCATABLE     :: vtime
  REAL(4),    DIMENSION(ni,nj)    :: Xr4, Xr4p, Xr4tn, Xr4to, Xr4q
  !!
  INTEGER, DIMENSION(ni,nj) :: mask
  !!
  INTEGER :: &
       & jlat, jy, jycpt,       &
       & id_f1, id_v1, &
       & id_f2, id_v2, &
       & id_f3, id_v3, &
       & id_f4, id_v4, &
       & id_f5, id_v5
  !!
  REAL(4) :: q_sat, qs_n , qs_o
  !!
  PRINT *, ''
  PRINT *, 'Attention, on a besoin des temperatures du DFS4 dabord!!!!'; PRINT *, ''
  !!
  !!
  CALL GETVAR_1D(cgrid, cv_lon, ni, vlon)
  CALL GETVAR_1D(cgrid, cv_lat, nj, vlat)
  !!
  DO jj=1,nj
     xlon(:,jj) = vlon(:)
  ENDDO

  DO ji=1,ni
     xlat(ji,:) = vlat(:)
  ENDDO

  cv_msc = 'lsm'; CALL GETMASK_2D(cf_msk, cv_msc, ni, nj, mask)
  !!
  !!
  DO jy = jy1, jy2
     !!

     IF ( MOD(jy - 1988,4) == 0 ) THEN
        nspy = freq * 366 
     ELSE
        nspy = freq * 365 
     ENDIF

     ALLOCATE( vtime(nspy) )

     PRINT *, 'Year =', jy
     !!
     !WRITE(cf_in, '(a,"/",a,"/",a,".",i4,".era40.nc")') &
     WRITE(cf_in, '(a,"/",a,"_INTERIM-512x256_y",i4,".nc")') &
          &         trim(cd_str), trim(cv_in), jy
     PRINT *, 'Treated file = ', trim(cf_in)
     !!
     !WRITE(cf_p, '(a,"/",a,"/",a,".",i4,".era40.nc")') &
     WRITE(cf_p, '(a,"/","msl_ERAinterim_y",i4,".nc")') &
          &         trim(cd_str), jy
     PRINT *, 'Treated pressure file = ', trim(cf_p)
     !!
     !WRITE(cf_tn, '(a,"/t2/",a,"_DFS4_",i4,".nc")') &
     WRITE(cf_tn, '(a,"/",a,"_DFS5.0_y",i4,".nc")') &
          &        trim(cd_dfs4), trim(cv_t), jy
     PRINT *, 'Treated file new temperature = ', trim(cf_tn)
     !!
     !WRITE(cf_to, '(a,"/",a,"/",a,".",i4,".era40.nc")') &
     WRITE(cf_to, '(a,"/",a,"_INTERIM-512x256_y",i4,".nc")') &
          &         trim(cd_str), trim(cv_t), jy
     PRINT *, 'Treated file old temperature = ', trim(cf_to)
     !!
     !!
     WRITE(cf_out, '(a,"/",a,"_DFS5.0_",i4,".nc")') &
          &        trim(cd_dfs4), trim(cv_in), jy
     PRINT *, 'Created file = ', trim(cf_out)
     !!
     DO js = 1, nspy
        vtime(js) = 0.25*js
     END DO
     !!
     !!
     DO js = 1, nspy
        !!
        CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, ni, nj, nspy, 0, js, Xr4q)
        CALL GETVAR_2D(id_f2, id_v2, cf_p, cv_p, ni, nj, nspy, 0, js, Xr4p)
        !!
        CALL GETVAR_2D(id_f3, id_v3, cf_tn, cv_t, ni, nj, nspy, 0, js, Xr4tn) ! new temp
        CALL GETVAR_2D(id_f4, id_v4, cf_to, cv_t, ni, nj, nspy, 0, js, Xr4to) ! old temp
        !!
        Xr4 = Xr4q
        !!
        !! Correction Nordique, due a la grosse modification de temperature :
        DO ji = 1, ni
           !DO jj = jstart, nj
           DO jj = 1, jstart
              !!
              !! Calcule du nouveau q_sat:
              qs_n = q_sat(Xr4tn(ji,jj), Xr4p(ji,jj))*1000.
              !! Calcule du vieux q_sat:
              qs_o = q_sat(Xr4to(ji,jj), Xr4p(ji,jj))*1000.
              !!
!              IF ( mask(ji,jj) /= 0 ) THEN
              !                 PRINT *, ''
              !                 PRINT *, 'P =', Xr4p(ji,jj)
              !                 PRINT *, 't_o, t_n = ', Xr4to(ji,jj), Xr4tn(ji,jj)
              !                 PRINT *, 'qs_o, qs_n', qs_o, qs_n
              !                 PRINT *, ''
              !              END IF
              !!
              Xr4(ji,jj) = Xr4q(ji,jj)*qs_n/qs_o
              !!
           END DO
        END DO
        !!
        WHERE ( mask == 0 )
           Xr4 = rmv
        END WHERE
        !!
        CALL P2D_T(id_f5, id_v5, nspy, js, REAL(xlon), REAL(xlat), vtime, &
             &     Xr4, cf_out, cv_lon, cv_lat, cv_tm, cv_in, cu, clnm, rmv)
        !!
     END DO
     !!
     
     DEALLOCATE( vtime )
     PRINT *, ''; PRINT *, ''; PRINT *, ''
     !!
  END DO
  !!
  !!
END PROGRAM MK_Q2_DFS4
!!
!!
!!
FUNCTION q_sat(rt, rp)
  !!
  IMPLICIT none
  !!
  REAL(4)             :: q_sat   !: vapour pressure at saturation  [Pa]
  REAL(4), INTENT(in) :: rt, rp  !: temperature (K), pression (Pa)
  !!
  !!
  REAL(4)  :: es
  REAL(4), PARAMETER :: eps = 0.62197
  !!
  es = 100*( 10**(10.79574*(1 - 273.16/rt) - 5.028*LOG10(rt/273.16)           &
       &       + 1.50475*10**(-4)*(1 - 10**(-8.2969*(rt/273.16 - 1)) )           &
       &       + 0.42873*10**(-3)*(10**(4.76955*(1 - 273.16/rt)) - 1) + 0.78614) )
  !!
  q_sat = eps*es/(rp - (1 - eps)*es)
  !!
END FUNCTION q_sat



