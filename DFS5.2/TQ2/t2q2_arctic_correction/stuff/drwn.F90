MODULE DRWN
  !!
  !!
  IMPLICIT none
  !!
  INTERFACE IDRWN
     module procedure DROWN
     module procedure TEST_XYZ
  END INTERFACE
  !!
  PRIVATE
  !!
  PUBLIC :: drown, fill_extra_bands, extra_2_east, extra_2_west, test_xyz
  !!
  LOGICAL, PARAMETER :: lsmooth = .TRUE.
  !!
CONTAINS
  !!
  !!
  SUBROUTINE DROWN(k_ew, X, mask)
    !!
    !!#############################################################################
    !!
    !!  PURPOSE : fill land values only with surrounding sea values
    !!  --------
    !!  X   : currently treated array                    (2D array)
    !!  mask : land-sea mask                              (2D array)
    !!  n_inc : band width (in mesh) to propagate sea values into a continent
    !!
    !!#############################################################################
    !!
    !!
    !! Arguments :
    !! -----------
    INTEGER,                 INTENT(in)    :: k_ew
    REAL(4), DIMENSION(:,:), INTENT(inout) :: X
    INTEGER, DIMENSION(:,:),    INTENT(in) :: mask
    !!
    !!
    !! Local :
    !! --------
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: maskv, mask2 
    !!
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: data_old  
    !!
    INTEGER :: &
         &      ni, nj,        &
         &      jinc,          &
         &      ji, jj, jt,    &
         &      ncpt, nsea,    &
         &      nd_north,      &
         &      nd_east,       &
         &      nd_south,      &
         &      nd_west, kk
    !!  
    CHARACTER(LEN=8) :: quer
    !!
    INTEGER, PARAMETER ::   &
         &      nband = 2,  &
         &      nprob = 1     !: iterative distance to probe for sea presenc
    !!
    !!
    IF ( (size(X,1) /= size(mask,1)).OR.(size(X,2) /= size(mask,2)) ) THEN
       PRINT *, 'ERROR, drwn.F90 => DROWN : size of data and mask do not match!!!'; STOP
    END IF
    !!
    ni = size(X,1) 
    nj = size(X,2)
    !!
    IF ( SUM(mask) == 0 ) THEN
       PRINT *, 'The mask does not have sea points! Skipping DROWN!'
       RETURN
    ELSE
    !!
    !! Backing up original mask into mask2(:,:)
    ALLOCATE ( maskv(ni,nj), mask2(ni,nj), data_old(ni,nj) )
    maskv = mask
    mask2 = mask
    !!
    !!
    jinc = 0
    DO WHILE ( ANY(maskv(2:ni-2,2:nj-2) == 0) ) ! as long as we still have continental points
       jinc = jinc + 1
       !!
       maskv    = mask2
       data_old = X
       !!
       DO jj=nband, nj-nband+1
          ji = 1
          DO WHILE ( ji <= ni )
             IF( maskv(ji,jj) == 0 ) THEN
                CALL TEST4SEA(ji,jj,ni,nj,data_old,maskv,nprob,1,nsea)
                nd_north = nsea
                !!              
                CALL TEST4SEA(ji,jj,ni,nj,data_old,maskv,nprob,2,nsea)
                nd_east = nsea
                !!              
                CALL TEST4SEA(ji,jj,ni,nj,data_old,maskv,nprob,3,nsea)
                nd_south = nsea
                !!              
                CALL TEST4SEA(ji,jj,ni,nj,data_old,maskv,nprob,4,nsea)
                nd_west = nsea
                !!              
                IF((nd_north==-1).and.(nd_east==-1).and.  &
                     &  (nd_south==-1).and.(nd_west ==-1)) THEN
                   ji = ji+1       ! We are on earth far from sea (at least nprob points distant)... 
                   ! Put initial land values
                ELSE
                   CALL DECISION(ji,jj,ni,nj,X,maskv,mask2) 
                 ! We are on earth but at least a sea point is closer than nprob
                   ji = ji+1
                ENDIF
             ELSE
                ji = ji+1              ! We are on sea, just moving on...
             ENDIF
          ENDDO
          !!        
       ENDDO
       !!
    ENDDO
    !!
    !! First and last line sometimes show weird features...
    DO ji = 1, ni
       IF ( mask(ji,1)  == 0 ) X(ji,1)  = X(ji,2)
       IF ( mask(ji,nj) == 0 ) X(ji,nj) = X(ji,nj-1)
    END DO
    !!
    !! Smoothing the what's been done on land:
    !!
    IF ( lsmooth ) THEN
       DO kk = 1, 50
          !!
          DO jj = 2, nj-1
             !!
             DO ji = 2, ni-1
                IF ( mask(ji,jj) == 0 )   X(ji,jj) = &
                     &    0.35*X(ji,jj) + 0.65*0.125*( &
                     &    X(ji+1,jj)   + X(ji,jj+1)   + X(ji-1,jj)   + X(ji,jj-1) + &
                     &    X(ji+1,jj+1) + X(ji-1,jj+1) + X(ji-1,jj-1) + X(ji+1,jj-1)   )
             END DO
             !!
             IF (k_ew /= -1) THEN   ! we can use east-west periodicity
                IF ( mask(1,jj) == 0 )   X(1,jj) = &
                     &    0.35*X(1,jj) + 0.65*0.125*( &
                     &    X(2,jj)   + X(1,jj+1)   + X(ni-k_ew,jj)   + X(1,jj-1) + &
                     &    X(2,jj+1) + X(ni-k_ew,jj+1) + X(ni-k_ew,jj-1) + X(2,jj-1)   )
                IF ( mask(ni,jj) == 0 )   X(ni,jj) = &
                     &    0.35*X(ni,jj) + 0.65*0.125*( &
                     &    X(1+k_ew,jj)   + X(ni,jj+1)   + X(ni-1,jj)   + X(ni,jj-1) + &
                     &    X(1+k_ew,jj+1) + X(ni-1,jj+1) + X(ni-1,jj-1) + X(1+k_ew,jj-1)   )
             END IF
             !!
          END DO
          !!
       END DO
       !!
    END IF
    !!
    !!
    DEALLOCATE ( maskv, mask2, data_old )
    !!
#if defined debug
    PRINT *, 'DROWN: jinc =', jinc; PRINT *, ''
#endif
    !!
    END IF
    !!
END SUBROUTINE DROWN
!!
!!
!!
!!
!!
!!
!!
!!
SUBROUTINE DECISION(ix, iy, ni, nj, X, mask, mask2)
  !!
  !! Arguments :
  !! -----------
  INTEGER,                   INTENT(in)    :: ix, iy, ni, nj
  REAL(4), DIMENSION(ni,nj), INTENT(inout) :: X
  INTEGER, DIMENSION(ni,nj), INTENT(in)    :: mask
  INTEGER, DIMENSION(ni,nj), INTENT(inout) :: mask2
  !!
  !!
  !! Local :
  !!--------
  !!  
  INTEGER :: ji, jj, jfi, jk
  !!
  REAL(kind=4), PARAMETER :: rr=0.707
  REAL(kind=4) :: Tn, Te, Ts, Tw, Tne, Tse, Tsw, Tnw
  INTEGER      :: Wn, We, Ws, Ww, Wne, Wse, Wsw, Wnw
  REAL(kind=8) :: rden
  !!
  !!
  !!
  !!
  !!     Initialisations
  !!-----------------------
  ji  = ix
  jfi = ix
  jj  = iy
  !!
  !!
  !!----------------------------
  !! North :
  !!----------------------------
  IF ( mask(ji,jj+1)==1 ) THEN
     Tn=X(ji,jj+1)
     Wn=1
  ELSE
     Tn=0.
     Wn=0
  ENDIF
  !!
  !!----------------------------
  !! North-east :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == ni ) THEN
     jfi=0
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi+1,jj+1)==1 ) THEN
     Tne=X(jfi+1,jj+1)
     Wne=1
  ELSE
     Tne=0.
     Wne=0
  ENDIF
  jfi=ji
  !!
  !!----------------------------
  !! East :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == ni ) THEN
     jfi=0
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi+1,jj)==1 ) THEN
     Te=X(jfi+1,jj)
     We=1
  ELSE
     Te=0.
     We=0
  ENDIF
  jfi=ji
  !!
  !!----------------------------
  ! South-east :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == ni ) THEN
     jfi=0
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi+1,jj-1)==1 ) THEN
     Tse=X(jfi+1,jj-1)
     Wse=1
  ELSE
     Tse=0.
     Wse=0
  ENDIF
  jfi=ji
  !!
  !!----------------------------
  ! South :
  !!----------------------------
  IF ( mask(ji,jj-1)==1 ) THEN
     Ts=X(ji,jj-1)
     Ws=1
  ELSE
     Ts=0.
     Ws=0
  ENDIF
  !!
  !!----------------------------
  ! South-west :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == 1 ) THEN
     jfi=ni+1
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi-1,jj-1)==1 ) THEN
     Tsw=X(jfi-1,jj-1)
     Wsw=1
  ELSE
     Tsw=0.
     Wsw=0
  ENDIF
  jfi=ji
  !!
  !!----------------------------
  ! West :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == 1 ) THEN
     jfi=ni+1
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi-1,jj)==1 ) THEN
     Tw=X(jfi-1,jj)
     Ww=1
  ELSE
     Tw=0.
     Ww=0
  ENDIF
  jfi=ji
  !!
  !!----------------------------
  ! North-west :
  !!----------------------------
  !! East-west periodicity :
  IF ( ji == 1 ) THEN
     jfi=ni+1
  ELSE
     jfi=ji
  ENDIF
  !!
  IF ( mask(jfi-1,jj+1)==1 ) THEN
     Tnw=X(jfi-1,jj+1)
     Wnw=1
  ELSE
     Tnw=0.
     Wnw=0
  ENDIF
  jfi=ji
  !!
  !!------------------------
  !!
  !! Computation of the sea value to give to this earth point
  !!----------------------------------------------------------
  rden     = (Wn+We+Ws+Ww+rr*(Wne+Wse+Wsw+Wnw))
  !!
  X(ji,jj) = (Wn*Tn+We*Te+Ws*Ts+Ww*Tw+rr*(Wne*Tne+Wse*Tse+Wsw*Tsw+Wnw*Tnw))/rden
  !!
  !! Former mask point becomes sea point :
  mask2(ji,jj) = 1
  !!
  !!
END SUBROUTINE DECISION
!!
!!
!!
!!----------------------------------------------------------
!!
!!
!!
!!
SUBROUTINE TEST4SEA(ix, iy, ni, nj, X, mask, ndist, ndir, nres)
  !!
  !!#########################################################################
  !!
  !!
  !!  ix  : current x position                (integer)
  !!  iy  : current y position                (integer)
  !!  ni  : x dimension of array X            (integer)     
  !!  nj  : y dimension of array Y            (integer)
  !!  X   : currently treated array           (2D array)
  !!  mask : mask array, 1=sea, 0=land         (2D array)
  !!  ndist : distance (in points) to check     (integer)
  !!
  !!                                              1
  !!  ndir: direction to check   (integer)    4 - | - 2
  !!                                              3
  !!
  !!  nres: result ->  -1 = no sea point was found on the way neither at position ndist 
  !!                  0<n<ndist = the first sea point was found at position n
  !!
  !!
  !!#########################################################################
  !!
  INTEGER,                        INTENT(in) :: ix, iy, ni, nj, ndist, ndir
  INTEGER, DIMENSION(ni,nj),      INTENT(in) :: mask
  REAL(kind=4), DIMENSION(ni,nj), INTENT(in) :: X
  INTEGER,                        INTENT(out):: nres
  !!
  !! Local :
  INTEGER :: ji, jj, ncpt, jfi
  !!
  ji  = ix
  jfi = ix
  jj  = iy 
  !!
  nres = -1  
  ncpt = 0
  !!
  !!
  SELECT CASE (ndir)
     !!
  CASE (1)
     DO WHILE ((jj<iy+ndist).and.(nres==-1))
        jj   = jj+1 ;  ncpt = ncpt+1
        IF(mask(ji,jj)/=0) nres = ncpt
     ENDDO
     !!
  CASE(2)
     DO WHILE ((ji<ix+ndist).and.(nres==-1))
        IF(ji==ni) jfi=0
        ji   = ji+1 ; jfi  = jfi+1 ; ncpt = ncpt+1
        IF (mask(jfi,jj)/=0) nres = ncpt
     ENDDO
     !!
  CASE(3)
     DO WHILE ((jj>iy-ndist).and.(nres==-1))
        jj   = jj-1 ;  ncpt = ncpt+1
        IF(mask(ji,jj)/=0) nres = ncpt
     ENDDO
     !!
  CASE(4)
     DO WHILE ((ji>ix-ndist).and.(nres==-1))
        IF(ji==1) jfi=ni+1
        ji   = ji-1 ;  jfi  = jfi-1 ;  ncpt = ncpt+1
        IF(mask(jfi,jj)/=0) nres = ncpt
     ENDDO
     !!
  CASE DEFAULT
     PRINT *, 'drwn.F90: you should not see this!'; STOP
  END SELECT
  !!
  !!
END SUBROUTINE TEST4SEA
!!
!!
!!
!!
!!##########################################################################

 SUBROUTINE FILL_EXTRA_BANDS(k_ew, lx, ly, X, Y, DAT, lxp4, lyp4, XP4, YP4, DATP4)
    !!
    !!============================================================================
    !! Extending input arrays with an extraband of two points at north,south,east 
    !! and west boundaries.
    !!
    !! The extension is done thanks to Akima's exptrapolation method.
    !!
    !! East-west periodicity of global map is taken into account through 'k_ew' :
    !!
    !!
    !!  k_ew : east-west periodicity on the input file/grid
    !!         k_ew = -1  --> no periodicity
    !!         k_ew >= 0  --> periodicity with overlap of k_ew points
    !! 
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    INTEGER ,                       INTENT(in)  :: k_ew
    !! 
    INTEGER,                        INTENT(in)  :: lx, ly, lxp4, lyp4
    REAL(8), DIMENSION(lx,ly),     INTENT(in)  :: X, Y, DAT
    !!
    REAL(8), DIMENSION(lxp4,lyp4), INTENT(out) :: XP4, YP4, DATP4
    !!
    !! Local
    INTEGER :: ji, jj
    !!
    !!
    !!
    !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
    !!   --------------------------------------------------
    !!
    !! Initialising :
    !! --------------
    XP4   = 0.
    YP4   = 0.
    DATP4 = 0.
    !!
    !! Filling centers :
    !! -----------------
    XP4(3:lxp4-2, 3:lyp4-2)     = X(:,:)
    YP4(3:lxp4-2, 3:lyp4-2)     = Y(:,:)
    DATP4(3:lxp4-2, 3:lyp4-2)   = DAT(:,:)
    !!
    !!
    !! X array :
    !! ---------
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity of input file to
       !!                   ! fill extra bands :
       XP4( 1     , 3:lyp4-2) = X(lx - 1 - k_ew , :) - 360.
       XP4( 2     , 3:lyp4-2) = X(lx - k_ew     , :) - 360.
       XP4(lxp4   , 3:lyp4-2) = X( 2 + k_ew     , :) + 360.
       XP4(lxp4-1 , 3:lyp4-2) = X( 1 + k_ew     , :) + 360.
       !!
    ELSE
       !!
       !! Left side :
       XP4(2, 3:lyp4-2) = X(2,:) - (X(3,:) - X(1,:))
       XP4(1, 3:lyp4-2) = X(1,:) - (X(3,:) - X(1,:))
       !!
       !! Right side :
       XP4(lxp4-1, 3:lyp4-2) = X(lx-1,:) + X(lx,:) - X(lx-2,:)
       XP4(lxp4  , 3:lyp4-2) = X(lx,:)   + X(lx,:) - X(lx-2,:)
       !!
    END IF
    !!
    !!
    !! Bottom side :
    XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
    XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))
    !!
    !! Top side :
    XP4(:,lyp4-1) = XP4(:,lyp4-3) + XP4(:,lyp4-2) - XP4(:,lyp4-4)
    XP4(:,lyp4)   = XP4(:,lyp4-2) + XP4(:,lyp4-2) - XP4(:,lyp4-4)
    !!
    !!
    !!
    !! Y array :
    !! ---------
    !!
    !! Top side :
    YP4(3:lxp4-2, lyp4-1) = Y(:, ly-1) + Y(:,ly) - Y(:,ly-2)
    YP4(3:lxp4-2, lyp4)   = Y(:, ly)   + Y(:,ly) - Y(:,ly-2)
    !! Bottom side :
    YP4(3:lxp4-2, 2) = Y(:,2) - (Y(:,3) - Y(:,1))
    YP4(3:lxp4-2, 1) = Y(:,1) - (Y(:,3) - Y(:,1))
    !!
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity
       !!                   ! fill extra bands :
       YP4( 1     , :) = YP4(lx - 1 - k_ew + 2, :)
       YP4( 2     , :) = YP4(lx - k_ew     + 2, :)
       YP4(lxp4   , :) = YP4( 2 + k_ew     + 2, :)
       YP4(lxp4-1 , :) = YP4( 1 + k_ew     + 2, :)
       !!
    ELSE
       !!
       !! Left side :
       YP4(2, :) = YP4(4,:) - (YP4(5,:) - YP4(3,:))
       YP4(1, :) = YP4(3,:) - (YP4(5,:) - YP4(3,:))
       !! Right side :
       YP4(lxp4-1,:) = YP4(lxp4-3,:) + YP4(lxp4-2, :) - YP4(lxp4-4, :)
       YP4(lxp4,:)   = YP4(lxp4-2,:) + YP4(lxp4-2,:)  - YP4(lxp4-4, :)
       !!
    END IF
    !!
    !!
    !! Data array :
    !! ------------
    !!
    IF (k_ew /= -1) THEN   ! we can use east-west periodicity of input file to
       !!                   ! fill extra bands :
       DATP4( 1     , 3:lyp4-2) = DAT(lx - 1 - k_ew , :)
       DATP4( 2     , 3:lyp4-2) = DAT(lx - k_ew     , :)
       DATP4(lxp4   , 3:lyp4-2) = DAT( 2 + k_ew     , :)
       DATP4(lxp4-1 , 3:lyp4-2) = DAT( 1 + k_ew     , :)
       !!
       !!
    ELSE
       !!
       !! Left side :
       DO jj = 3, lyp4-2
          CALL extra_2_east(XP4(lxp4-4,jj),XP4(lxp4-3,jj),XP4(lxp4-2,jj),        &
               &          XP4(lxp4-1,jj),XP4(lxp4,jj),                         &
               &          DATP4(lxp4-4,jj),DATP4(lxp4-3,jj),DATP4(lxp4-2,jj),  &
               &          DATP4(lxp4-1,jj),DATP4(lxp4,jj) )  
       END DO
       !!
       !! Right side :
       DO jj = 3, lyp4-2
          CALL extra_2_west(XP4(5,jj),XP4(4,jj),XP4(3,jj),                    &
               &          XP4(2,jj),XP4(1,jj),                               &
               &          DATP4(5,jj),DATP4(4,jj),DATP4(3,jj),               &
               &          DATP4(2,jj),DATP4(1,jj) )  
       END DO
       !!
       !!
    END IF
    !!
    !!
    !! Top side :
    DO ji = 1, lxp4
       CALL extra_2_east(YP4(ji,lyp4-4),YP4(ji,lyp4-3),YP4(ji,lyp4-2),        &
            &          YP4(ji,lyp4-1),YP4(ji,lyp4),                         &
            &          DATP4(ji,lyp4-4),DATP4(ji,lyp4-3),DATP4(ji,lyp4-2),  &
            &          DATP4(ji,lyp4-1),DATP4(ji,lyp4) )  
    END DO
    !!
    !! Bottom side :
    DO ji = 1, lxp4
       CALL extra_2_west(YP4(ji,5),YP4(ji,4),YP4(ji,3),        &
            &          YP4(ji,2),YP4(ji,1),                    &
            &          DATP4(ji,5),DATP4(ji,4),DATP4(ji,3),    &
            &          DATP4(ji,2),DATP4(ji,1) )  
    END DO
    !!
    !!
  END SUBROUTINE FILL_EXTRA_BANDS
  !!
  !!
  !!
  SUBROUTINE extra_2_east(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
    !!
    !!============================================================================
    !!
    !! Extrapolates 2 extra east (or north) points of a curve with Akima's 1D method
    !!
    !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
    !! Output : y4, y5
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y1, y2, y3
    REAL(8), INTENT(out) :: y4, y5
    !!
    !! Local :
    REAL(8) :: A, B, C, D, ALF, BET
    !!
    !!
    A    = x2 - x1
    B    = x3 - x2
    C    = x4 - x3
    D    = x5 - x4
    !!
    ALF  = y2 - y1
    BET  = y3 - y2
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       y4 = y3 ; y5 = y3
    ELSE
       y4   = C*(2*BET/B - ALF/A) + y3
       y5   = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C 
    END IF
    !!
    !!
  END SUBROUTINE extra_2_east
  !!
  !!
  !!
  !!
  SUBROUTINE extra_2_west(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
    !!
    !!============================================================================
    !!
    !! Extrapolates 2 extra west (or south) points of a curve with Akima's 1D method
    !!
    !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
    !! Output : y4, y5
    !!
    !!                       Author : Laurent BRODEAU, 2007
    !!============================================================================
    !!
    !!
    REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y5, y4, y3
    REAL(8), INTENT(out) :: y1, y2
    REAL(8) :: A, B, C, D, ALF, BET
    !!
    !! x1 -> x5
    !! x2 -> x4
    !! x3 -> x3
    !! x4 -> x2
    !! x5 -> x1
    !!
    A    = x4 - x5
    B    = x3 - x4
    C    = x2 - x3
    D    = x1 - x2
    !!
    ALF  = y4 - y5
    BET  = y3 - y4
    !!
    IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
       y2 = y3; y1 = y3
    ELSE
       y2   = C*(2*BET/B - ALF/A) + y3
       y1   = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C 
    END IF
    !!
    !!
  END SUBROUTINE extra_2_west
  !!
  !!
  FUNCTION TEST_XYZ(rx, ry, rz)
    !!
    !! Testing if 2D coordinates or 1D, and if match shape of data...
    !!
    CHARACTER(len=2) :: TEST_XYZ
    !!
    REAL(4), DIMENSION(:,:), INTENT(in) :: rx, ry, rz
    !!
    INTEGER :: ix1, ix2, iy1, iy2, iz1, iz2
    !!
    ix1 = size(rx,1) ; ix2 = size(rx,2)
    iy1 = size(ry,1) ; iy2 = size(ry,2)
    iz1 = size(rz,1) ; iz2 = size(rz,2)
    !!
    IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN
       !!
       IF ( (ix1 == iz1).AND.(iy1 == iz2) ) THEN
          TEST_XYZ = '1d'
       ELSE
          PRINT *, 'ERROR, drwn.F90 = >TEST_XYZ 1 : longitude and latitude array do not match data!'
          PRINT *, ''; STOP
       END IF
       !!
    ELSE
       IF ( (ix1 == iz1).AND.(iy1 == iz1).AND.(ix2 == iz2).AND.(iy2 == iz2) ) THEN
          TEST_XYZ = '2d'
       ELSE
          PRINT *, 'ERROR, drwn.F90 = >TEST_XYZ 2 : longitude and latitude array do not match data!'
          PRINT *, ''; STOP
       END IF
    END IF
    !!
    !!
  END FUNCTION TEST_XYZ
  !!
  !!
  !!
END MODULE DRWN
