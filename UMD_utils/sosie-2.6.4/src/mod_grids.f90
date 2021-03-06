MODULE MOD_GRIDS

   USE mod_conf

   IMPLICIT none

   PRIVATE

   PUBLIC :: SRC_DOMAIN, &
      &    TRG_DOMAIN, &
      &    TERMINATE

   INTEGER :: &
      &     ji, jj, jt0, jz0, &
      &     n1, n2, n3, nrec,   &
      &     if0, iv0

   REAL(wpl) :: rmv, dx, dy

   LOGICAL :: lmval


CONTAINS




   SUBROUTINE SRC_DOMAIN()

      USE io_ezcdf  !: => only to get time array !LB???

      INTEGER :: jt

      !! Determine input dimensions from input file :
      CALL know_dim_in()

      !! Allocate input arrays with input dimensions :
      ALLOCATE ( data_in(ni_in,nj_in), mask_in(ni_in,nj_in,nk_in), data_in_b(ni_in,nj_in),    &
         &     mask_in_b(ni_in,nj_in,nk_in), vt0(Ntr0), vt(Ntr) )

      IF ( lregin ) THEN
         ALLOCATE ( lon_in(ni_in,1),     lat_in(nj_in,1)     )
      ELSE
         ALLOCATE ( lon_in(ni_in,nj_in), lat_in(ni_in,nj_in) )
      END IF

      IF ( l_int_3d ) ALLOCATE ( data3d_in(ni_in,nj_in,nk_in), depth_in(nk_in) )
      !! Filling time array :
      IF ( ltime )  THEN
         IF ( lct ) THEN       ! time is being controlled
            DO jt = 1, Ntr
               vt(jt) = t0 + t_stp*REAL(jt)
            END DO
         ELSE        ! we use time from input file
            CALL GETVAR_1D(cf_in, cv_t_in, vt0) ;  vt(:) = vt0(j_start:j_stop)
         END IF
      END IF

      CALL get_src_conf()

      max_lon_in = maxval(lon_in);   max_lat_in = maxval(lat_in)
      min_lon_in = minval(lon_in);   min_lat_in = minval(lat_in)

   END SUBROUTINE SRC_DOMAIN




   SUBROUTINE TRG_DOMAIN()

      CALL know_dim_out()

      !! Allocate output arrays with output dimensions :
      ALLOCATE ( mask_out(ni_out,nj_out,nk_out), data_out(ni_out,nj_out) )

      IF ( lregout ) THEN
         ALLOCATE ( lon_out(ni_out, 1) , lat_out(nj_out, 1), lon_out_b(ni_out, 1) )
      ELSE
         ALLOCATE ( lon_out(ni_out, nj_out) , lat_out(ni_out, nj_out), lon_out_b(ni_out, nj_out) )
      END IF

      IF (l_int_3d) THEN
         ALLOCATE ( depth_out(nk_out) )
         IF ( cmethod == 'no_xy' ) THEN
            ALLOCATE ( data3d_tmp(ni_in,  nj_in,  nk_in), data3d_out(ni_in,nj_in,nk_out) )
         ELSE
            ALLOCATE ( data3d_tmp(ni_out, nj_out, nk_in), data3d_out(ni_out,nj_out,nk_out) )
         END IF
      END IF

      CALL get_trg_conf()

      max_lat_out = maxval(lat_out) ;   min_lat_out = minval(lat_out) ;
      PRINT*,'' ; WRITE(6,*) 'Latitude max on input grid =', max_lat_in
      WRITE(6,*) 'Latitude max on output grid =', max_lat_out ; PRINT*,''

      !! Is output latitude increasing with j : 1 = yes | -1 = no
      nlat_inc_out = 1

      IF ( lregout ) THEN
         IF (lat_out(1,1) > lat_out(nj_out,1)) nlat_inc_out = -1
      ELSE
         IF (lat_out(1,1) > lat_out(1,nj_out)) nlat_inc_out = -1
      END IF

      jj_ex_top = 0 ; jj_ex_btm = 0

      !! Find first jj that exits max_lat_in --> a mettre ailleurs!
      IF ( lregout ) THEN
         IF  ( max_lat_out >  max_lat_in ) THEN
            jj_ex_top = (nlat_inc_out + 1)/2 + (1 - nlat_inc_out)/2*nj_out
            DO WHILE ( lat_out(jj_ex_top,1) < max_lat_in )
               jj_ex_top = jj_ex_top + nlat_inc_out
            END DO
         END IF
         !! Find first ji that exits max_lat_in --> a mettre ailleurs!
         IF  ( min_lat_out <  min_lat_in ) THEN
            jj_ex_btm = (nlat_inc_out + 1)/2*nj_out + (1 - nlat_inc_out)/2
            DO WHILE ( lat_out(jj_ex_btm,1) > min_lat_in )
               jj_ex_btm = jj_ex_btm - nlat_inc_out
            END DO
         END IF
      END IF

      WRITE(6,*) 'nj_out    =', nj_out    ; WRITE(6,*) 'jj_ex_top =', jj_ex_top
      WRITE(6,*) 'jj_ex_btm =', jj_ex_btm ; PRINT*,''

      IF (jj_ex_top > 0) jj_ex_top = jj_ex_top - nlat_inc_out
      IF (jj_ex_btm > 0) jj_ex_btm = jj_ex_btm + nlat_inc_out

   END SUBROUTINE TRG_DOMAIN





   SUBROUTINE TERMINATE()

      !data_in = 0.0; mask_in = 0.0; data_in_b = 0.0; lon_in = 0.0; lat_in = 0.0
      !mask_in_b = 0.0 ; vt0 = 0.0; vt = 0.0

      DEALLOCATE ( data_in, mask_in, data_in_b, lon_in, lat_in, mask_in_b, vt0, vt )

      !mask_out = 0.0; data_out = 0.0; lat_out = 0.0; lon_out = 0.0; lon_out_b = 0.0

      DEALLOCATE ( mask_out, data_out, lat_out, lon_out, lon_out_b )

      IF (l_int_3d) THEN
         !data3d_in = 0.0 ; depth_in = 0.0
         DEALLOCATE ( data3d_in, depth_in )
         !data3d_tmp = 0.0;  depth_out = 0.0;  data3d_out = 0.0
         DEALLOCATE ( data3d_tmp, depth_out, data3d_out )
      END IF

   END SUBROUTINE TERMINATE






   !! LOCAL SUBROUTINES
   !! ~~~~~~~~~~~~~~~~~

   SUBROUTINE get_src_conf()

      !! LOLO: Filling following array:
      !!       (depth_in)

      USE io_ezcdf
      USE mod_manip

      !! Local :
      REAL(wpl), DIMENSION(:,:,:), ALLOCATABLE :: z3d_tmp

      !! Getting grid on source domain:
      CALL rd_grid(-1, lregin, cf_x_in, cv_lon_in, cv_lat_in, lon_in, lat_in)

      IF ( l_int_3d ) THEN
         CALL rd_vgrid(nk_in, cf_z_in, cv_z_in, depth_in)
         WRITE(6,*) ''; WRITE(6,*) 'Input Depths ='; PRINT *, depth_in ; WRITE(6,*) ''
      END IF

      !! What about scale_factor / add_offset
      CALL GET_SF_AO(cf_in, cv_in, rsf, rao)
      WRITE(6,*) 'Scale factor =', rsf; WRITE(6,*) 'Add   offset =', rao; PRINT*,''

      IF ( lregin ) THEN
         !! Fixing input 1D longitude:
         CALL FIX_LONG_1D(ni_in, lon_in(:,1), nlon_inc_in, i_chg_lon)
         !! Fixing input 1D latitude:
         CALL FIX_LATD_1D(nj_in, lat_in(:,1), nlat_inc_in)
      ELSE
         WHERE ( lon_in < 0. )  lon_in = lon_in + 360.
      END IF



      !! Getting land-sea mask on source domain
      IF ( ldrown ) THEN

         IF ( trim(cf_lsm_in) == 'missing_value' ) THEN

            WRITE(6,*) 'Opening land-sea mask from missing_value on input data!'

            CALL WHO_IS_MV(cf_in, cv_in, ca_missval, rmv)

            ALLOCATE ( z3d_tmp(ni_in,nj_in,nk_in) )

            !! Read data field (at time 1 if time exists) :
            IF ( ltime  ) jt0 = 1

            IF ( l_int_3d ) THEN
               CALL GETVAR_3D(if0, iv0, cf_in, cv_in, Ntr,      jt0, z3d_tmp)
            ELSE
               IF ( l3d ) jz0 = jplev
               CALL GETVAR_2D(if0, iv0, cf_in, cv_in, Ntr, jz0, jt0, z3d_tmp(:,:,1))
            END IF

            mask_in = 1
            WHERE ( z3d_tmp == rmv ) mask_in = 0

            DEALLOCATE ( z3d_tmp )

         ELSE

            !! We are reading source land-sea mask from a file:
            !! -----------------------------------------------
            !! checking dimension:
            CALL DIMS(cf_lsm_in, cv_lsm_in, n1, n2, n3, nrec)

            IF ( l3d .AND. ( jplev > 1 ) ) THEN
               IF ( n3 == nk_in ) THEN
                  WRITE(6,*) 'Opening 3D land-sea mask on source grid for level', jplev
                  PRINT *, trim(cv_lsm_in)
                  CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1), jlev=jplev)
               ELSE
                  WRITE(6,*) 'PROBLEM! You want to interpolate level', jplev
                  WRITE(6,*) 'but your source land-sea mask is not 3D!'
                  WRITE(6,*) '=> set ldrown to false in the namelist'
                  WRITE(6,*) 'If you want to "drown" a level other than the surface,'
                  WRITE(6,*) 'please provide a 3D input land-sea mask'
                  STOP
               END IF
            ELSEIF ( l_int_3d ) THEN
               IF ( n3 == nk_in ) THEN
                  WRITE(6,*) 'Opening 3D land-sea mask file on source grid, ', trim(cv_lsm_in)
                  CALL GETMASK_3D(cf_lsm_in, cv_lsm_in, mask_in)
               ELSE
                  WRITE(6,*) 'We need to open the 3D source land-sea mask,'
                  WRITE(6,*) 'but the vertical dimension of it does not match!'
                  STOP
               END IF
            ELSE
               WRITE(6,*) 'Opening land-sea mask file on source grid, ', trim(cv_lsm_in)
               CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1))
            END IF
         END IF
      END IF ! IF ( ldrown )

      !! Need to modify the mask if lon or lat have been modified :
      IF ( nlat_inc_in == -1 ) CALL FLIP_UD_3D(mask_in)
      IF ( nlon_inc_in == -1 ) CALL LONG_REORG_3D(i_chg_lon, mask_in)

   END SUBROUTINE get_src_conf


   SUBROUTINE get_trg_conf()

      USE io_ezcdf

      REAL(wpl), DIMENSION(:,:,:), ALLOCATABLE :: z3d_tmp


      IF ( (lregout).AND.(trim(cf_x_out) == 'spheric') ) THEN

         !! Building target grid:
         READ(cv_lon_out,'(f5.2)') dx ; READ(cv_lat_out,'(f5.2)') dy
         WRITE(6,*) '  * dx, dy =', dx, dy
         WRITE(6,*) '  * ni_out, nj_out =', ni_out, nj_out ;  PRINT*,''
         DO ji = 1, ni_out
            lon_out(ji,1) = dx/2.0 + dx*(ji - 1)
         END DO
         DO jj = 1, nj_out
            lat_out(jj,1) = -90 + dy/2.0 + dy*(jj - 1)
         END DO

         WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_out; WRITE(6,*) ''
         WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_out; WRITE(6,*) ''; WRITE(6,*) ''


      ELSE
         !! Getting target grid from netcdf file:
         CALL rd_grid(ivect, lregout, cf_x_out, cv_lon_out, cv_lat_out, lon_out, lat_out)
         !!
         IF ( lregout ) THEN
            WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_out; WRITE(6,*) ''
            WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_out; WRITE(6,*) ''; WRITE(6,*) ''
         END IF
      END IF


      lon_out_b = lon_out
      WHERE ( lon_out < 0. ) lon_out = lon_out + 360.

      IF ( l_int_3d ) THEN
         IF ( trim(cf_x_out)  == 'spheric') THEN
            cf_z_out = cf_z_in ;  cv_z_out = cv_z_in         !Important
         END IF
         CALL rd_vgrid(nk_out, cf_z_out, cv_z_out, depth_out)
         WRITE(6,*) ''; WRITE(6,*) 'Output Depths ='; PRINT *, depth_out ; WRITE(6,*) ''
      END IF



      !!  Getting target mask (mandatory doing 3D interpolation!)
      IF ( lmout .OR. l_int_3d ) THEN

         IF ( (l3d).and.(jplev > 1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) 'We do not know output mask yet, since it is at a given depth!'
            WRITE(6,*) '--> next version of SOSIE!'
            WRITE(6,*) 'So we do not mask output file'
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) ''

            mask_out = 1

         ELSE

            IF ( trim(cf_lsm_out) == 'missing_value' ) THEN

               WRITE(6,*) 'Opening target land-sea mask from missing_value!'
               CALL WHO_IS_MV(cf_x_out, cv_lsm_out, ca_missval, rmv)

               ALLOCATE ( z3d_tmp(ni_out,nj_out,nk_out) )
               if0 = 0 ; iv0 = 0   ! Read data field (at time 1 if time exists)
               IF ( l_int_3d ) THEN
                  CALL GETVAR_3D(if0, iv0, cf_x_out, cv_lsm_out, 1, 1,    z3d_tmp)
               ELSE
                  CALL GETVAR_2D(if0, iv0, cf_x_out, cv_lsm_out, 1, 1, 1, z3d_tmp(:,:,1))
               END IF
               mask_out = 1
               WHERE ( z3d_tmp == rmv ) mask_out = 0
               DEALLOCATE ( z3d_tmp )

            ELSE

               IF ( l_int_3d ) THEN
                  IF ( ((trim(cf_lsm_out) == '').OR.(trim(cv_lsm_out) == '')) ) THEN
                     WRITE(6,*) 'WARNING: no target 3D land-sea mask provided (cf_lsm_out)!'
                     mask_out = 1
                  ELSE
                     WRITE(6,*) 'Opening 3D land-sea mask file on target grid: ',trim(cf_lsm_out)
                     WRITE(6,*) '             => name mask : ',trim(cv_lsm_out)
                     CALL GETMASK_3D(cf_lsm_out, cv_lsm_out, mask_out)
                     WRITE(6,*) ''
                  END IF
               ELSE
                  WRITE(6,*) 'Opening 2D land-sea mask file on target grid: ', trim(cf_lsm_out)
                  CALL GETMASK_2D(cf_lsm_out, cv_lsm_out, mask_out(:,:,1))
               END IF
               WRITE(6,*) ''
            END IF

         END IF

      END IF

   END SUBROUTINE get_trg_conf







   SUBROUTINE rd_grid(iv, lreg, cfgrd, cvx, cvy, rlon, rlat)

      USE io_ezcdf

      !!   AUTHOR:
      !!   -------
      !!   Laurent Brodeau, 2015
      !!
      !!   may 2007: modified by Pierre Mathiot to detect and handle
      !!               regular 2D input longitude and regular 2D input
      !!               latitude when the input grid is declared as irregular.
      !!               Allowing akima interpolation.

      !! Arguments:
      INTEGER,                   INTENT(in)  :: iv !: iv = -1 means we're handling input grid
      !!                                                /= -1 means target grid
      LOGICAL,                   INTENT(in)  :: lreg
      CHARACTER(len=400),        INTENT(in)  :: cfgrd
      CHARACTER(len=80),         INTENT(in)  :: cvx, cvy
      REAL(8),   DIMENSION(:,:), INTENT(out) :: rlon, rlat

      !! Local
      LOGICAL :: lreg2d, l2dyreg_x, l2dyreg_y
      CHARACTER(len=8) :: cdomain, clreg
      INTEGER :: &
         &     idx, lx, ly, &
         &     ii, ij, if1, iv1, &
         &     ilx1, ily1, ilz1,   &
         &     ilx2, ily2, ilz2
      REAL(8),   DIMENSION(:,:), ALLOCATABLE :: zrlon, zrlat

      IF ( iv == -1 ) THEN
         cdomain = 'source'
         clreg   = 'lregin'
         idx = 1
      ELSE
         cdomain = 'target'
         clreg   = 'lregout'
         idx = 2
      END IF

      IF ( lreg ) THEN
         IF ( (size(rlon,2) /= size(rlat,2)).OR.(size(rlon,2) /= 1) ) THEN
            WRITE(6,*) 'ERROR 1: rd_grid (prepare.F90)!' ; STOP
         END IF
         lx = size(rlon,1) ; ly = size(rlat,1)
      ELSE
         IF ( (size(rlon,1) /= size(rlat,1)).OR.(size(rlon,1) /= size(rlat,1)) ) THEN
            WRITE(6,*) 'ERROR 2: rd_grid (prepare.F90)!' ; STOP
         END IF
         lx = size(rlon,1) ; ly = size(rlon,2)
      END IF
      rlon = 0. ; rlat = 0.

      !! Checking the dimension of longitude variable
      CALL DIMS(cfgrd, cvx, ilx1, ily1, ilz1, nrec)
      CALL DIMS(cfgrd, cvy, ilx2, ily2, ilz2, nrec)

      WRITE(6,*) ''


      IF (lreg) THEN
         !!         -------------------------------
         !!            R E G U L A R   G R I D :
         !!         -------------------------------
         WRITE(6,*) 'Opening regular grid ', trim(cvx), ' , ', trim(cvy), &
            &                 ' in file ', trim(cfgrd), ' .'

         l2dyreg_x = .FALSE.
         IF ( (ilz1 /= -1).or.(ily1 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' longitude ', trim(cvx), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_x = .TRUE.
         END IF

         l2dyreg_y = .FALSE.
         IF ( (ilz2 /= -1).or.(ily2 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' latitude ', trim(cvy), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_y = .TRUE.
         END IF

         IF ( l2dyreg_x .AND. (.NOT. l2dyreg_y) ) THEN
            WRITE(6,*) 'Error! '//trim(cdomain)//' longitude and latidude do not agree in shape (1D vs 2D)!' ; STOP
         END IF

         IF ( l2dyreg_x .AND. l2dyreg_y ) THEN
            WRITE(6,*) ' '
            WRITE(6,*) ' =================================================================================================='
            WRITE(6,*) ' *** Assuming that '//trim(cdomain)//' grid is regular even though longitude and latidude are 2D!'
            WRITE(6,*) '                     (because you set '//trim(clreg)//'=.TRUE. in the namelist)'
            WRITE(6,*) ' =================================================================================================='
            l_2d_grid_yet_regular(idx) = .TRUE.
            WRITE(6,*) ' '
         END IF


         IF ( l_2d_grid_yet_regular(idx) ) THEN
            !! Getting regular 2D grid :
            ALLOCATE ( zrlon(lx,ly), zrlat(lx,ly) )
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 0, 0, zrlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 0, 0, zrlat)
            !!
            !! Checking for Regular 2D longitude and latitude
            lreg2d = .TRUE.
            !! Checking if longitude array changes with latitude
            DO ij = 2, ly
               IF ( ABS(SUM(zrlon(:,ij)-zrlon(:,1))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            !! Checking if latitude array changes with longitude
            DO ii = 2, lx
               IF ( ABS(SUM(zrlat(ii,:)-zrlat(1,:))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            IF ( lreg2d ) THEN
               WRITE(6,*) ' *** OK! You were right, '//trim(cdomain)//' longitude and latitude are 2D but the grid is regular!'
               WRITE(6,*) ''
            ELSE
               WRITE(6,*) 'Error! We checked, and 2D '//trim(cdomain)//' longitude/latitude do not qualify for a regular grid!!!!'
               WRITE(6,*) '       => so set '//trim(clreg)//' to .FALSE. !'
               STOP
            END IF
            !!
            !! Giving values to 1D arrays:
            rlon(:,1) = zrlon(:,2)
            rlat(:,1) = zrlat(2,:)
            DEALLOCATE ( zrlon, zrlat)
            !!
         ELSE
            !!
            !! Normal case: Getting regular 1D grid :
            IF ( (ilx1 /= lx).or.(ilx2 /= ly) ) THEN
               WRITE(6,*) 'Error! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), &
                  &   ' do not agree in dimension with configuration dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
            CALL GETVAR_1D(cfgrd, cvx, rlon(:,1))
            CALL GETVAR_1D(cfgrd, cvy, rlat(:,1))
         END IF


      ELSEIF (.NOT. lreg) THEN
         !!         ----------------------------------
         !!            I R R E G U L A R   G R I D :
         !!         ----------------------------------
         WRITE(6,*) 'Opening irregular grid ', trim(cvx), ' , ', trim(cvy), &
            &                 ' in file ', trim(cfgrd), ' .'
         WRITE(6,*) ''

         IF (ilx1 /= ilx2) THEN

            IF ( (ily1 == ily2).and.(ilz1 == ilz2).and.(ily1 == -1).and.(ilz1 == -1) ) THEN
               WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
                  &       trim(cvy), ' seem to be regular (1D), check namelist!'
               WRITE(6,*) 'In the file ', trim(cfgrd)
               WRITE(6,*) '       => so maybe try to set '//trim(clreg)//' to .TRUE. in the namelist!'
               STOP
            ELSE
               WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), ' do not agree in dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
         END IF

         IF ( (ily1 /= ily2).or.(ilz1 /= ilz2) ) THEN
            WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',    &
               &   trim(cvy), &
               &   ' do not agree in dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF
         IF ( (ilx1 /= lx).or.(ily1 /= ly) ) THEN
            WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
               &   trim(cvy),   &
               &   ' do not agree in dimension with configuration dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF


         !! Getting source longitude array at level=1 (surface) and time =1 :
         !! -----------------------------------------------------------------
         !! If 3d dimension, we chose first level
         IF ( ilz1 /= -1 ) THEN
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 1, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 1, 0, rlat)
         ELSE
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 0, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 0, 0, rlat)
         END IF

      END IF


   END SUBROUTINE rd_grid



   SUBROUTINE rd_vgrid(nk, cfgrd, cvz, vdepth)

      USE io_ezcdf

      INTEGER,                   INTENT(in) :: nk
      CHARACTER(len=400),        INTENT(in) :: cfgrd
      CHARACTER(len=80),         INTENT(in) :: cvz
      REAL(8), DIMENSION(nk) ,  INTENT(out) :: vdepth

      CALL GETVAR_1D(cfgrd, cvz, vdepth)

   END SUBROUTINE rd_vgrid





   SUBROUTINE know_dim_in()

      USE io_ezcdf

      INTEGER :: jk0, jrec

      nk_in = 1

      !! Determine input dimensions from input file :
      CALL DIMS(cf_in, cv_in, ni_in, nj_in, jk0, jrec)

      !PRINT *, 'mod_grids.f90, cf_in, cv_in =>', TRIM(cf_in), '', TRIM(cv_in)
      !PRINT *, 'ni_in, nj_in, jk0, jrec =>', ni_in, nj_in, jk0, jrec

      IF ( nj_in == -1 ) THEN
         WRITE(6,*) 'prepare.know_dim_in => ERROR! variable ',trim(cv_in),' should be at least 2D!!!'
         STOP
      END IF


      IF ( jplev == -1 ) THEN
         !! This is the overide case! Means 2D+T !
         !! Case when we want to read a 2D+T field but someone screwed up and the record
         !! dimension is not of 'UNLIMITED' type inside the netcdf file...
         !! So it just  overides good sence and force sosie to understand that
         !! your field to interpolate is 2D with a time record
         !! (usually the case if the time record dimension in your
         !! input file is not declared as UNLIMITED => bad! :(
         WRITE(6,*) ''
         Ntr = jk0
         WRITE(6,*) 'WARNING: know_dim_in of mod_grids.f90 !!!'
         WRITE(6,*) '   => we force input field "'//TRIM(cv_in)//'" to be 2D + time !!!'
         WRITE(6,*) '   => because you specified "jplev = -1" in the namelist!'
         WRITE(6,*) '   => the time-record dimension is therefore:', Ntr
         WRITE(6,*) ''

      ELSE

         IF ( jrec == -1) ltime=.FALSE.  ! no time records

         !! 3D variable
         IF ( jk0 > 0 ) THEN

            l3d = .TRUE.
            nk_in = jk0 ; WRITE(6,*) 'Number of level into input file =', nk_in

            IF (jplev /= 0) THEN
               IF ( (jplev <= 0).OR.(jplev > nk_in) ) THEN
                  WRITE(6,*) 'Level jplev is wrong! -> number of level =', nk_in
                  WRITE(6,*) 'In the namelist, jplev =', jplev ;  STOP
               END IF
               WRITE(6,*) 'We are going to interpolate level', jplev, ' of input file'
            ELSE
               WRITE(6,*) ''
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) 'WARNING: We are going to perform a 3D interpolation !!!'
               WRITE(6,*) '      => if this is not what you want, please specify '
               WRITE(6,*) '         which level (jplev) should be interpolated.'
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) ''
               l_int_3d = .TRUE.
            END IF

            IF ( jrec > 0 ) THEN
               !! 3D variable + time
               WRITE(6,*) 'Variable is 3D + time !'
               Ntr0 = jrec
            ELSE
               !! 3D variable without time
               WRITE(6,*) 'Variable is 3D without time axis'
               Ntr0 = 1
            END IF

         END IF


         !! 2D
         !! ~~~
         IF ( jk0 == -1 ) THEN

            l3d = .FALSE.

            IF ( jrec > 0 ) THEN
               !! 2D variable with time
               WRITE(6,*) 'Variable is 2D with time axis'
               Ntr0 = jrec
            ELSE
               !! 2D variable without time
               WRITE(6,*) 'Variable is 2D without time axis'
               Ntr0 = 1
            END IF

         END IF

         IF (.NOT. ltime) Ntr0 = 1

         Ntr = Ntr0;  j_start = 1 ; j_stop = Ntr0

      END IF

      j_start = 1 ; j_stop = Ntr


      IF ( (jt1 > 0).AND.(jt2 > 0) ) THEN
         Ntr = jt2 - jt1 + 1 ;  ;  j_start = jt1 ;  j_stop = jt2
         !! jrec is the time dimension of the input file, Ntr is the length requested by the user :
         IF ( ( Ntr > jrec ).OR.(jt1 < 1).OR.(jt1 > jt2).OR.(jt2 > jrec) ) THEN
            WRITE(6,*) ''; WRITE(6,*) 'Check jt1 and jt2 in the namelist:'
            WRITE(6, '("=> the time dimension of ",a," is ", i5)') TRIM(cv_in), jrec
            WRITE(6, '("   and you specified jt1, jt2 =",i5," ,",i5)') jt1, jt2
            STOP
         END IF
      END IF

   END SUBROUTINE know_dim_in






   SUBROUTINE know_dim_out
      !!
      USE io_ezcdf
      !!
      nk_out = 1
      !!
      !! If 3D interpolation and building spherical grid, we use levels from source grid:
      IF ( l_int_3d .AND. (trim(cf_x_out)  == 'spheric') ) THEN
         cf_z_out = cf_z_in ;  cv_z_out = cv_z_in
      END IF
      !!
      IF (lregout) THEN
         IF ( trim(cf_x_out) == 'spheric') THEN
            WRITE(6,*) ''; WRITE(6,*) 'Building regular spherical output grid!'
            READ(cv_lon_out,'(f5.2)') dx ; READ(cv_lat_out,'(f5.2)') dy
            ni_out = INT(360./dx) ; nj_out = INT(180./dy)
            GOTO 100
         ELSE
            CALL DIMS(cf_x_out, cv_lon_out, ni_out, n1, n2, nrec)
            CALL DIMS(cf_x_out, cv_lat_out, nj_out, n1, n2, nrec)
         END IF
      ELSE
         CALL DIMS(cf_x_out, cv_lon_out, ni_out, nj_out, n1, nrec)
      END IF
      !!
100   CONTINUE
      !!
      IF ( l_int_3d ) THEN
         !!
         WRITE(6,*) ''
         IF ( trim(cf_x_out)  == 'spheric' ) THEN
            WRITE(6,*) 'Since we are building our own spherical target grid,'
            WRITE(6,*) 'we are going to use source levels as target levels!'
         END IF
         WRITE(6,*) ''
         WRITE(6,*) ' => we read target levels in the following file:'
         PRINT *, trim(cf_z_out); WRITE(6,*) ''
         CALL DIMS(cf_z_out, cv_z_out, nk_out, n1, n2, nrec)
         WRITE(6,*) 'nk_out = ', nk_out ; WRITE(6,*) ''
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_out,'x',nj_out,'x',nk_out
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !!
      ELSE
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_out,'x', nj_out
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      END IF
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
      !!
   END SUBROUTINE know_dim_out
   !!
   !!
   !!
   !!
   !!
   SUBROUTINE FIX_LONG_1D(nx, vlon, n_inc, i_chg_x)
      !!
      !!
      !! Making 1D longitude understandable
      !! ==================================
      !!
      INTEGER,                  INTENT(in)    :: nx
      REAL(8), DIMENSION(nx),  INTENT(inout) :: vlon
      INTEGER,                  INTENT(out)   :: n_inc
      INTEGER,                  INTENT(out)   :: i_chg_x
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE    :: x1d_b
      INTEGER :: jx



      IF ( (vlon(1) == -180.0).AND.(vlon(nx) == 180.0) ) THEN
         !! In that very special case we do the following trick to avoid having
         !! twice the point "lon. 180" in the final longitude array
         vlon(1)  = -179.99
         vlon(nx) =  179.99
         IF ( ewper /= 0 ) THEN
            ewper = 0
            WRITE(6,*) 'prepare.F90: FIX_LONG_1D => "ewper" forced to 0 !'; WRITE(6,*) ''
         END IF
      END IF


      !! We want positive and increasing longitude !
      WHERE ( vlon < 0. )  vlon = vlon + 360.
      !!
      i_chg_x = 0 ; n_inc = 1
      DO jx = 2, nx
         IF ( vlon(jx) < vlon(jx-1) ) THEN
            IF ( i_chg_x /= 0) THEN  ! another changing sign point has already been found!
               WRITE(6,*) 'Your longitude array is a mess!'
               WRITE(6,*) ' --> Fix it (positive and increasing)!'; STOP
            ELSE
               i_chg_x = jx
               n_inc = -1
            END IF
         END IF
      END DO
      !!
      !!
      IF ( n_inc == -1 ) THEN
         !!
         ALLOCATE( x1d_b(nx) )
         x1d_b = vlon
         !!
         DO jx = i_chg_x, nx
            vlon(jx - i_chg_x + 1) = x1d_b(jx)
         END DO
         !!
         DO jx = 1, i_chg_x - 1
            vlon(nx - i_chg_x + 1 + jx) = x1d_b(jx)
         END DO
         !!
         DEALLOCATE( x1d_b )
         !!
      END IF
      !!
   END SUBROUTINE FIX_LONG_1D
   !!
   !!
   SUBROUTINE FIX_LATD_1D(ny, vlat, n_inc)
      !!
      !! Making 1D latitude increasing with jy
      !! =====================================
      !!
      INTEGER,                  INTENT(in)    :: ny
      REAL(8), DIMENSION(ny), INTENT(inout) :: vlat
      INTEGER,                  INTENT(out)   :: n_inc
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE :: y1d_b
      INTEGER :: jy
      !!
      n_inc = 1
      !!
      IF ( lregin .AND. ( vlat(1) > vlat(ny) ) ) THEN
         !!
         n_inc = -1 ; ALLOCATE( y1d_b(ny) )
         !!
         WRITE(6,*) ''; WRITE(6,*) 'Latitude does not seem to increase with j'
         WRITE(6,*) '--> We reverse 1D input latitude array!'
         y1d_b = vlat
         DO jy = 1, ny
            vlat(jy) = y1d_b(ny-jy+1)
         END DO
         !!
         DEALLOCATE( y1d_b )
         !!
      END IF
      !!
   END SUBROUTINE FIX_LATD_1D
   !!
   !!
   !!
END MODULE MOD_GRIDS
!!
