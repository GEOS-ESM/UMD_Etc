MODULE MOD_INIT

   USE mod_conf

   IMPLICIT NONE

   PUBLIC :: GET_ARGUMENTS, READ_NMLST, REMINDER


   !! Declaration of namelist :
   !! -------------------------

   NAMELIST /ninput/  ivect, lregin, cf_in, cv_in, cv_t_in, jt1, jt2, &
      &             jplev, cf_x_in, cv_lon_in, cv_lat_in, cf_lsm_in, cv_lsm_in, &
      &             ldrown, ewper, vmax, vmin, cf_coor_in

   NAMELIST /n3d/     cf_z_in, cv_z_in, cf_z_out, cv_z_out, cv_z_out_name

   NAMELIST /noutput/ lregout, cf_x_out, cv_lon_out, cv_lat_out, cf_lsm_out,   &
      &             cv_lsm_out, lmout, rmaskvalue, lct, t0, t_stp, ewper_out

   NAMELIST /nnetcdf/ cmethod, cv_l_out, cv_p_out, cv_t_out,    &
      &             cv_out, cu_out, cu_t, cln_out, cd_out,  &
      &             csource, ctarget, cextra, lpcknc4

   PRIVATE usage

CONTAINS


   SUBROUTINE GET_ARGUMENTS()

      INTEGER            :: iargc, jarg
      CHARACTER(len=400) :: cr
      CHARACTER(LEN=2), DIMENSION(3), PARAMETER :: clist_opt = (/ '-h','-p','-f' /)

      PRINT *, ''

      jarg = 0

      DO WHILE ( jarg < iargc() )

         jarg = jarg + 1
         CALL getarg(jarg,cr)

         SELECT CASE (trim(cr))

         CASE('-h')
            call usage()

         CASE('-f')

            IF ( jarg + 1 > iargc() ) THEN
               PRINT *, 'ERROR: Missing namelist name!' ; call usage()
            ELSE

               jarg = jarg + 1
               CALL getarg(jarg,cr)

               IF ( ANY(clist_opt == trim(cr)) ) THEN
                  PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
                  call usage()
               ELSE
                  cf_nml_sosie = trim(cr)
               END IF

            END IF

         CASE DEFAULT
            PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
            CALL usage()

         END SELECT

      END DO

      PRINT *, ''
      PRINT *, ' * namelist = ', trim(cf_nml_sosie)
      PRINT *

   END SUBROUTINE GET_ARGUMENTS


   SUBROUTINE READ_NMLST(iex)

      INTEGER, INTENT(in) :: iex !: 1 => sosie, 2 => corr_vect.x

      LOGICAL :: lexist

      INQUIRE(FILE=trim(cf_nml_sosie), EXIST=lexist )
      IF ( .NOT. lexist ) THEN
         WRITE(*,'("Namelist ",a," cannot be found!")') trim(cf_nml_sosie)
         CALL usage()
      END IF

      PRINT *, ''
      PRINT *, 'Opening namelist "', trim(cf_nml_sosie), '"'

      OPEN( UNIT=numnam, FILE=trim(cf_nml_sosie), FORM='FORMATTED', STATUS='OLD' )

      !! Reading input section:
      READ(numnam,ninput)

      !! Only reading 3D conf if 3D interpolation is required:
      IF (jplev == 0) THEN
         READ(numnam,n3d)
         IF ( trim(cv_z_out_name) == '' ) THEN
            cv_z_out_name = trim(cv_z_out) ; ! keep same name
         ELSE
            PRINT *, '*** Target depth vector in output file will be renamed to: ', trim(cv_z_out_name)
         END IF
      END IF

      !! Reading output section:
      READ(numnam,noutput)

      !! Reading netcdf section:
      READ(numnam,nnetcdf)

      CLOSE(numnam)



      IF ( iex == 1 ) THEN



         !! If this is a vector component, then interpolated files won't be correct
         !! until grid distorsion correction, so changing name of output variable :
         IF ( ivect /= 0 ) THEN
            lmout = .FALSE.         ! never mask for 1st stage of interpolation of the vector
            IF  ( ivect == 1 )  THEN
               cv_out = 'uraw'
            ELSEIF ( ivect == 2 )  THEN
               cv_out = 'vraw'
            ELSE
               PRINT *, 'Error: while interpolating a component of a vector field: specify&
                  & ivect=1 or ivect=2 into the namelist!'
            END IF
         END IF

         !! Output file name:
         !! -----------------

         IF ( lpcknc4 ) THEN
            cf_short= TRIM(cv_out)//'_'//TRIM(csource)//'-'//TRIM(ctarget)//   &
               &    '_'//TRIM(cextra)//'.nc4'
         ELSE
            cf_short= TRIM(cv_out)//'_'//TRIM(csource)//'-'//TRIM(ctarget)//   &
               &    '_'//TRIM(cextra)//'.nc'
         END IF

         cf_out = trim(cd_out)//'/'//trim(cf_short)

         !! Printing user's settings
         !!-------------------------
         !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)
         !!
         WRITE(6,*)''
         IF ( ivect /= 0 ) THEN
            WRITE(6,*)'The current field is a component of a vector!'
         ELSE
            WRITE(6,*)'The current field is NOT a component of a vector!'
         END IF
         WRITE(6,*)''
         !!
         WRITE(6,*)''
         IF ( lregin ) THEN
            WRITE(6,*)'Source grid is regular'
         ELSE
            WRITE(6,*)'Source grid is irregular'
         END IF
         WRITE(6,*)''
         !!
         !!
         WRITE(6,*)''
         IF ( lregout ) THEN
            WRITE(6,*)'Target grid is regular'
         ELSE
            WRITE(6,*)'Target grid is irregular'
         END IF
         WRITE(6,*)''
         !!
         WRITE(6,*)''
         WRITE(6,*)'Input source file: ', trim(cf_in)
         WRITE(6,*)''
         WRITE(6,*)'Method for 2D interoplation: ', cmethod
         WRITE(6,*)''
         WRITE(6,*)'Time record name in source file: ', trim(cv_t_in)
         WRITE(6,*)''
         WRITE(6,*)'File containing the source grid: ', trim(cf_x_in)
         WRITE(6,*)''
         WRITE(6,*)'Longitude name: ', trim(cv_lon_in)
         WRITE(6,*)''
         WRITE(6,*)'Latitude name: ', trim(cv_lat_in)
         WRITE(6,*)''
         WRITE(6,*)'Source grid mask: ', trim(cf_lsm_in)
         WRITE(6,*)''
         WRITE(6,*)'Source mask variable: ', trim(cv_lsm_in)
         WRITE(6,*)''
         WRITE(6,*)'Variable to be treated:', trim(cv_in)
         WRITE(6,*)''
         WRITE(6,*)'Level to treat:', jplev
         WRITE(6,*)''
         WRITE(6,*)'Name of variable on the output file:', trim(cv_out)
         WRITE(6,*)''
         WRITE(6,*)'Units of variable to be treated:', trim(cu_out)
         WRITE(6,*)''
         WRITE(6,*)'Long name of variable to be treated:', trim(cln_out)
         WRITE(6,*)''
         WRITE(6,*)'Target grid: ', trim(cf_x_out)
         !!
         IF ( (.NOT. lregout).and.(trim(cf_x_out) == 'spheric') ) THEN
            PRINT *, 'Problem! If you set "cf_x_out" to "spheric", then set "lregout" to ".TRUE."'
            PRINT *, ''
            STOP
         END IF
         !!
         WRITE(6,*)''
         WRITE(6,*)'Longitude on target grid: ', trim(cv_lon_out)
         WRITE(6,*)''
         WRITE(6,*)'Latitude on target grid: ', trim(cv_lat_out)
         WRITE(6,*)''
         WRITE(6,*)'Longitude name on output file: ', trim(cv_lon_out)
         WRITE(6,*)''
         WRITE(6,*)'Latitude name on output file: ', trim(cv_lat_out)
         WRITE(6,*)''
         WRITE(6,*)'Output target land-sea mask file: ', trim(cf_lsm_out)
         PRINT *, ''
         WRITE(6,*)'Name of land-sea mask variable on target grid: ', trim(cv_lsm_out)
         WRITE(6,*)''
         WRITE(6,*)'Output file: ', trim(cf_out)
         WRITE(6,*)''
         WRITE(6,*)'Shall we use DROWN: ', ldrown
         WRITE(6,*)''
         WRITE(6,*)'East west periodicity of source grid: ', ewper
         WRITE(6,*)''
         WRITE(6,*)'Masking output file: ', lmout
         WRITE(6,*)''
         WRITE(6,*)'Value for masked points in output file: ', rmaskvalue
         WRITE(6,*)''
         WRITE(6,*)'East west periodicity of target grid: ', ewper_out
         WRITE(6,*)''


         IF (jplev == 0) THEN
            WRITE(6,*)''
            WRITE(6,*)' Going to perform 3D interpolation (jplev = 0): '
            WRITE(6,*)''
            WRITE(6,*)' => file containing source depth vector: ', trim(cf_z_in)
            WRITE(6,*)' => name of source depth vector: ',         trim(cv_z_in)
            WRITE(6,*)' => file containing target depth vector: ', trim(cf_z_out)
            WRITE(6,*)' => name of target depth vector: ',         trim(cv_z_out)
            PRINT *,   ' => Target depth vector in output file will be renamed to: ', trim(cv_z_out_name)
            WRITE(6,*)''
         END IF

         IF (lct) THEN
            WRITE(6,*)'Time is controled ! '
            WRITE(6,*)'Start at: ', t0
            WRITE(6,*)'Time step is: ', t_stp
            WRITE(6,*)''
         ELSE
            WRITE(6,*)'Time is not controlled and will be the same than in input file! '
         END IF


         WRITE(6,*)''
         WRITE(6,*)' * Output file to compressed netcdf4 : ', lpcknc4


         WRITE(6,*)''
         WRITE(6,*)'                ---------------------'
         WRITE(6,*)''
         WRITE(6,*) ''
         WRITE(6,*) ''
         !CLOSE(6)


         !! Checking for some "non-sense"
         IF ( (cmethod == 'akima').and.(.NOT. lregin) ) THEN
            PRINT *, 'The Akima method only supports regular input grids!'
            PRINT *, '--> If the grid of the source domain is irregular, '
            PRINT *, '    please change "cmethod" from akima to "bilin" or "bicub" into the namelist!'
            STOP
         END IF

         WRITE(cpat,'(a,"-",a)') trim(csource), trim(ctarget)
         PRINT *, '';  PRINT *, 'Starting interpolation for config "',trim(cpat),'"...'
         PRINT *, '';  PRINT *, ''; PRINT *, ''

      END IF

   END SUBROUTINE READ_NMLST



   SUBROUTINE REMINDER()
      !!
      WRITE(6,*) '' ;  WRITE(6,*) ''
      WRITE(6,*) '====================================================='
      WRITE(6,*) '                    Current config:'
      WRITE(6,*) '                    ---------------'  !
      WRITE(6,'(" Input domain dimension          : ",i5,"  x",i5,"  x",i4)') &
         &        ni_in , nj_in, nk_in
      WRITE(6,'(" Output domain dimension         : ",i5,"  x",i5,"  x",i4)') &
         &        ni_out, nj_out, nk_out
      WRITE(6,'(" Number of time steps to proceed : ",i5)') Ntr
      IF (l_int_3d) THEN
         WRITE(6,*) 'Type of interpolation           :    3D'
      ELSE
         IF (l3d) THEN
            WRITE(6,'(" Type of interpolation           : 2D on level ",i3," of ",i4)') &
               jplev, nk_in
         ELSE
            WRITE(6,*) 'Type of interpolation           :    2D'
         END IF
      END IF
      WRITE(6,*) '====================================================='
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE REMINDER



   !! Local routine:

   SUBROUTINE usage()
      !!
      PRINT *,''
      PRINT *,'List of command line options:'
      PRINT *, ''
      PRINT *,' -f  <namelist_file>  => Specify which namelist file to use'
      PRINT *, ''
      PRINT *,' -h                   => Show this message'
      PRINT *, ''
      STOP
      !!
   END SUBROUTINE usage

END MODULE MOD_INIT
