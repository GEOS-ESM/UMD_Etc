program rst2socaice

  use netcdf
  use ocnice_utils

  implicit none

  integer, parameter :: real_kind = selected_real_kind(6)
  integer, parameter :: dbl_kind = selected_real_kind(13)
  integer :: i, j, n, nt, k, argc, fid_in, arg_index, Nargs
  character*1028 :: seaice_rst_fname  = 'seaicethermo_internal_rst'
  character*400  :: tilefile_fname  = 'tile.data'
  character*1028 :: SOCA_grid_fname  = 'soca_gridspec.nc'
  character*1028 :: Outfile  = 'cice_rst_OUT.nc4'

  character*128 :: fname
  integer :: varid, varidvi, varidlon, varidlat
  character*128 :: varname

  real*8, allocatable, dimension(  :, : , :)    :: DUMVAR

  integer, parameter   :: ncat = 5
  integer :: ncid
  integer :: x_dimid, y_dimid, z_dimid, t_dimid, dimids2d(3), dimids3d(4)

  real*8, allocatable, dimension(:, :) :: aicen, vicen ,vsnon
  real*8, allocatable, dimension(:, :, :) :: aicen2, vicen2 ,vsnon2
  real*8, allocatable, dimension(:) :: latq, lonq

  logical        :: append

  type (OceanRaster4SOCA)                       :: oraster   ! Ocean tiles
  type (ModelGrid4SOCA)                         :: ogrid     ! Ocean grid (soca_gridspec)
 

  !RASTER INFO
  !===========
  call readraster4SOCA(oraster)

  !OCEAN GRID
  !==========
  call read_grid_SOCA(ogrid)


!! READ RESTART FILE
 allocate( aicen(oraster%Ntiles,ncat) )
  varname = 'FR'
  call readrst2d(seaice_rst_fname, varname, aicen)

 allocate(vicen(oraster%Ntiles,ncat))
  varname = 'VOLICE'
  call readrst2d(seaice_rst_fname, varname, vicen)

allocate(vsnon(oraster%Ntiles,ncat))
  varname = 'VOLSNO'
  call readrst2d(seaice_rst_fname, varname, vsnon)

allocate( aicen2(oraster%IM2,oraster%JM2,ncat) )   !!In cice.rst.rc in SOCA: aicen(ncat,lat,lon)
allocate( vicen2(oraster%IM2,oraster%JM2,ncat) )
allocate( vsnon2(oraster%IM2,oraster%JM2,ncat) )

allocate( latq(ogrid%Ny) )
allocate( lonq(ogrid%Nx) )

latq = ogrid%yq
lonq = ogrid%xq

call remap_t2g_SOCA( oraster,aicen, aicen2 )
call remap_t2g_SOCA( oraster,vicen, vicen2 )
call remap_t2g_SOCA( oraster,vsnon, vsnon2 )
  print*, minval(aicen), maxval(aicen)
  print*, minval(aicen2), maxval(aicen2)

! WRITE OUTPUT: NEW CICE RESTART FILE

varname = 'aicen'
append=.false.
call write_soca_seaice_rst3d(Outfile, varname, aicen2, append, oraster)

varname = 'vicen'
append=.true.
call write_soca_seaice_rst3d(Outfile, varname, vicen2, append, oraster)

varname = 'vsnon'
append=.true.
call write_soca_seaice_rst3d(Outfile, varname, vsnon2, append, oraster)

varname = 'latq'
append=.true.
call write_soca_seaice_rst1dy(Outfile, varname, latq, append, oraster)

varname = 'lonq'
append=.true.
call write_soca_seaice_rst1dx(Outfile, varname, lonq, append, oraster)

end program rst2socaice
