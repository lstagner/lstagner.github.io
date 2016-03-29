!!FIDASIM Version 1.0.0

!The main routine (fidasim) is at the end of the file!
module simulation

USE H5LT !! High level HDF5 Interface
USE HDF5 !! Base HDF5
USE hdf5_extra !! Additional HDF5 routines
USE eigensystem, ONLY : eigen, matinv
USE parallel_rng

implicit none

integer, parameter :: long   = 4 !bytes = 32 bits (-2,147,483,648 to 2,147,483,647)
integer, parameter :: long64 = 8 !bytes = 64 bits (-9,223,372,036,854,775,808 to 9,223,372,036,854,775,807)
integer, parameter :: float  = 4 !bytes = 32 bits (1.2E-38 to 3.4E+38) at 6 decimal places
integer, parameter :: double = 8 !bytes = 64 bits (2.3E-308 to 1.7E+308) at 15 decimal places

character(120)     :: namelist_file
integer, parameter :: nbif_type  = 1 ! full energy NBI spectra/density
integer, parameter :: nbih_type  = 2 ! half energy NBI spectra/density
integer, parameter :: nbit_type  = 3 ! third energy NBI spectra/density
integer, parameter :: halo_type  = 4 ! halo spectra/density
integer, parameter :: fida_type  = 5 ! fida spectra/density
integer, parameter :: brems_type = 6 ! brems-strahlung
integer, parameter :: ntypes     = 6 ! number of different types of neutrals

integer, parameter :: beam_ion = 1
integer, parameter :: thermal_ion = 2

!!Physical units:
real(double), parameter :: e_amu = 5.485799093287202d-4
real(double), parameter :: H_1_amu = 1.00782504d0
real(double), parameter :: H_2_amu = 2.0141017778d0
real(double), parameter :: B5_amu = 10.81d0
real(double), parameter :: C6_amu = 12.011d0

real(double), parameter :: mass_u    = 1.6605402d-27  ! [kg]
real(double), parameter :: e0        = 1.60217733d-19 ! [C]
real(double), parameter :: pi        = 3.14159265358979323846264d0
real(double), parameter :: c0        = 2.99792458d+08 !! [m/s]
real(double), parameter :: h_planck  = 4.135667516d-15 !![eV/s]
real(double), parameter :: lambda0   = 6561.d0        !!D-alpha [A]
real(double), parameter :: v2_to_E_per_amu = mass_u/(2.*e0*1.d3)*1.d-4 !!conversion cm^2/s^2 to keV

!! ---- Stark splitting, wavelength and intenisty of all 15 lines ---- !!
integer, parameter ::n_stark = 15
real(double), parameter, dimension(n_stark) :: stark_wavel = &
     [-2.20200d-06,-1.65200d-06,-1.37700d-06,-1.10200d-06, &
      -8.26400d-07,-5.51000d-07,-2.75600d-07, 0.00000d0,   &
       2.75700d-07, 5.51500d-07, 8.27400d-07, 1.10300d-06, &
       1.38000d-06, 1.65600d-06, 2.20900d-06               ]
real(double), parameter, dimension(n_stark) :: stark_intens= &
     [ 1.000d0, 18.00d0, 16.00d0, 1681.d0, 2304.d0, &
       729.0d0, 1936.d0, 5490.d0, 1936.d0, 729.0d0, &
       2304.d0, 1681.d0, 16.00d0, 18.00d0, 1.000d0  ]
integer, parameter, dimension(n_stark) :: stark_pi= &
     [1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1]
integer, parameter, dimension(n_stark) :: stark_sigma=1 - stark_pi

!!Numerical Settings
integer, parameter :: nlevs=6             !!nr of quantum states
real(double), parameter :: n_halo_neutrate=20. !! to average halo neut-rate
real(double) :: colrad_threshold=1.d6 !! to speed up simulation!
real(double), dimension(ntypes) :: halo_iter_dens = 0.d0
integer :: nbi_outside=0

type BeamGrid
    !! Beam Grid settings
    integer(long) :: nx     !! Nr. of cells in x direction
    integer(long) :: ny     !! Nr. of cells in y direction
    integer(long) :: nz     !! Nr. of cells in z direction
    real(double)  :: xmin, xmax
    real(double)  :: ymin, ymax
    real(double)  :: zmin, zmax
    
    !! Tait-bryan angles for z-y'-x" rotation
    real(double)  :: alpha  !! rotation about z
    real(double)  :: beta   !! rotation about y'
    real(double)  :: gamma  !! rotation about x"
    real(double)  :: drmin  !! min(dx,dy,dz)
    real(double)  :: dv     !! volume of cells
    real(double)  :: volume !! Grid volume
    integer(long) :: ntrack !! Maximum Nr. of cells for tracking
    integer(long) :: ngrid  !! Nr. of cells
    real(double), dimension(3)   :: origin
    real(double), dimension(3)   :: center !! Center of grid
    real(double), dimension(3)   :: dr     !! dx, dy, dz
    real(double), dimension(3)   :: lwh    !! Grid length(x),width(y),height(z)
    real(double), dimension(3,3) :: basis
    real(double), dimension(3,3) :: inv_basis
    real(double), dimension(:), allocatable :: xc !! X centers
    real(double), dimension(:), allocatable :: yc !! Y centers
    real(double), dimension(:), allocatable :: zc !! Z centers
end type BeamGrid

type InterpolationGrid
    integer(long) :: nr !! Number of radii
    integer(long) :: nz !! Number of Z
    real(double)  :: dr !! Radial spacing
    real(double)  :: dz !! Vertical spacing
    real(double)  :: da !! Area
    real(double), dimension(:),   allocatable :: r  !! Radius
    real(double), dimension(:),   allocatable :: z   !! W/Z
    real(double), dimension(:,:), allocatable :: r2d !! 2D radius grid
    real(double), dimension(:,:), allocatable :: z2d !! 2D W grid
end type InterpolationGrid

type Profiles
    real(double) :: dene = 0.d0
    real(double) :: denp = 0.d0
    real(double) :: denimp = 0.d0
    real(double) :: denf = 0.d0
    real(double) :: te = 0.d0
    real(double) :: ti = 0.d0
    real(double) :: zeff = 0.d0
    real(double) :: vr = 0.d0
    real(double) :: vt = 0.d0
    real(double) :: vz = 0.d0
end type Profiles

type, extends( Profiles ) :: LocalProfiles
    logical :: in_plasma = .False.
    real(double), dimension(3) :: pos = [0.d0, 0.d0, 0.d0]
    real(double), dimension(3) :: vrot = [0.d0, 0.d0, 0.d0] !xyz
end type LocalProfiles

type EMFields
    real(double) :: br = 0.d0
    real(double) :: bt = 0.d0
    real(double) :: bz = 0.d0
    real(double) :: er = 0.d0
    real(double) :: et = 0.d0
    real(double) :: ez = 0.d0
end type EMFields

type, extends( EMFields ) :: LocalEMFields
    logical      :: in_plasma = .False.
    real(double) :: b_abs = 0.d0
    real(double) :: e_abs = 0.d0
    real(double), dimension(3) :: pos = [0.d0, 0.d0, 0.d0]
    real(double), dimension(3) :: a_norm = [0.d0, 0.d0, 0.d0]
    real(double), dimension(3) :: b_norm = [0.d0, 0.d0, 0.d0]
    real(double), dimension(3) :: c_norm = [0.d0, 0.d0, 0.d0]
    real(double), dimension(3) :: e_norm = [0.d0, 0.d0, 0.d0]
end type LocalEMFields

type Equilibrium
    type(EMFields), dimension(:,:), allocatable :: fields
    type(Profiles), dimension(:,:), allocatable :: plasma
    real(double), dimension(:,:), allocatable   :: mask
end type Equilibrium

type FastIonDistribution
    integer(long) :: nenergy !! Number of energies
    integer(long) :: npitch  !! Number of pitches
    real(double)  :: dE      !! Energy spacing
    real(double)  :: dp      !! Pitch spacing
    real(double)  :: emin
    real(double)  :: emax
    real(double)  :: e_range
    real(double)  :: pmin
    real(double)  :: pmax
    real(double)  :: p_range
    real(double), dimension(:), allocatable       :: energy  !! Energy array
    real(double), dimension(:), allocatable       :: pitch   !! Pitch array
    real(double), dimension(:,:,:,:), allocatable :: f       !! F(E,p,r,w)
end type FastIonDistribution

type FastIon
    logical       :: cross_grid = .False.
    real(double)  :: r = 0.d0 !cm
    real(double)  :: z = 0.d0 !cm
    real(double)  :: phi_enter = 0.d0 !radians
    real(double)  :: delta_phi = 0.d0 !radians
    real(double)  :: energy = 0.d0 !keV
    real(double)  :: pitch = 0.d0
    real(double)  :: vabs = 0.d0
    real(double)  :: vr = 0.d0 ! Radial velocity
    real(double)  :: vt = 0.d0 ! Torodial velocity
    real(double)  :: vz = 0.d0 ! Z velocity
    real(double)  :: weight = 0.d0
    integer(long) :: class = 0
end type FastIon

type FastIonParticles
    logical       :: guiding_center = .True.
    integer(long) :: nparticle = 0
    integer(long) :: nclass = 1
    type(FastIon), dimension(:), allocatable :: fast_ion
end type FastIonParticles

type NeutralBeam
    character(25) :: name = '' !! Beam name
    integer       :: shape     !! beam source shape 1="rectangular", 2="circular" 
    real(double)  :: widy      !! half height in y direction
    real(double)  :: widz      !! half width in z direction
    real(double)  :: focy      !! focal length in y direction
    real(double)  :: focz      !! focal length in z direction
    real(double)  :: einj      !! NBI voltage  [kV]
    real(double)  :: pinj      !! NBI power    [MW]
    real(double)  :: vinj      !! NBI velocity [cm/s]
    real(double)  :: alpha     !! Z rotation not same as beam_grid%alpha
    real(double)  :: beta      !! Tilt rotation not same as beam_grid%beta
    real(double), dimension(3)   :: divy    !! divergence in y direction
    real(double), dimension(3)   :: divz    !! divergence in z direction
    real(double), dimension(3)   :: species_mix
    real(double), dimension(3)   :: src     !! position of source
    real(double), dimension(3)   :: axis    !! position along beam sightline
    real(double), dimension(3,3) :: basis   !! Basis of neutral beam
    real(double), dimension(3,3) :: inv_basis !! Inverse of basis
end type NeutralBeam

type AtomicCrossSection
    integer      :: nenergy = 1
    real(double) :: logemin = 0.d0
    real(double) :: logemax = 0.d0 
    integer      :: n_max = nlevs
    integer      :: m_max = nlevs
    real(double) :: dlogE = 0.d0
    real(double) :: minlog_cross
    real(double), dimension(:,:,:), allocatable :: log_cross
end type AtomicCrossSection

type AtomicRates
    integer      :: nenergy = 1
    real(double) :: logemin = 0.d0
    real(double) :: logemax = 0.d0
    integer      :: ntemp = 1
    real(double) :: logtmin = 0.d0
    real(double) :: logtmax = 0.d0
    integer      :: n_max = nlevs
    integer      :: m_max = nlevs
    real(double) :: dlogE = 0.d0
    real(double) :: dlogT = 0.d0
    real(double) :: minlog_pop = 0.d0 
    real(double) :: minlog_depop = 0.d0 
    real(double), dimension(2) :: ab = 0.d0
    real(double), dimension(:,:,:,:,:), allocatable :: log_pop
    real(double), dimension(:,:,:,:), allocatable :: log_depop
end type AtomicRates

type AtomicTables
    type(AtomicCrossSection) :: H_H_cx
    type(AtomicRates)        :: H_H
    type(AtomicRates)        :: H_e
    type(AtomicRates)        :: H_Aq
    real(double), dimension(nlevs,nlevs) :: einstein
end type AtomicTables

type LineOfSight
    real(double) :: sigma_pi = 1.d0
    real(double) :: spot_size = 0.d0 ! Radius of spot size
    real(double), dimension(3) :: lens = [0.d0, 0.d0, 0.d0] ! xyz pos
    real(double), dimension(3) :: axis = [0.d0, 0.d0, 0.d0] ! xyz dir
end type LineOfSight

type SpectralChords
    integer :: nchan = 0
    type(LineOfSight), dimension(:), allocatable  :: los
    real(double), dimension(:), allocatable       :: radius
    logical, dimension(:,:,:), allocatable        :: los_inter
    real(double), dimension(:,:,:,:), allocatable :: dlength
end type SpectralChords

type BoundedPlane
    integer                      :: shape    = 0    ! 1="Rectangular", 2="circular"
    real(double)                 :: hh       = 0.d0 ! Half height [cm]
    real(double)                 :: hw       = 0.d0 ! Half width [cm]
    real(double), dimension(3)   :: origin   = 0.d0 ! Origin of plane in machine coordinates
    real(double), dimension(3,3) :: basis    = 0.d0 ! Basis vectors basis(:,1) = u_1 is plane normal
    real(double), dimension(3,3) :: inv_basis= 0.d0 ! Inverse basis
end type BoundedPlane

type NPADetector
    type(BoundedPlane) :: detector
    type(BoundedPlane) :: aperture
end type NPADetector

type NPAProbability
    real(double) :: p = 0.d0
    real(double), dimension(3) :: eff_rd = [0.d0, 0.d0, 0.d0] !Effective position of detector
end type NPAProbability

type NPAChords
    integer :: nchan = 0
    type(NPADetector), dimension(:), allocatable          :: det
    real(double), dimension(:), allocatable               :: radius
    logical, dimension(:,:,:), allocatable                :: hit
    type(NPAProbability), dimension(:,:,:,:), allocatable :: phit ! probability of hitting detector
end type NPAChords

type NPAParticle
    integer      :: detector = 0
    real(double) :: xi = 0.d0
    real(double) :: yi = 0.d0
    real(double) :: zi = 0.d0
    real(double) :: xf = 0.d0
    real(double) :: yf = 0.d0
    real(double) :: zf = 0.d0
    real(double) :: weight = 0.d0
    real(double) :: energy = 0.d0
    real(double) :: pitch = 0.d0
end type NPAParticle

type NPAResults
    integer(long) :: nchan = 0
    integer(long) :: npart = 0
    integer(long) :: nmax = 1000000
    integer(long) :: nloop = 1000
    type(NPAParticle), dimension(:), allocatable :: part
    real(double), dimension(:), allocatable      :: energy !! energy array
    real(double), dimension(:,:), allocatable    :: flux !! flux
end type NPAResults

type BirthProfile
    integer :: ind = 1
    real(double), dimension(:,:), allocatable     :: ri
    real(double), dimension(:,:), allocatable     :: vi
    real(double), dimension(:,:,:,:), allocatable :: dens
end type BirthProfile

type Spectra
    real(double), dimension(:,:), allocatable   :: brems
    real(double), dimension(:,:,:), allocatable :: bes
    real(double), dimension(:,:,:), allocatable :: fida
end type Spectra

type NeutralDensity
    real(double), dimension(:,:,:,:,:), allocatable :: dens
end type NeutralDensity

type FIDAWeights
    real(double), dimension(:,:), allocatable     :: fida
    real(double), dimension(:,:,:), allocatable   :: mean_f
    real(double), dimension(:,:,:,:), allocatable :: weight
end type FIDAWeights

type NPAWeights
    real(double), dimension(:,:,:,:,:), allocatable :: attenuation
    real(double), dimension(:,:,:,:,:), allocatable :: cx
    real(double), dimension(:,:,:,:), allocatable   :: emissivity
    real(double), dimension(:,:,:), allocatable     :: weight
    real(double), dimension(:,:), allocatable       :: flux
end type NPAWeights

type SimulationInputs
    integer(long)  :: shot_number
    real(double)   :: time
    character(120) :: runid = ''
    character(10)  :: version = ''
    character(120) :: result_dir = ''
    character(120) :: tables_file = ''
    character(120) :: geometry_file = ''
    character(120) :: equilibrium_file = ''
    character(120) :: distribution_file = ''
    character(120) :: neutrals_file = ''
    !! Monte Carlo settings
    integer(long)  :: n_fida
    integer(long)  :: n_npa
    integer(long)  :: n_nbi
    integer(long)  :: n_dcx
    integer(long)  :: n_halo
    integer(long)  :: n_birth
    !! Simulation switches
    integer(long)  :: calc_spec
    integer(long)  :: calc_brems
    integer(long)  :: calc_bes
    integer(long)  :: calc_fida
    integer(long)  :: load_neutrals
    integer(long)  :: calc_npa
    integer(long)  :: calc_fida_wght
    integer(long)  :: calc_npa_wght
    integer(long)  :: calc_birth
    integer(long)  :: dump_dcx
    integer(long)  :: verbose
    !! Neutral Beam Settings
    real(double)   :: ab   !! atomic mass of beam neutrals
    !! Plasma parameters
    integer(long)  :: impurity_charge !! charge of impurity
    real(double)   :: ai   !! atomic mass of plasma ions
    !! Distribution settings
    integer(long)  :: dist_type
    !! Spectrum parameters
    integer(long)  :: nlambda
    real(double)   :: dlambda
    real(double)   :: lambdamin
    real(double)   :: lambdamax
    !! Weight function settings
    integer(long)  :: ne_wght
    integer(long)  :: np_wght
    integer(long)  :: nphi_wght
    integer(long)  :: nlambda_wght
    real(double)   :: emax_wght
    real(double)   :: lambdamin_wght
    real(double)   :: lambdamax_wght
end type SimulationInputs

type ParticleTrack
    real(double) :: time = 0.d0
    real(double) :: flux = 0.d0
    integer(long), dimension(3) :: ind = [0,0,0]
    real(double), dimension(3)  :: pos = [0.d0,0.d0,0.d0]
end type ParticleTrack

!! operator overloading interface definitions
!! this allows for adding,multipying,and dividing 
!! LocalProfiles/LocalEMFields types
interface assignment(=)
    module procedure pp_assign, lpp_assign, plp_assign, lplp_assign, &
                     ff_assign, lff_assign, flf_assign, lflf_assign, &
                     fast_ion_assign,npa_part_assign
end interface

interface operator(+)
    module procedure pp_add,lplp_add,ff_add,lflf_add
end interface

interface operator(-)
    module procedure pp_subtract,lplp_subtract,ff_subtract,lflf_subtract
end interface

interface operator(*)
    module procedure sp_multiply, ps_multiply, lps_multiply, slp_multiply, &
                     sf_multiply, fs_multiply, lfs_multiply, slf_multiply
end interface

interface operator(/)
    module procedure ps_divide, lps_divide, fs_divide, lfs_divide
end interface

interface interpol_coeff
    module procedure interpol1D_coeff, interpol1D_coeff_arr
    module procedure interpol2D_coeff, interpol2D_coeff_arr
end interface

interface interpol
    module procedure interpol1D_arr
    module procedure interpol2D_arr, interpol2D_2D_arr
end interface

!! definition of the structures:
type(BeamGrid)            :: beam_grid
type(InterpolationGrid)   :: inter_grid
type(FastIonDistribution) :: fbm
type(FastIonParticles)    :: particles
type(Equilibrium)         :: equil
type(NeutralBeam)         :: nbi
type(AtomicTables)        :: tables
type(NPAResults)          :: npa
type(SpectralChords)      :: spec_chords
type(NPAChords)           :: npa_chords
type(SimulationInputs)    :: inputs
type(BirthProfile)        :: birth
type(NeutralDensity)      :: neut
type(Spectra)             :: spec
type(FIDAWeights)         :: fweight
type(NPAWeights)          :: nweight

contains

subroutine print_banner()
    write(*,'(a)') "   ____ ____ ___   ___    ____ ____ __  ___"
    write(*,'(a)') "  / __//  _// _ \ / _ |  / __//  _//  |/  /"
    write(*,'(a)') " / _/ _/ / / // // __ | _\ \ _/ / / /|_/ / "
    write(*,'(a)') "/_/  /___//____//_/ |_|/___//___//_/  /_/  "
    write(*,'(a)') "                                           "
    write(*,'(a)') "Version: 1.0.0"
    write(*,'(a)') ""
    write(*,'(a)') "FIDASIM is released as open source code under the MIT Licence."
    write(*,'(a)') "For more information visit ENTER WEBSITE HERE"
    write(*,'(a)') ""
end subroutine print_banner
                                     
!============================================================================
!---------------------------Operator Overloading-----------------------------
!============================================================================
subroutine fast_ion_assign(p1, p2)
    type(FastIon), intent(in)  :: p2
    type(FastIon), intent(out) :: p1
  
    p1%cross_grid = p2%cross_grid  
    p1%r          = p2%r
    p1%z          = p2%z
    p1%phi_enter  = p2%phi_enter
    p1%delta_phi  = p2%delta_phi
    p1%energy     = p2%energy
    p1%pitch      = p2%pitch
    p1%vabs       = p2%vabs
    p1%vr         = p2%vr 
    p1%vt         = p2%vt 
    p1%vz         = p2%vz 
    p1%weight     = p2%weight
    p1%class      = p2%class

end subroutine fast_ion_assign

subroutine npa_part_assign(p1, p2)
    type(NPAParticle), intent(in)  :: p2
    type(NPAParticle), intent(out) :: p1
  
    p1%xi = p2%xi 
    p1%yi = p2%yi
    p1%zi = p2%zi 
    p1%xf = p2%xf
    p1%yf = p2%yf
    p1%zf = p2%zf
    p1%weight = p2%weight
    p1%energy = p2%energy
    p1%pitch = p2%pitch
    p1%detector = p2%detector

end subroutine npa_part_assign

subroutine pp_assign(p1, p2)
    type(Profiles), intent(in)  :: p2
    type(Profiles), intent(out) :: p1
  
    p1%dene   = p2%dene  
    p1%ti     = p2%ti    
    p1%te     = p2%te    
    p1%denp   = p2%denp  
    p1%denf   = p2%denf  
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff  
    p1%vr     = p2%vr 
    p1%vt     = p2%vt 
    p1%vz     = p2%vz

end subroutine pp_assign

subroutine lpp_assign(p1, p2)
    type(Profiles), intent(in)       :: p2
    type(LocalProfiles), intent(out) :: p1
  
    p1%dene   = p2%dene  
    p1%ti     = p2%ti    
    p1%te     = p2%te    
    p1%denp   = p2%denp  
    p1%denf   = p2%denf  
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff  
    p1%vr     = p2%vr 
    p1%vt     = p2%vt 
    p1%vz     = p2%vz

end subroutine lpp_assign

subroutine plp_assign(p1, p2)
    type(LocalProfiles), intent(in) :: p2
    type(Profiles), intent(out)     :: p1
  
    p1%dene   = p2%dene  
    p1%ti     = p2%ti    
    p1%te     = p2%te    
    p1%denp   = p2%denp  
    p1%denf   = p2%denf  
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff  
    p1%vr     = p2%vr 
    p1%vt     = p2%vt 
    p1%vz     = p2%vz

end subroutine plp_assign

subroutine lplp_assign(p1, p2)
    type(LocalProfiles), intent(in)  :: p2
    type(LocalProfiles), intent(out) :: p1
  
    p1%pos    = p2%pos
    p1%dene   = p2%dene  
    p1%ti     = p2%ti    
    p1%te     = p2%te    
    p1%denp   = p2%denp  
    p1%denf   = p2%denf  
    p1%denimp = p2%denimp
    p1%zeff   = p2%zeff  
    p1%vr     = p2%vr 
    p1%vt     = p2%vt 
    p1%vz     = p2%vz
    p1%vrot   = p2%vrot

end subroutine lplp_assign

subroutine ff_assign(p1, p2)
    type(EMFields), intent(in)  :: p2
    type(EMFields), intent(out) :: p1
  
    p1%br   = p2%br  
    p1%bt   = p2%bt  
    p1%bz   = p2%bz  
    p1%er   = p2%er  
    p1%et   = p2%et  
    p1%ez   = p2%ez  

end subroutine ff_assign

subroutine lff_assign(p1, p2)
    type(EMFields), intent(in)       :: p2
    type(LocalEMFields), intent(out) :: p1
  
    p1%br   = p2%br  
    p1%bt   = p2%bt  
    p1%bz   = p2%bz  
    p1%er   = p2%er  
    p1%et   = p2%et  
    p1%ez   = p2%ez  

end subroutine lff_assign

subroutine flf_assign(p1, p2)
    type(LocalEMFields), intent(in) :: p2
    type(EMFields), intent(out)     :: p1
  
    p1%br   = p2%br  
    p1%bt   = p2%bt  
    p1%bz   = p2%bz  
    p1%er   = p2%er  
    p1%et   = p2%et  
    p1%ez   = p2%ez  

end subroutine flf_assign

subroutine lflf_assign(p1, p2)
    type(LocalEMFields), intent(in)  :: p2
    type(LocalEMFields), intent(out) :: p1
  
    p1%pos  = p2%pos
    p1%br   = p2%br  
    p1%bt   = p2%bt  
    p1%bz   = p2%bz  
    p1%er   = p2%er  
    p1%et   = p2%et  
    p1%ez   = p2%ez
    p1%b_abs = p2%b_abs
    p1%e_abs = p2%e_abs
    p1%a_norm = p2%a_norm
    p1%b_norm = p2%b_norm
    p1%c_norm = p2%c_norm
    p1%e_norm = p2%e_norm  

end subroutine lflf_assign

function pp_add(p1, p2) result (p3)
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3
  
    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%denp   = p1%denp   + p2%denp
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz

end function pp_add

function pp_subtract(p1, p2) result (p3)
    type(Profiles), intent(in) :: p1,p2
    type(Profiles)             :: p3
  
    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%denp   = p1%denp   - p2%denp
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz

end function pp_subtract

function lplp_add(p1, p2) result (p3)
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3
  
    p3%pos    = p1%pos    + p2%pos
    p3%dene   = p1%dene   + p2%dene
    p3%ti     = p1%ti     + p2%ti
    p3%te     = p1%te     + p2%te
    p3%denp   = p1%denp   + p2%denp
    p3%denf   = p1%denf   + p2%denf
    p3%denimp = p1%denimp + p2%denimp
    p3%zeff   = p1%zeff   + p2%zeff
    p3%vr     = p1%vr     + p2%vr
    p3%vt     = p1%vt     + p2%vt
    p3%vz     = p1%vz     + p2%vz
    p3%vrot   = p1%vrot   + p2%vrot

end function lplp_add

function lplp_subtract(p1, p2) result (p3)
    type(LocalProfiles), intent(in) :: p1,p2
    type(LocalProfiles)             :: p3
  
    p3%pos    = p1%pos    - p2%pos
    p3%dene   = p1%dene   - p2%dene
    p3%ti     = p1%ti     - p2%ti
    p3%te     = p1%te     - p2%te
    p3%denp   = p1%denp   - p2%denp
    p3%denf   = p1%denf   - p2%denf
    p3%denimp = p1%denimp - p2%denimp
    p3%zeff   = p1%zeff   - p2%zeff
    p3%vr     = p1%vr     - p2%vr
    p3%vt     = p1%vt     - p2%vt
    p3%vz     = p1%vz     - p2%vz
    p3%vrot   = p1%vrot   - p2%vrot

end function lplp_subtract

function ps_multiply(p1, real_scalar) result (p3)
    type(Profiles), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(Profiles)             :: p3
  
    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%denp   = p1%denp   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar 
    p3%vt     = p1%vt     * real_scalar 
    p3%vz     = p1%vz     * real_scalar 

end function ps_multiply

function sp_multiply(real_scalar, p1) result (p3)
    type(Profiles), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(Profiles)             :: p3
  
    p3 = p1*real_scalar

end function sp_multiply

function ps_divide(p1, real_scalar) result (p3)
    type(Profiles), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(Profiles)             :: p3
    
    p3 = p1*(1.d0/real_scalar)

end function ps_divide

function lps_multiply(p1, real_scalar) result (p3)
    type(LocalProfiles), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalProfiles)             :: p3
  
    p3%pos    = p1%pos    * real_scalar
    p3%dene   = p1%dene   * real_scalar
    p3%ti     = p1%ti     * real_scalar
    p3%te     = p1%te     * real_scalar
    p3%denp   = p1%denp   * real_scalar
    p3%denf   = p1%denf   * real_scalar
    p3%denimp = p1%denimp * real_scalar
    p3%zeff   = p1%zeff   * real_scalar
    p3%vr     = p1%vr     * real_scalar 
    p3%vt     = p1%vt     * real_scalar 
    p3%vz     = p1%vz     * real_scalar 
    p3%vrot   = p1%vrot   * real_scalar 

end function lps_multiply

function slp_multiply(real_scalar, p1) result (p3)
    type(LocalProfiles), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalProfiles)             :: p3
  
    p3 = p1*real_scalar

end function slp_multiply

function lps_divide(p1, real_scalar) result (p3)
    type(LocalProfiles), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalProfiles)             :: p3
    
    p3 = p1*(1.d0/real_scalar)

end function lps_divide

function ff_add(p1, p2) result (p3)
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3
  
    p3%br   = p1%br   + p2%br
    p3%bt   = p1%bt   + p2%bt
    p3%bz   = p1%bz   + p2%bz
    p3%er   = p1%er   + p2%er
    p3%et   = p1%et   + p2%et
    p3%ez   = p1%ez   + p2%ez

end function ff_add

function ff_subtract(p1, p2) result (p3)
    type(EMFields), intent(in) :: p1,p2
    type(EMFields)             :: p3
  
    p3%br   = p1%br   - p2%br
    p3%bt   = p1%bt   - p2%bt
    p3%bz   = p1%bz   - p2%bz
    p3%er   = p1%er   - p2%er
    p3%et   = p1%et   - p2%et
    p3%ez   = p1%ez   - p2%ez

end function ff_subtract

function fs_multiply(p1, real_scalar) result (p3)
    type(EMFields), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(EMFields)             :: p3
  
    p3%br   = p1%br   * real_scalar 
    p3%bt   = p1%bt   * real_scalar 
    p3%bz   = p1%bz   * real_scalar 
    p3%er   = p1%er   * real_scalar 
    p3%et   = p1%et   * real_scalar 
    p3%ez   = p1%ez   * real_scalar 

end function fs_multiply

function sf_multiply(real_scalar, p1) result (p3)
    type(EMFields), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(EMFields)             :: p3
  
    p3 = p1*real_scalar

end function sf_multiply

function fs_divide(p1, real_scalar) result (p3)
    type(EMFields), intent(in) :: p1
    real(double), intent(in)   :: real_scalar
    type(EMFields)             :: p3
  
    p3 = p1*(1.d0/real_scalar)

end function fs_divide

function lflf_add(p1, p2) result (p3)
    type(LocalEMFields), intent(in) :: p1,p2
    type(LocalEMFields)             :: p3
  
    p3%pos    = p1%pos    + p2%pos
    p3%br     = p1%br     + p2%br
    p3%bt     = p1%bt     + p2%bt
    p3%bz     = p1%bz     + p2%bz
    p3%er     = p1%er     + p2%er
    p3%et     = p1%et     + p2%et
    p3%ez     = p1%ez     + p2%ez
    p3%b_abs  = p1%b_abs  + p2%b_abs
    p3%e_abs  = p1%e_abs  + p2%e_abs
    p3%a_norm = p1%a_norm + p2%a_norm
    p3%b_norm = p1%b_norm + p2%b_norm
    p3%c_norm = p1%c_norm + p2%c_norm
    p3%e_norm = p1%e_norm + p2%e_norm

end function lflf_add

function lflf_subtract(p1, p2) result (p3)
    type(LocalEMFields), intent(in) :: p1,p2
    type(LocalEMFields)             :: p3
  
    p3%pos    = p1%pos    - p2%pos
    p3%br     = p1%br     - p2%br
    p3%bt     = p1%bt     - p2%bt
    p3%bz     = p1%bz     - p2%bz
    p3%er     = p1%er     - p2%er
    p3%et     = p1%et     - p2%et
    p3%ez     = p1%ez     - p2%ez
    p3%b_abs  = p1%b_abs  - p2%b_abs
    p3%e_abs  = p1%e_abs  - p2%e_abs
    p3%a_norm = p1%a_norm - p2%a_norm
    p3%b_norm = p1%b_norm - p2%b_norm
    p3%c_norm = p1%c_norm - p2%c_norm
    p3%e_norm = p1%e_norm - p2%e_norm  

end function lflf_subtract

function lfs_multiply(p1, real_scalar) result (p3)
    type(LocalEMFields), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalEMFields)             :: p3
  
    p3%pos  = p1%pos  * real_scalar
    p3%br   = p1%br   * real_scalar 
    p3%bt   = p1%bt   * real_scalar 
    p3%bz   = p1%bz   * real_scalar 
    p3%er   = p1%er   * real_scalar 
    p3%et   = p1%et   * real_scalar 
    p3%ez   = p1%ez   * real_scalar 
    p3%b_abs  = p1%b_abs  * real_scalar
    p3%e_abs  = p1%e_abs  * real_scalar
    p3%a_norm = p1%a_norm * real_scalar
    p3%b_norm = p1%b_norm * real_scalar
    p3%c_norm = p1%c_norm * real_scalar
    p3%e_norm = p1%e_norm * real_scalar  

end function lfs_multiply

function slf_multiply(real_scalar, p1) result (p3)
    type(LocalEMFields), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalEMFields)             :: p3
  
    p3 = p1*real_scalar

end function slf_multiply

function lfs_divide(p1, real_scalar) result (p3)
    type(LocalEMFields), intent(in) :: p1
    real(double), intent(in)        :: real_scalar
    type(LocalEMFields)             :: p3
  
    p3 = p1*(1.d0/real_scalar)

end function lfs_divide

!============================================================================
!-------------------------------I/O Routines---------------------------------
!============================================================================
subroutine read_inputs
    character(120) :: runid,version,result_dir, tables_file
    character(120) :: distribution_file, equilibrium_file
    character(120) :: geometry_file, neutrals_file
    integer        :: calc_brems,calc_bes,calc_fida,calc_npa
    integer        :: calc_birth,calc_fida_wght,calc_npa_wght
    integer        :: load_neutrals,verbose,dump_dcx
    integer(long)  :: shot,n_fida,n_npa,n_nbi,n_halo,n_dcx,n_birth
    integer(long)  :: nlambda,ne_wght,np_wght,nphi_wght,nlambda_wght
    real(double)   :: time,lambdamin,lambdamax,emax_wght
    real(double)   :: lambdamin_wght,lambdamax_wght
    real(double)   :: ai,ab,pinj,einj,species_mix(3)
    integer(long)  :: impurity_charge
    integer(long)  :: nx,ny,nz
    real(double)   :: xmin,xmax,ymin,ymax,zmin,zmax
    real(double)   :: alpha,beta,gamma,origin(3)
    logical        :: exis, error
  
    NAMELIST /fidasim_inputs/ result_dir, tables_file, &
        distribution_file, geometry_file, equilibrium_file, neutrals_file, &
        shot, time, runid, version, &
        calc_brems, calc_bes, calc_fida, calc_npa, calc_birth, &
        calc_fida_wght, calc_npa_wght, load_neutrals, dump_dcx, verbose, &
        n_fida,n_npa, n_nbi, n_halo, n_dcx, n_birth, &
        ab, pinj, einj, species_mix, ai, impurity_charge, &
        nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, &
        origin, alpha, beta, gamma, &
        ne_wght, np_wght, nphi_wght, &
        nlambda, lambdamin,lambdamax,emax_wght, &
        nlambda_wght,lambdamin_wght,lambdamax_wght
  
    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'READ_INPUTS: Input file does not exist: ', trim(namelist_file)
        stop
    endif
  
    open(13,file=namelist_file)
    read(13,NML=fidasim_inputs)
    close(13)
  
    !!General Information
    inputs%shot_number=shot
    inputs%time=time
    inputs%runid=runid
    inputs%result_dir=result_dir
  
    !!Input Files
    inputs%tables_file=tables_file
    inputs%geometry_file=geometry_file
    inputs%equilibrium_file=equilibrium_file
    inputs%distribution_file=distribution_file
    inputs%neutrals_file=neutrals_file
  
    !!Simulation Switches
    if((calc_brems+calc_bes+calc_fida).gt.0) then
        inputs%calc_spec=1
    else
        inputs%calc_spec=0
    endif
    inputs%calc_brems=calc_brems
    inputs%calc_bes=calc_bes
    inputs%calc_fida=calc_fida
    inputs%calc_npa=calc_npa
    inputs%calc_birth=calc_birth
    inputs%calc_fida_wght=calc_fida_wght
    inputs%calc_npa_wght=calc_npa_wght
    inputs%load_neutrals=load_neutrals
    inputs%dump_dcx=dump_dcx
    inputs%verbose=verbose
  
    !!Monte Carlo Settings
    inputs%n_fida=max(10,n_fida)
    inputs%n_npa=max(10,n_npa)
    inputs%n_nbi=max(10,n_nbi)
    inputs%n_halo=max(10,n_halo)
    inputs%n_dcx=max(10,n_dcx)
    inputs%n_birth= max(1,nint(n_birth/real(n_nbi)))
  
    !!Plasma Settings
    inputs%ai=ai
    inputs%impurity_charge=impurity_charge
  
    !!Neutral Beam Settings
    inputs%ab=ab
    nbi%species_mix=species_mix
    nbi%einj=einj
    nbi%pinj=pinj
  
    !!Weight Function Settings
    inputs%ne_wght=ne_wght
    inputs%np_wght=np_wght
    inputs%nphi_wght=nphi_wght
    inputs%emax_wght=emax_wght
    inputs%nlambda_wght = nlambda_wght
    inputs%lambdamin_wght=lambdamin_wght
    inputs%lambdamax_wght=lambdamax_wght
  
    !!Wavelength Grid Settings
    inputs%nlambda=nlambda
    inputs%lambdamin=lambdamin*10. !convert to angstroms
    inputs%lambdamax=lambdamax*10. !convert to angstroms
    inputs%dlambda=(inputs%lambdamax-inputs%lambdamin)/inputs%nlambda
  
    !!Beam Grid Settings
    beam_grid%nx=nx
    beam_grid%ny=ny
    beam_grid%nz=nz
    beam_grid%xmin=xmin
    beam_grid%xmax=xmax
    beam_grid%ymin=ymin
    beam_grid%ymax=ymax
    beam_grid%zmin=zmin
    beam_grid%zmax=zmax
    beam_grid%alpha=alpha
    beam_grid%beta=beta
    beam_grid%gamma=gamma
    beam_grid%origin=origin
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Shot settings ----"
        write(*,'(T2,"Shot: ",i8)') inputs%shot_number
        write(*,'(T2,"Time: ",i4," [ms]")') int(inputs%time*1.d3)
        write(*,'(T2,"Runid: ",a)') trim(adjustl(inputs%runid))
        write(*,*) ''
        write(*,'(a)') "---- input files ----"
    endif
  
    inquire(file=inputs%tables_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Tables file: ",a)') trim(inputs%tables_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Tables file does not exist: ', &
                         trim(inputs%tables_file)
        error = .True.
    endif
  
    inquire(file=inputs%geometry_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Geometry file: ",a)') trim(inputs%geometry_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Geometry file does not exist: ', &
                         trim(inputs%geometry_file)
        error = .True.
    endif
  
    inquire(file=inputs%equilibrium_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Equilibrium file: ",a)') trim(inputs%equilibrium_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Equilibrium file does not exist: ', &
                         trim(inputs%equilibrium_file)
      error = .True.
    endif
  
    inquire(file=inputs%distribution_file,exist=exis)
    if(exis) then
        if(inputs%verbose.ge.1) then
            write(*,'(T2,"Distribution file: ",a)') trim(inputs%distribution_file)
        endif
    else
        write(*,'(a,a)') 'READ_INPUTS: Distribution file does not exist: ', &
                         trim(inputs%distribution_file)
        error = .True.
    endif
    if(inputs%verbose.ge.1) then
        write(*,*) ''
    endif
  
    if(error) then 
        stop
    endif

end subroutine read_inputs

subroutine make_beam_grid
    integer(long) :: i
    real(double) :: dx, dy, dz
  
    allocate(beam_grid%xc(beam_grid%nx),  &
             beam_grid%yc(beam_grid%ny),  &
             beam_grid%zc(beam_grid%nz))
  
    dx = (beam_grid%xmax - beam_grid%xmin)/beam_grid%nx
    dy = (beam_grid%ymax - beam_grid%ymin)/beam_grid%ny
    dz = (beam_grid%zmax - beam_grid%zmin)/beam_grid%nz
  
    do i=1, beam_grid%nx
        beam_grid%xc(i) = beam_grid%xmin + (i-0.5)*dx
    enddo    
    do i=1, beam_grid%ny
        beam_grid%yc(i) = beam_grid%ymin + (i-0.5)*dy
    enddo    
    do i=1, beam_grid%nz
        beam_grid%zc(i) = beam_grid%zmin + (i-0.5)*dz
    enddo    
  
    beam_grid%dr(1) = abs(beam_grid%xc(2)-beam_grid%xc(1))
    beam_grid%dr(2) = abs(beam_grid%yc(2)-beam_grid%yc(1))
    beam_grid%dr(3) = abs(beam_grid%zc(2)-beam_grid%zc(1))
  
    beam_grid%lwh(1) = abs(beam_grid%xc(beam_grid%nx) - beam_grid%xc(1)) + beam_grid%dr(1)
    beam_grid%lwh(2) = abs(beam_grid%yc(beam_grid%ny) - beam_grid%yc(1)) + beam_grid%dr(2)
    beam_grid%lwh(3) = abs(beam_grid%zc(beam_grid%nz) - beam_grid%zc(1)) + beam_grid%dr(3)
  
    beam_grid%volume = beam_grid%lwh(1)*beam_grid%lwh(3)*beam_grid%lwh(3)
  
    beam_grid%center(1) = (minval(beam_grid%xc) - 0.5*beam_grid%dr(1)) + 0.5*beam_grid%lwh(1)
    beam_grid%center(2) = (minval(beam_grid%yc) - 0.5*beam_grid%dr(2)) + 0.5*beam_grid%lwh(2)
    beam_grid%center(3) = (minval(beam_grid%zc) - 0.5*beam_grid%dr(3)) + 0.5*beam_grid%lwh(3)
  
    beam_grid%drmin  = minval(beam_grid%dr)
    beam_grid%dv     = beam_grid%dr(1)*beam_grid%dr(2)*beam_grid%dr(3)
    beam_grid%ntrack = beam_grid%nx+beam_grid%ny+beam_grid%nz
    beam_grid%ngrid  = beam_grid%nx*beam_grid%ny*beam_grid%nz
  
    call tb_zyx(beam_grid%alpha,beam_grid%beta,beam_grid%gamma, &
                beam_grid%basis, beam_grid%inv_basis)
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Beam grid settings ----"
        write(*,'(T2,"Nx: ", i3)') beam_grid%nx
        write(*,'(T2,"Ny: ", i3)') beam_grid%ny
        write(*,'(T2,"Nz: ", i3)') beam_grid%nz
        write(*,'(T2,"dV: ", f5.2," [cm^3]")') beam_grid%dv
        write(*,'(T2,"alpha: ",f5.2," [rad]")') beam_grid%alpha
        write(*,'(T2,"beta:  ",f5.2," [rad]")') beam_grid%beta
        write(*,'(T2,"gamma: ",f5.2," [rad]")') beam_grid%gamma
        write(*,'(T2,"origin: [",f7.2,",",f7.2,",",f7.2,"] [cm]")') beam_grid%origin
        write(*,*) ''
    endif

end subroutine make_beam_grid

subroutine read_beam
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(1) :: dims
  
    real(double), dimension(3) ::uvw_src,uvw_axis,pos
    real(double) :: dis
  
    integer :: error
  
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDONLY_F, fid, error)
  
    !!Open NBI group
    call h5gopen_f(fid, "/nbi", gid, error)
  
    !!Read in datasets
    call h5ltread_dataset_string_f(gid, "/nbi/name",nbi%name, error)
    dims(1) = 3
    call h5ltread_dataset_double_f(gid, "/nbi/src", uvw_src, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/axis", uvw_axis, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/divy", nbi%divy, dims, error)
    call h5ltread_dataset_double_f(gid, "/nbi/divz", nbi%divz, dims, error)
    call h5ltread_dataset_int_scalar_f(gid, "/nbi/shape", nbi%shape, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/focy", nbi%focy, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/focz", nbi%focz, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/widy", nbi%widy, error)
    call h5ltread_dataset_double_scalar_f(gid, "/nbi/widz", nbi%widz, error)
  
    !!Close NBI group
    call h5gclose_f(gid, error)
  
    !!Close file id
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
    !!Convert to beam grid coordinates
    call uvw_to_xyz(uvw_src,nbi%src)
    nbi%axis = matmul(beam_grid%inv_basis,uvw_axis)
  
    nbi%vinj=sqrt(2.d0*nbi%einj*1.d3 *e0/(inputs%ab*mass_u))*1.d2 !! [cm/s]
    pos = nbi%src + 200.0*nbi%axis
    dis = sqrt(sum((pos - nbi%src)**2.0))
    nbi%beta = asin((nbi%src(3) - pos(3))/dis)
    nbi%alpha = atan2(pos(2)-nbi%src(2),pos(1)-nbi%src(1))
  
    call tb_zyx(nbi%alpha,nbi%beta,0.d0,nbi%basis,nbi%inv_basis)
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Neutral beam settings ----'
        write(*,'(T2,"Beam: ",a)') nbi%name
        write(*,'(T2,"Power:   ",f5.2," [MW]")') nbi%pinj
        write(*,'(T2,"Voltage: ",f5.2," [keV]")') nbi%einj
        write(*,*) ''
    endif

end subroutine read_beam

subroutine read_chords
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
    logical :: path_valid
  
    real(double), dimension(:,:), allocatable :: lenses
    real(double), dimension(:,:), allocatable :: axes
    real(double), dimension(:), allocatable :: spot_size, sigma_pi
    real(double) :: r0(3), v0(3), r_enter(3), r_exit(3)
    real(double) :: xyz_lens(3), xyz_axis(3), length
    real(double), dimension(3,3) :: basis
    real(double), dimension(2) :: randomu
    real(double) :: theta, sqrt_rho
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    character(len=20) :: system = ''
  
    integer :: i, j, ic, nc, ncell, ind(3)
    integer :: error
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- FIDA/BES settings ----'
    endif
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDONLY_F, fid, error)
  
    !!Check if SPEC group exists
    call h5ltpath_valid_f(fid, "/spec", .True., path_valid, error)
    if(.not.path_valid) then
        write(*,'(a)') 'FIDA/BES geometry is not in the geometry file'
        write(*,'(a)') 'Continuing without spectral diagnostics' 
        inputs%calc_spec = 0
        inputs%calc_fida = 0
        inputs%calc_bes = 0
        inputs%calc_brems = 0
        inputs%calc_fida_wght = 0
        return
    endif
  
    !!Open SPEC group
    call h5gopen_f(fid, "/spec", gid, error)
  
    call h5ltread_dataset_string_f(gid, "/spec/system", system, error)
    call h5ltread_dataset_int_scalar_f(gid, "/spec/nchan", spec_chords%nchan, error)
  
    allocate(lenses(3, spec_chords%nchan))
    allocate(axes(3, spec_chords%nchan))
    allocate(spot_size(spec_chords%nchan))
    allocate(sigma_pi(spec_chords%nchan))
    allocate(spec_chords%los(spec_chords%nchan))
    allocate(spec_chords%radius(spec_chords%nchan))
    allocate(spec_chords%dlength(spec_chords%nchan, &
                                 beam_grid%nx, & 
                                 beam_grid%ny, &
                                 beam_grid%nz) )
  
    allocate(spec_chords%los_inter(beam_grid%nx, & 
                                 beam_grid%ny, &
                                 beam_grid%nz) )
  
    spec_chords%dlength = 0.d0
    spec_chords%los_inter = .False.
  
    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid, "/spec/lens", lenses, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/axis", axes, dims, error)
    call h5ltread_dataset_double_f(gid, "/spec/spot_size", spot_size, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/sigma_pi", sigma_pi, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/spec/radius", spec_chords%radius, dims(2:2), error)
  
    !!Close SPEC group
    call h5gclose_f(gid, error)
  
    !!Close file id
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
    chan_loop: do i=1,spec_chords%nchan
        call uvw_to_xyz(lenses(:,i),xyz_lens)
        xyz_axis = matmul(beam_grid%inv_basis,axes(:,i))
        spec_chords%los(i)%lens = xyz_lens
        spec_chords%los(i)%axis = xyz_axis
        spec_chords%los(i)%sigma_pi = sigma_pi(i)
        spec_chords%los(i)%spot_size = spot_size(i)
  
        r0 = xyz_lens
        v0 = xyz_axis
        v0 = v0/normp(v0)
        call line_basis(r0,v0,basis)
  
        call grid_intersect(r0,v0,length,r_enter,r_exit)
        if(length.le.0.d0) then
            WRITE(*,'("Channel ",i3," missed the beam grid")'),i
            cycle chan_loop 
        endif
  
        if(spot_size(i).le.0.d0) then 
            nc = 1
        else
            nc = 100
        endif
  
        !$OMP PARALLEL DO schedule(guided) private(ic,randomu,sqrt_rho,theta,r0, &
        !$OMP& length, r_enter, r_exit, j, tracks, ncell, ind)
        do ic=1,nc
            ! Uniformally sample within spot size
            call randu(randomu)
            sqrt_rho = sqrt(randomu(1))
            theta = 2*pi*randomu(2)
            r0(1) = 0.d0
            r0(2) = spot_size(i)*sqrt_rho*cos(theta)
            r0(3) = spot_size(i)*sqrt_rho*sin(theta)
            r0 = matmul(basis,r0) + xyz_lens
  
            call grid_intersect(r0, v0, length, r_enter, r_exit)
            call track(r_enter, v0, tracks, ncell)
            track_loop: do j=1, ncell
                ind = tracks(j)%ind
                !inds can repeat so add rather than assign
                !$OMP CRITICAL(read_chords_1)
                spec_chords%dlength(i,ind(1),ind(2),ind(3)) = & 
                spec_chords%dlength(i,ind(1),ind(2),ind(3)) + tracks(j)%time/real(nc) !time == distance
                spec_chords%los_inter(ind(1),ind(2),ind(3)) = .True.
                !$OMP END CRITICAL(read_chords_1)
            enddo track_loop
        enddo
        !$OMP END PARALLEL DO
    enddo chan_loop
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"FIDA/BES System: ",a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i3)') spec_chords%nchan
        write(*,*) ''
    endif
  
    deallocate(lenses,axes,spot_size,sigma_pi)

end subroutine read_chords

subroutine read_npa
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
    logical :: path_valid
  
    real(double), dimension(:,:), allocatable :: a_tedge,a_redge,a_cent
    real(double), dimension(:,:), allocatable :: d_tedge,d_redge,d_cent
    integer, dimension(:), allocatable :: a_shape, d_shape
    character(len=20) :: system = ''
    
    real(double), dimension(3) :: xyz_a_tedge,xyz_a_redge,xyz_a_cent
    real(double), dimension(3) :: xyz_d_tedge,xyz_d_redge,xyz_d_cent
    real(double), dimension(3) :: eff_rd, rd, rd_d, r0, r0_d, v0
    real(double), dimension(3,3) :: basis, inv_basis
    real(double), dimension(50) :: xd, yd
    real(double), dimension(:,:,:,:,:), allocatable :: effrd
    real(double), dimension(:,:,:,:), allocatable :: phit
    real(double) :: total_prob, hh, hw, dprob, dx, dy, r
    integer :: ichan,i,j,k,ix,iy,d_index,nd,cnt
    integer :: error
    
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%geometry_file, H5F_ACC_RDWR_F, fid, error)
  
    !!Check if NPA group exists
    call h5ltpath_valid_f(fid, "/npa", .True., path_valid, error)
    if(.not.path_valid) then
        write(*,'(a)') 'NPA geometry is not in the geometry file'
        write(*,'(a)') 'Continuing without NPA diagnostics' 
        inputs%calc_npa = 0
        inputs%calc_npa_wght = 0
        return
    endif
  
    !!Open NPA group
    call h5gopen_f(fid, "/npa", gid, error)
  
    call h5ltread_dataset_string_f(gid, "/npa/system", system, error)
    call h5ltread_dataset_int_scalar_f(gid, "/npa/nchan", npa_chords%nchan, error)
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- NPA settings ----"
        write(*,'(T2,"NPA System: ", a)') trim(adjustl(system))
        write(*,'(T2,"Number of channels: ",i3)') npa_chords%nchan
    endif
  
    allocate(a_tedge(3, npa_chords%nchan))
    allocate(a_redge(3, npa_chords%nchan))
    allocate(a_cent(3,  npa_chords%nchan))
    allocate(a_shape(npa_chords%nchan))    
    allocate(d_tedge(3, npa_chords%nchan))
    allocate(d_redge(3, npa_chords%nchan))
    allocate(d_cent(3,  npa_chords%nchan))
    allocate(d_shape(npa_chords%nchan))
    allocate(npa_chords%radius(npa_chords%nchan))
    allocate(npa_chords%det(npa_chords%nchan))
    allocate(npa_chords%phit(beam_grid%nx, & 
                             beam_grid%ny, &
                             beam_grid%nz, &
                             npa_chords%nchan) )
    allocate(npa_chords%hit(beam_grid%nx, & 
                            beam_grid%ny, &
                            beam_grid%nz) )
    npa_chords%hit = .False.
  
    allocate(effrd(3,beam_grid%nx,&
                     beam_grid%ny,&
                     beam_grid%nz,&
                     npa_chords%nchan) )
    allocate(phit(beam_grid%nx,&
                  beam_grid%ny,&
                  beam_grid%nz,&
                  npa_chords%nchan) )
    effrd = 0.d0
    phit = 0.d0
  
    dims = [3,spec_chords%nchan]
    call h5ltread_dataset_double_f(gid,"/npa/radius", npa_chords%radius, dims(2:2), error)
    call h5ltread_dataset_int_f(gid, "/npa/a_shape", a_shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/a_tedge", a_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_redge", a_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/a_cent",  a_cent, dims, error)
  
    call h5ltread_dataset_int_f(gid, "/npa/d_shape", d_shape, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/npa/d_tedge", d_tedge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_redge", d_redge, dims, error)
    call h5ltread_dataset_double_f(gid, "/npa/d_cent",  d_cent, dims, error)
  
    !!Close NPA group
    call h5gclose_f(gid, error)
  
    !!Close file id
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
    ! Define detector/aperture shape
    npa_chords%det%detector%shape = d_shape
    npa_chords%det%aperture%shape = a_shape
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,a)') "Calculating hit probabilities for NPA channels"
    endif
    chan_loop: do ichan=1,npa_chords%nchan
        ! Convert to beam grid coordinates
        call uvw_to_xyz(a_cent(:,ichan), xyz_a_cent)
        call uvw_to_xyz(a_redge(:,ichan),xyz_a_redge)
        call uvw_to_xyz(a_tedge(:,ichan),xyz_a_tedge)
        call uvw_to_xyz(d_cent(:,ichan), xyz_d_cent)
        call uvw_to_xyz(d_redge(:,ichan),xyz_d_redge)
        call uvw_to_xyz(d_tedge(:,ichan),xyz_d_tedge)
  
        ! Define detector/aperture hh/hw
        npa_chords%det(ichan)%detector%hw = normp(xyz_d_redge - xyz_d_cent)
        npa_chords%det(ichan)%aperture%hw = normp(xyz_a_redge - xyz_a_cent)
  
        npa_chords%det(ichan)%detector%hh = normp(xyz_d_tedge - xyz_d_cent)
        npa_chords%det(ichan)%aperture%hh = normp(xyz_a_tedge - xyz_a_cent)
  
        ! Define detector/aperture origin
        npa_chords%det(ichan)%detector%origin = xyz_d_cent
        npa_chords%det(ichan)%aperture%origin = xyz_a_cent
  
        ! Define detector/aperture basis
        call plane_basis(xyz_d_cent, xyz_d_redge, xyz_d_tedge, &
             npa_chords%det(ichan)%detector%basis, &
             npa_chords%det(ichan)%detector%inv_basis)
        call plane_basis(xyz_a_cent, xyz_a_redge, xyz_a_tedge, &
             npa_chords%det(ichan)%aperture%basis, &
             npa_chords%det(ichan)%aperture%inv_basis)
   
        hw = npa_chords%det(ichan)%detector%hw
        hh = npa_chords%det(ichan)%detector%hh
        nd = size(xd)
        do i=1,nd
            xd(i) = -hw + 2*hw*(i-0.5)/real(nd)
            yd(i) = -hh + 2*hh*(i-0.5)/real(nd)
        enddo
        dx = abs(xd(2) - xd(1))
        dy = abs(yd(2) - yd(1))
        basis = npa_chords%det(ichan)%detector%basis
        inv_basis = npa_chords%det(ichan)%detector%inv_basis
        cnt = 0
        ! For each grid point find the probability of hitting the detector given an isotropic source
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,ix,iy,total_prob,eff_rd,r0,r0_d, &
        !$OMP& rd_d,rd,d_index,v0,dprob,r)
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    cnt = cnt+1
                    total_prob = 0.d0
                    eff_rd = eff_rd*0.d0
                    r0 = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)]
                    r0_d = matmul(inv_basis,r0-xyz_d_cent)
                    do ix = 1, nd
                        do iy = 1, nd
                            rd_d = [xd(ix),yd(iy),0.d0]
                            rd = matmul(basis,rd_d) + xyz_d_cent
                            v0 = rd - r0
                            d_index = 0
                            call hit_npa_detector(r0,v0,d_index)
                            if(d_index.ne.0) then
                                r = normp(rd_d - r0_d)**2.0
                                dprob = (dx*dy) * ((4*pi)**(-1.0)) * r0_d(3) * (r**(-1.5))
                                eff_rd = eff_rd + dprob*rd
                                total_prob = total_prob + dprob
                            endif
                        enddo !yd loop
                    enddo !xd loop
                    if(total_prob.gt.0.0) then
                        eff_rd = eff_rd/total_prob
                        phit(i,j,k,ichan) = total_prob
                        call xyz_to_uvw(eff_rd,effrd(:,i,j,k,ichan))
                        npa_chords%phit(i,j,k,ichan)%p = total_prob
                        npa_chords%phit(i,j,k,ichan)%eff_rd = eff_rd
                        npa_chords%hit(i,j,k) = .True.
                    endif
                    if(inputs%verbose.ge.2) then
                        WRITE(*,'(T4,"Channel: ",i3," ",f7.2,"% completed",a,$)') &
                                 ichan, cnt/real(beam_grid%ngrid)*100,char(13)
                    endif
                enddo !x loop
            enddo !y loop
        enddo !z loop
        !$OMP END PARALLEL DO
  
        total_prob = sum(npa_chords%phit(:,:,:,ichan)%p)
        if(total_prob.le.0.d0) then
            WRITE(*,'("Channel ",i3," missed the beam grid")'),ichan
            cycle chan_loop 
        endif
  
    enddo chan_loop
    if(inputs%verbose.ge.1) write(*,'(50X,a)') ""
  
    deallocate(phit,effrd)
    deallocate(a_shape,a_cent,a_redge,a_tedge)
    deallocate(d_shape,d_cent,d_redge,d_tedge)

end subroutine read_npa

subroutine read_equilibrium
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dims
  
    integer :: impc
    integer :: error
    
    integer, dimension(:,:), allocatable :: p_mask,f_mask
  
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%equilibrium_file, H5F_ACC_RDONLY_F, fid, error)
  
    !!Open PLASMA group
    call h5gopen_f(fid, "/plasma", gid, error)
  
    !!Read in interpolation grid
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nr", inter_grid%nr, error)
    call h5ltread_dataset_int_scalar_f(gid, "/plasma/nz", inter_grid%nz, error)
  
    allocate(inter_grid%r(inter_grid%nr),inter_grid%z(inter_grid%nz))
    allocate(inter_grid%r2d(inter_grid%nr,inter_grid%nz))
    allocate(inter_grid%z2d(inter_grid%nr,inter_grid%nz))
    allocate(p_mask(inter_grid%nr,inter_grid%nz))
    allocate(f_mask(inter_grid%nr,inter_grid%nz))
  
    dims = [inter_grid%nr, inter_grid%nz]
    call h5ltread_dataset_double_f(gid, "/plasma/r", inter_grid%r, dims(1:1), error)
    call h5ltread_dataset_double_f(gid, "/plasma/z", inter_grid%z, dims(2:2), error)
    call h5ltread_dataset_double_f(gid, "/plasma/r2d", inter_grid%r2d, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/z2d", inter_grid%z2d, dims, error)
  
    inter_grid%dr = abs(inter_grid%r(2)-inter_grid%r(1))
    inter_grid%dz = abs(inter_grid%z(2)-inter_grid%z(1))
    inter_grid%da = inter_grid%dr*inter_grid%dz
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Interpolation grid settings ----'
        write(*,'(T2,"Nr: ",i3)') inter_grid%nr
        write(*,'(T2,"Nz: ",i3)') inter_grid%nz
        write(*,'(T2,"dA: ", f4.2," [cm^2]")') inter_grid%da
        write(*,*) ''
    endif
  
    !!Read in plasma parameters
    allocate(equil%plasma(inter_grid%nr,inter_grid%nz))
  
    call h5ltread_dataset_double_f(gid, "/plasma/dene", equil%plasma%dene, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/te", equil%plasma%te, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/ti", equil%plasma%ti, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/zeff", equil%plasma%zeff, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vr", equil%plasma%vr, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vt", equil%plasma%vt, dims, error)
    call h5ltread_dataset_double_f(gid, "/plasma/vz", equil%plasma%vz, dims, error)
    call h5ltread_dataset_int_f(gid, "/plasma/mask", p_mask, dims,error)
  
    impc = inputs%impurity_charge
    equil%plasma%denimp = ((equil%plasma%zeff-1.d0)/(impc*(impc-1.d0)))*equil%plasma%dene
    equil%plasma%denp = equil%plasma%dene - impc*equil%plasma%denimp
  
    !!Close PLASMA group
    call h5gclose_f(gid, error)
  
    !!Open FIELDS group
    call h5gopen_f(fid, "/fields", gid, error)
    
    allocate(equil%fields(inter_grid%nr,inter_grid%nz))
  
    !!Read in electromagnetic fields
    call h5ltread_dataset_double_f(gid, "/fields/br", equil%fields%br, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bt", equil%fields%bt, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/bz", equil%fields%bz, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/er", equil%fields%er, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/et", equil%fields%et, dims, error)
    call h5ltread_dataset_double_f(gid, "/fields/ez", equil%fields%ez, dims, error)
    call h5ltread_dataset_int_f(gid, "/fields/mask", f_mask, dims,error)
  
    !!Close FIELDS group
    call h5gclose_f(gid, error)
  
    !!Close file
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
    allocate(equil%mask(inter_grid%nr,inter_grid%nz))
    equil%mask = 0.d0
    where ((p_mask.eq.1).and.(f_mask.eq.1)) equil%mask = 1.d0
  
end subroutine read_equilibrium

subroutine read_f(fid, error)
    integer(HID_T), intent(inout) :: fid
    integer, intent(out)          :: error
  
    integer(HSIZE_T), dimension(4) :: dims
    real(double) :: dummy(1)
  
    if(inputs%verbose.ge.1) then 
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif
  
    call h5ltread_dataset_int_scalar_f(fid,"/nenergy", fbm%nenergy, error)
    call h5ltread_dataset_int_scalar_f(fid,"/npitch", fbm%npitch, error)
  
    allocate(fbm%energy(fbm%nenergy), fbm%pitch(fbm%npitch))
    allocate(fbm%f(fbm%nenergy, fbm%npitch, inter_grid%nr, inter_grid%nz))
  
    dims = [fbm%nenergy, fbm%npitch, inter_grid%nr, inter_grid%nz]
    call h5ltread_dataset_double_f(fid, "/energy", fbm%energy, dims(1:1), error)
    call h5ltread_dataset_double_f(fid, "/pitch", fbm%pitch, dims(2:2), error)
    call h5ltread_dataset_double_f(fid, "/f", fbm%f, dims, error)
    call h5ltread_dataset_double_f(fid, "/denf",equil%plasma%denf, dims(3:4), error)
  
    fbm%dE = abs(fbm%energy(2) - fbm%energy(1))
    fbm%dp = abs(fbm%pitch(2) - fbm%pitch(1))
  
    dummy = minval(fbm%energy)
    fbm%emin = dummy(1)
    dummy = maxval(fbm%energy)
    fbm%emax = dummy(1)
    fbm%e_range = fbm%emax - fbm%emin
    dummy = minval(fbm%pitch)
    fbm%pmin = dummy(1)
    dummy = maxval(fbm%pitch)
    fbm%pmax = dummy(1)
    fbm%p_range = fbm%pmax - fbm%pmin
  
    if(inputs%verbose.ge.1) then 
        write(*,'(T2,"Distribution type: ",a)') "Fast-ion Density Function F(energy,pitch,R,Z)"
        write(*,'(T2,"Nenergy = ",i3)'),fbm%nenergy
        write(*,'(T2,"Npitch  = ",i3)'),fbm%npitch 
        write(*,'(T2,"Energy range = [",f5.2,",",f6.2,"]")'),fbm%emin,fbm%emax
        write(*,'(T2,"Pitch  range = [",f5.2,",",f5.2,"]")'),fbm%pmin,fbm%pmax
        write(*,*) ''
    endif

end subroutine read_f

subroutine read_mc(fid, error)
    integer(HID_T), intent(inout) :: fid
    integer, intent(out)          :: error
  
    integer(HSIZE_T), dimension(1) :: dims
  
    integer(long)         :: i,ii,ir,iz,nphi
    real(double)    :: phi,xmin,xmax,ymin,ymax,zmin,zmax
    real(double)    :: phi_enter,delta_phi
    real(double), dimension(3) :: uvw,xyz
    integer(long), dimension(1) :: minpos
    logical :: in_grid
    real(double), dimension(:), allocatable :: weight
    real(double), dimension(:), allocatable :: r, z, vr, vt, vz
    real(double), dimension(:), allocatable :: energy, pitch
    integer(long), dimension(:), allocatable :: orbit_class
    integer :: cnt,num
    character(len=32) :: dist_type_name = ''
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') '---- Fast-ion distribution settings ----'
    endif
  
    call h5ltread_dataset_int_scalar_f(fid, "/nparticle", particles%nparticle, error)
    call h5ltread_dataset_int_scalar_f(fid, "/nclass", particles%nclass, error)
  
    !!ALLOCATE SPACE
    allocate(particles%fast_ion(particles%nparticle))
    allocate(r(particles%nparticle))
    allocate(z(particles%nparticle))
    allocate(weight(particles%nparticle))
    allocate(orbit_class(particles%nparticle))
  
    dims(1) = particles%nparticle
    call h5ltread_dataset_double_f(fid, "/r", r, dims, error)
    call h5ltread_dataset_double_f(fid, "/z", z, dims, error)
    call h5ltread_dataset_double_f(fid, "/weight", weight, dims, error)
    call h5ltread_dataset_int_f(fid, "/class", orbit_class, dims, error)
  
    if(inputs%dist_type.eq.2) then
        dist_type_name = "Guiding Center Monte Carlo"
        allocate(energy(particles%nparticle))
        allocate(pitch(particles%nparticle))
        call h5ltread_dataset_double_f(fid, "/energy", energy, dims, error)
        call h5ltread_dataset_double_f(fid, "/pitch", pitch, dims, error)
        particles%guiding_center = .True.
    else
        dist_type_name = "Full Orbit Monte Carlo"
        allocate(vr(particles%nparticle))
        allocate(vt(particles%nparticle))
        allocate(vz(particles%nparticle))
        call h5ltread_dataset_double_f(fid, "/vr", vr, dims, error)
        call h5ltread_dataset_double_f(fid, "/vt", vt, dims, error)
        call h5ltread_dataset_double_f(fid, "/vz", vz, dims, error)
        particles%guiding_center = .False.
    endif
  
    xmin = beam_grid%xmin
    xmax = beam_grid%xmax
    ymin = beam_grid%ymin
    ymax = beam_grid%ymax
    zmin = beam_grid%zmin
    zmax = beam_grid%zmax
  
    cnt=0
    num=0
    nphi = 1000
    !$OMP PARALLEL DO schedule(guided) private(i,ii,ir,iz,phi,uvw,xyz, &
    !$OMP& in_grid,minpos,delta_phi,phi_enter) 
    do i=1,particles%nparticle
        if(inputs%verbose.ge.2) then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt/real(particles%nparticle)*100,char(13)
        endif
        delta_phi = 0.0
        phi_enter = 0.0
        phi = 0.0
        do ii=1,nphi
            phi = phi + 2*pi/dble(nphi)
            uvw(1) = r(i)*cos(phi)
            uvw(2) = r(i)*sin(phi)
            uvw(3) = z(i)
            call uvw_to_xyz(uvw,xyz)
            
            if(((xyz(1).ge.xmin.and.xyz(1).le.xmax).and. &
                (xyz(2).ge.ymin.and.xyz(2).le.ymax).and. &
                (xyz(3).ge.zmin.and.xyz(3).le.zmax))) then
                delta_phi = delta_phi + 2*pi/dble(nphi)
                if(ii.ne.1)then
                  if(.not.in_grid) phi_enter = phi
                endif
                in_grid = .True.
            else
                in_grid = .False.
            endif
        enddo
  
        minpos = minloc(abs(inter_grid%r - r(i)))
        ir = minpos(1)
        minpos = minloc(abs(inter_grid%z - z(i)))
        iz = minpos(1)
        equil%plasma(ir,iz)%denf = equil%plasma(ir,iz)%denf + weight(i)/(2*pi*r(i)*inter_grid%da)
  
        particles%fast_ion(i)%r = r(i)
        particles%fast_ion(i)%z = z(i)
        particles%fast_ion(i)%phi_enter = phi_enter
        particles%fast_ion(i)%delta_phi = delta_phi
        particles%fast_ion(i)%weight = weight(i)*(delta_phi/(2*pi))/beam_grid%dv
        particles%fast_ion(i)%class = orbit_class(i)
        if(delta_phi.gt.0) then 
            particles%fast_ion(i)%cross_grid = .True.
        else
            particles%fast_ion(i)%cross_grid = .False.
        endif
        if(particles%guiding_center) then
            particles%fast_ion(i)%energy = energy(i)
            particles%fast_ion(i)%pitch = pitch(i)
            particles%fast_ion(i)%vabs = sqrt(energy(i)/(v2_to_E_per_amu*inputs%ab))
        else
            particles%fast_ion(i)%vr = vr(i)
            particles%fast_ion(i)%vt = vt(i)
            particles%fast_ion(i)%vz = vz(i)
            particles%fast_ion(i)%vabs = sqrt(vr(i)**2 + vt(i)**2 + vz(i)**2)
        endif
        if(delta_phi.gt.0.d0) num = num + 1
        cnt=cnt+1
    enddo
    !$OMP END PARALLEL DO
  
    if(num.le.0) then
        write(*,'(a)') 'READ_MC: No mc particles in beam grid'
        stop
    endif
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Distribution type: ",a)') dist_type_name
        write(*,'(T2,"Number of mc particles: ",i9)') particles%nparticle
        write(*,'(T2,"Number of orbit classes: ",i6)') particles%nclass
        write(*,*) ''
    endif
  
end subroutine read_mc

subroutine read_distribution
    integer(HID_T) :: fid
    integer :: error
    
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%distribution_file, H5F_ACC_RDONLY_F, fid, error)
  
    !!Get distribution type
    call h5ltread_dataset_int_scalar_f(fid, "/type", inputs%dist_type, error)
    
    if(inputs%dist_type.eq.1) then
        call read_f(fid, error)
    else !2 or 3
        call read_mc(fid, error)
    endif
  
    !!Close file
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
end subroutine read_distribution

subroutine read_cross(fid, grp, cross)
    integer(HID_T), intent(in)              :: fid
    character(len=*), intent(in)            :: grp
    type(AtomicCrossSection), intent(inout) :: cross
  
    integer(HSIZE_T), dimension(3) :: dim3
    real(double) :: emin, emax, rmin
    integer :: i, n_max, m_max, error
    real(double), dimension(:,:,:), allocatable :: dummy3
  
    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", cross%nenergy, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/m_max", m_max, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid,grp//"/dlogE", cross%dlogE, error)
    cross%logemin = log10(emin)
    cross%logemax = log10(emax)
  
    allocate(dummy3(n_max, m_max, cross%nenergy))
  
    allocate(cross%log_cross(cross%m_max,cross%n_max, cross%nenergy))
  
    dim3 = [n_max, m_max, tables%H_H_cx%nenergy]
    call h5ltread_dataset_double_f(fid,grp//"/cx", dummy3, dim3, error)
    rmin = minval(dummy3,dummy3.gt.0.d0)
    where (dummy3.le.0.0)
        dummy3 = 0.9*rmin
    end where
    cross%minlog_cross = log10(rmin)
    do i=1, cross%nenergy
        cross%log_cross(:,:,i) = log10(transpose(dummy3(1:nlevs,1:nlevs,i)))
    enddo
    deallocate(dummy3)

end subroutine read_cross

subroutine read_rates(fid, grp, b_amu, t_amu, rates)
    integer(HID_T), intent(in)             :: fid
    character(len=*), intent(in)           :: grp
    real(double), dimension(2), intent(in) :: b_amu
    real(double), intent(in)               :: t_amu
    type(AtomicRates), intent(inout)       :: rates
  
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
    logical :: path_valid
    integer :: i, j, n, n_max, m_max, error
    integer :: n_bt_amu, tt_ind, bt_ind, drank
    real(double) :: emin,emax,tmin,tmax,rmin
    real(double) :: bt_min, tt_min, tt_dum, bt_dum
    real(double), dimension(2) :: bt_amu, tt_amu
    real(double), dimension(:,:), allocatable :: dummy2
    real(double), dimension(:,:,:,:), allocatable :: dummy4
    real(double), dimension(:,:,:,:,:), allocatable :: dummy5
  
    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_bt_amu", n_bt_amu, error)
    allocate(dummy2(2, n_bt_amu))
    dim2 = [2, n_bt_amu]
    call h5ltread_dataset_double_f(fid, grp//"/bt_amu", dummy2, dim2, error)
  
    call h5ltread_dataset_int_scalar_f(fid, grp//"/n_max", n_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/m_max", m_max, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/nenergy", rates%nenergy, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emin", emin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/emax", emax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogE", rates%dlogE, error)
    call h5ltread_dataset_int_scalar_f(fid, grp//"/ntemp", rates%ntemp, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmin", tmin, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/tmax", tmax, error)
    call h5ltread_dataset_double_scalar_f(fid, grp//"/dlogT", rates%dlogT, error)
    rates%logemin = log10(emin)
    rates%logemax = log10(emax)
    rates%logtmin = log10(tmin)
    rates%logtmax = log10(tmax)
  
    bt_ind = 1
    tt_ind = 1
    bt_amu = [b_amu(1), t_amu]
    tt_amu = [b_amu(2), t_amu]
    bt_min = normp(bt_amu - dummy2(:,1))
    tt_min = normp(tt_amu - dummy2(:,1))
    do i=2, n_bt_amu
        bt_dum = normp(bt_amu - dummy2(:,i))
        tt_dum = normp(tt_amu - dummy2(:,i))
        if(bt_dum.lt.bt_min) then
            bt_min = bt_dum
            bt_ind = i
        endif
        if(tt_dum.lt.tt_min) then
            tt_min = tt_dum
            tt_ind = i
        endif
    enddo
    rates%ab(1) = dummy2(1,bt_ind)
    rates%ab(2) = dummy2(1,tt_ind)
  
    deallocate(dummy2)
  
    allocate(rates%log_pop(&
                    rates%m_max, &
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp, 2))
    allocate(rates%log_depop(&
                    rates%n_max, &
                    rates%nenergy, &
                    rates%ntemp, 2))
    rates%log_pop = 0.d0
    rates%log_depop = 0.d0
  
    !!Read CX
    call h5ltpath_valid_f(fid, grp//"/cx", .True., path_valid, error)
    if(path_valid) then
        call h5ltget_dataset_ndims_f(fid, grp//"/cx", drank, error)
        if(drank.eq.4) then
            allocate(dummy4(n_max, &
                           rates%nenergy, &
                           rates%ntemp, n_bt_amu))
            dim4 = [n_max, rates%nenergy, rates%ntemp,n_bt_amu]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy4, dim4, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j,1) = dummy4(n,i,j,bt_ind)
                        rates%log_depop(n,i,j,2) = dummy4(n,i,j,tt_ind)
                    enddo
                enddo
            enddo
            deallocate(dummy4)
        endif
        if(drank.eq.5) then
            allocate(dummy5(n_max, m_max, &
                           rates%nenergy, &
                           rates%ntemp, n_bt_amu))
            dim5 = [n_max, m_max, rates%nenergy, rates%ntemp,n_bt_amu]
            call h5ltread_dataset_double_f(fid, grp//"/cx", dummy5, dim5, error)
            do j=1,rates%ntemp
                do i=1,rates%nenergy
                    do n=1,rates%n_max
                        rates%log_depop(n,i,j,1) = sum(dummy5(n,:,i,j,bt_ind))
                        rates%log_depop(n,i,j,2) = sum(dummy5(n,:,i,j,tt_ind))
                    enddo
                enddo
            enddo
            deallocate(dummy5)
        endif
    endif
    
    !!Read ionization
    call h5ltpath_valid_f(fid, grp//"/ionization", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy4(n_max, &
                       rates%nenergy, &
                       rates%ntemp, n_bt_amu))
        dim4 = [n_max, rates%nenergy, rates%ntemp,n_bt_amu]
        call h5ltread_dataset_double_f(fid, grp//"/ionization", dummy4, dim4, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                do n=1,rates%n_max
                    rates%log_depop(n,i,j,1) = rates%log_depop(n,i,j,1) + &
                                               dummy4(n,i,j,bt_ind)
                    rates%log_depop(n,i,j,2) = rates%log_depop(n,i,j,2) + &
                                               dummy4(n,i,j,tt_ind)
                enddo
            enddo
        enddo
        deallocate(dummy4)
    endif
  
    !!Read excitation
    call h5ltpath_valid_f(fid, grp//"/excitation", .True., path_valid, error)
    if(path_valid) then
        allocate(dummy5(n_max, m_max,&
                       rates%nenergy, &
                       rates%ntemp, n_bt_amu))
        dim5 = [n_max, m_max, rates%nenergy, rates%ntemp,n_bt_amu]
        call h5ltread_dataset_double_f(fid, grp//"/excitation", dummy5, dim5, error)
        do j=1,rates%ntemp
            do i=1,rates%nenergy
                rates%log_pop(:,:,i,j,1) = transpose(dummy5(1:nlevs,1:nlevs,i,j,bt_ind))
                rates%log_pop(:,:,i,j,2) = transpose(dummy5(1:nlevs,1:nlevs,i,j,tt_ind))
                do n=1,rates%n_max
                    rates%log_depop(n,i,j,1) = rates%log_depop(n,i,j,1) + &
                                               sum(dummy5(n,:,i,j,bt_ind))
                    rates%log_depop(n,i,j,2) = rates%log_depop(n,i,j,2) + &
                                               sum(dummy5(n,:,i,j,tt_ind))
                enddo
            enddo
        enddo
        deallocate(dummy5)
    endif
  
    rmin = minval(rates%log_depop, rates%log_depop.gt.0.d0)
    where (rates%log_depop.le.0.d0)
        rates%log_depop = 0.9*rmin
    end where
    rates%minlog_depop = log10(rmin)
    rates%log_depop = log10(rates%log_depop)
  
    rmin = minval(rates%log_pop, rates%log_pop.gt.0.d0)
    where (rates%log_pop.le.0.d0)
        rates%log_pop = 0.9*rmin
    end where
    rates%minlog_pop = log10(rmin)
    rates%log_pop = log10(rates%log_pop)
  
end subroutine read_rates

subroutine read_tables
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(2) :: dim2
    integer :: error
   
    integer :: n_max, m_max
    character(len=4) :: impname
    real(double) :: imp_amu
    real(double), dimension(2) :: b_amu
    real(double), dimension(:,:), allocatable :: dummy2
  
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- Atomic tables settings ----"
    endif
  
    !!Initialize HDF5 interface
    call h5open_f(error)
  
    !!Open HDF5 file
    call h5fopen_f(inputs%tables_file, H5F_ACC_RDONLY_F, fid, error)
  
    !!Read Hydrogen-Hydrogen CX Cross Sections
    call read_cross(fid,"/cross/H_H",tables%H_H_cx)
  
    !!Read Hydrogen-Hydrogen Rates
    b_amu = [inputs%ab, inputs%ai]
    call read_rates(fid,"/rates/H_H",b_amu, inputs%ai, tables%H_H) 
    inputs%ab = tables%H_H%ab(1)
    inputs%ai = tables%H_H%ab(2)
  
    !!Read Hydrogen-Electron Rates
    call read_rates(fid,"/rates/H_e",b_amu, e_amu, tables%H_e)
  
    !!Read Hydrogen-Impurity rates
    impname = ''
    select case (inputs%impurity_charge)
        case (5)
            impname = "B5"
            imp_amu = B5_amu
        case (6) 
            impname = "C6"
            imp_amu = C6_amu
        case DEFAULT
            impname = "Aq"
            imp_amu = 2.d0*inputs%impurity_charge
    end select
  
    call read_rates(fid,"/rates/H_"//trim(adjustl(impname)), b_amu, imp_amu, tables%H_Aq)
  
    !!Read Einstein coefficients
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/n_max", n_max, error) 
    call h5ltread_dataset_int_scalar_f(fid,"/rates/spontaneous/m_max", m_max, error) 
    allocate(dummy2(n_max,m_max))
    dim2 = [n_max, m_max]
    call h5ltread_dataset_double_f(fid,"/rates/spontaneous/einstein",dummy2, dim2, error)
    tables%einstein(:,:) = transpose(dummy2(1:nlevs,1:nlevs))
    deallocate(dummy2)
  
    !!Close file
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Maximum n/m: ",i2)') nlevs
        write(*,'(T2,"Beam/Fast-ion mass: ",f6.3," [amu]")') inputs%ab 
        write(*,'(T2,"Thermal/Bulk-ion mass: ",f6.3," [amu]")') inputs%ai
        write(*,'(T2,"Impurity mass: ",f6.3," [amu]")') imp_amu
        write(*,*) ''
    endif

end subroutine read_tables

subroutine write_beam_grid(id, error)
    integer(HID_T), intent(inout) :: id
    integer, intent(out)          :: error
    
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(3) :: dims
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: u_grid, v_grid, w_grid
    real(double) :: xyz(3),uvw(3)
    integer :: i,j,k
  
    !Create uvw grids
    do k=1, beam_grid%nz
        do j=1, beam_grid%ny
            do i=1, beam_grid%nx
                xyz = [beam_grid%xc(i), &
                       beam_grid%yc(j), &
                       beam_grid%zc(k)]
                call xyz_to_uvw(xyz,uvw)
                u_grid(i,j,k) = uvw(1)
                v_grid(i,j,k) = uvw(2)
                w_grid(i,j,k) = uvw(3)
           enddo
        enddo
    enddo
  
    !Create grid group
    call h5gcreate_f(id, "grid", gid, error)
  
    !Write variables
    dims(1) = 1
    call h5ltmake_dataset_int_f(gid,"nx", 0, dims(1:1), [beam_grid%nx], error)
    call h5ltmake_dataset_int_f(gid,"ny", 0, dims(1:1), [beam_grid%ny], error)
    call h5ltmake_dataset_int_f(gid,"nz", 0, dims(1:1), [beam_grid%nz], error)
  
    dims = [beam_grid%nx, beam_grid%ny, beam_grid%nz]
    call h5ltmake_compressed_dataset_double_f(gid,"x", 1, dims(1:1), beam_grid%xc, error)
    call h5ltmake_compressed_dataset_double_f(gid,"y", 1, dims(2:2), beam_grid%yc, error)
    call h5ltmake_compressed_dataset_double_f(gid,"z", 1, dims(3:3), beam_grid%zc, error)
  
    call h5ltmake_compressed_dataset_double_f(gid,"x_grid", 3, dims, u_grid, error)
    call h5ltmake_compressed_dataset_double_f(gid,"y_grid", 3, dims, v_grid, error)
    call h5ltmake_compressed_dataset_double_f(gid,"z_grid", 3, dims, w_grid, error)
  
    !Write attributes
    call h5ltset_attribute_string_f(gid,"nx","description", &
         "Number of cells in the X direction", error)
    call h5ltset_attribute_string_f(gid,"ny","description", &
         "Number of cells in the Y direction", error)
    call h5ltset_attribute_string_f(gid,"nz","description", &
         "Number of cells in the Z direction", error)
    call h5ltset_attribute_string_f(gid,"x","description", &
         "X value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"x","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"y","description", &
         "Y value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"y","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"z","description", &
         "Z value of cell center in beam grid coordinates", error)
    call h5ltset_attribute_string_f(gid,"z","units", "cm", error)
  
    call h5ltset_attribute_string_f(gid,"x_grid","description", &
         "X value of cell center in machine coordinates: x_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"x_grid","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"y_grid","description", &
         "Y value of cell center in machine coordinates: y_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"y_grid","units", "cm", error)
    call h5ltset_attribute_string_f(gid,"z_grid","description", &
         "Z value of cell center in machine coordinates: z_grid(x,y,z)", error)
    call h5ltset_attribute_string_f(gid,"z_grid","units", "cm", error)
  
    call h5ltset_attribute_string_f(id,"grid","coordinate_system", &
         "Right-handed cartesian",error)
  
    !Close grid group 
    call h5gclose_f(gid, error)
  
end subroutine write_beam_grid

subroutine write_birth_profile
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error, i
  
    character(120)    :: filename
    real(double), dimension(:,:), allocatable :: ri
    real(double), dimension(:,:), allocatable :: vi
    real(double), dimension(3) :: xyz,uvw,v_uvw
  
    allocate(ri(3,size(birth%ri,2)))
    allocate(vi(3,size(birth%vi,2)))
  
    do i=1,size(birth%ri,2)
        ! Convert position to rzphi
        xyz = birth%ri(:,i)
        call xyz_to_uvw(xyz,uvw)
        ri(1,i) = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        ri(2,i) = uvw(3)
        ri(3,i) = atan2(uvw(2),uvw(1))
        ! Convert velocity to rzphi
        v_uvw = matmul(beam_grid%basis, birth%vi(:,i))
        vi(1,i) = v_uvw(1)*cos(ri(3,i)) + v_uvw(2)*sin(ri(3,i))
        vi(2,i) = v_uvw(3)
        vi(3,i) = -v_uvw(1)*sin(ri(3,i)) + v_uvw(2)*cos(ri(3,i))
    enddo
  
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_birth.h5"
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    call write_beam_grid(fid, error)
    d(1) = 1
    call h5ltmake_dataset_int_f(fid, "/n_nbi", 0, d, [inputs%n_nbi], error)
    call h5ltmake_dataset_int_f(fid, "/n_birth", 0, d, [inputs%n_birth], error)
    dim4 = shape(birth%dens)
    call h5ltmake_compressed_dataset_double_f(fid,"/dens", 4, dim4, birth%dens, error)
    dim2 = shape(birth%ri)
    call h5ltmake_compressed_dataset_double_f(fid,"/ri", 2, dim2, ri, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vi", 2, dim2, vi, error)
  
    !Add attributes
    call h5ltset_attribute_string_f(fid, "/n_nbi", "description", &
         "Number of beam mc particles", error)
    call h5ltset_attribute_string_f(fid, "/n_birth","description", &
         "Number of birth mc particles deposited per beam mc particle", error)
    call h5ltset_attribute_string_f(fid, "/dens", "description", &
         "Birth density: dens(beam_component,x,y,z)", error)
    call h5ltset_attribute_string_f(fid, "/dens", "units", &
         "fast-ions/(s*cm^3*dE*dP)", error)
    call h5ltset_attribute_string_f(fid, "/ri", "description", &
         "Fast-ion birth position in R-Z-Phi: ri([r,z,phi],particle)", error)
    call h5ltset_attribute_string_f(fid, "/ri", "units", "cm, radians", error)
    call h5ltset_attribute_string_f(fid, "/vi", "description", &
         "Fast-ion birth velocity in R-Z-Phi: vi([r,z,phi],particle)", error)
    call h5ltset_attribute_string_f(fid, "/vi", "units", "cm/s", error)
    call h5ltset_attribute_string_f(fid, "/", "coordinate_system", &
         "Cylindrical (R,Z,Phi)",error)
    call h5ltset_attribute_string_f(fid, "/", "description", &
         "Birth density and particles calculated by FIDASIM", error)
  
    !!Close file
    call h5fclose_f(fid, error)
  
    !!Close HDF5 interface
    call h5close_f(error)
    
    deallocate(ri,vi)
    if(inputs%verbose.ge.1) then 
        write(*,'(T4,a,a)') 'birth profile written to: ',trim(filename)
    endif

end subroutine write_birth_profile

subroutine write_dcx
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error
  
    character(120)    :: filename
  
    integer         :: i
    real(double), dimension(:)  , allocatable :: lambda_arr
    real(double), dimension(:,:), allocatable :: dcx_spec
  
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_dcx.h5"
  
    allocate(lambda_arr(inputs%nlambda))
    do i=1,inputs%nlambda 
        lambda_arr(i) = (i-0.5)*inputs%dlambda*0.1d0 &
                      + inputs%lambdamin*0.1d0 ! [nm]
    enddo 
  
    allocate(dcx_spec(inputs%nlambda,spec_chords%nchan))
  
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    dcx_spec = spec%bes(:,:,halo_type)/(0.1d0*inputs%dlambda)/(4.d0*pi)*1.d4
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    call write_beam_grid(fid, error)
    d(1) =1 
    call h5ltmake_dataset_int_f(fid,"/nlevel", 0, d, [nlevs], error)
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [spec_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nlambda", 0, d, [inputs%nlambda], error)
    dims(1) = inputs%nlambda
    dims(2) = spec_chords%nchan
    call h5ltmake_compressed_dataset_double_f(fid, "/spec", 2, dims(1:2), &
         dcx_spec, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/lambda", 1, dims(1:1), &
         lambda_arr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, dims(2:2), &
         spec_chords%radius, error)
    dims = [nlevs, beam_grid%nx, beam_grid%ny, beam_grid%nz ] 
    call h5ltmake_compressed_dataset_double_f(fid, "/dens", 4, dims, &
         neut%dens(:,halo_type,:,:,:), error)
  
    !Add attributes
    call h5ltset_attribute_string_f(fid,"/nlevel","description", &
         "Number of atomic energy levels", error)
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nlambda", "description", &
         "Number of wavelengths", error)
    call h5ltset_attribute_string_f(fid,"/dens", "description", &
         "Direct Charge Exchange (DCX) neutral density: dcx(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/dens","units","neutrals*cm^-3", error)
    call h5ltset_attribute_string_f(fid,"/spec","description", &
         "Direct Charge Exchange (DCX) beam emission: spec(lambda, chan)", error)
    call h5ltset_attribute_string_f(fid,"/spec","units","Ph/(s*nm*sr*m^2)",error)
    call h5ltset_attribute_string_f(fid,"/lambda","description", &
         "Wavelength array", error) 
    call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius","units","cm", error)
  
    call h5ltset_attribute_string_f(fid,"/","description", &
         "Direct Charge Exchange spectra and neutral density calculated by FIDASIM", error)
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
  
    deallocate(dcx_spec,lambda_arr)
    
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'dcx written to: ',trim(filename)
    endif

end subroutine write_dcx

subroutine write_neutrals
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(inputs%neutrals_file, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    call write_beam_grid(fid, error)
    dims = [nlevs, beam_grid%nx, beam_grid%ny, beam_grid%nz]
    d(1) =1 
    call h5ltmake_dataset_int_f(fid,"/nlevel", 0, d, [nlevs], error)
    call h5ltmake_compressed_dataset_double_f(fid, "/fdens", 4, dims, &
         neut%dens(:,nbif_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/hdens", 4, dims, &
         neut%dens(:,nbih_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/tdens", 4, dims, &
         neut%dens(:,nbit_type,:,:,:), error)
    call h5ltmake_compressed_dataset_double_f(fid, "/halodens", 4, dims, &
         neut%dens(:,halo_type,:,:,:), error)
  
    !Write attributes
    call h5ltset_attribute_string_f(fid,"/nlevel","description", &
         "Number of atomic energy levels", error)
    call h5ltset_attribute_string_f(fid,"/fdens","description", &
         "Neutral density for the full energy component of the beam: fdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/fdens","units","neutrals*cm^-3",error)
    
    call h5ltset_attribute_string_f(fid,"/hdens","description", &
         "Neutral density for the half energy component of the beam: hdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/hdens","units","neutrals*cm^-3",error)
  
    call h5ltset_attribute_string_f(fid,"/tdens","description", &
         "Neutral density for the third energy component of the beam: tdens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/tdens","units","neutrals*cm^-3",error)
  
    call h5ltset_attribute_string_f(fid,"/halodens","description", &
         "Neutral density of the beam halo(including dcx): halodens(level,x,y,z)", error)
    call h5ltset_attribute_string_f(fid,"/halodens","units","neutrals*cm^-3",error)
  
    call h5ltset_attribute_string_f(fid,"/","description", &
         "Beam neutral density calculated by FIDASIM", error)
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
  
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'neutral density written to: ',trim(inputs%neutrals_file)
    endif

end subroutine write_neutrals

subroutine write_npa
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(2) :: dim2 
    integer(HSIZE_T), dimension(1) :: d
    integer :: error
  
    integer, dimension(:), allocatable :: dcount
    real(double), dimension(:,:), allocatable :: ri,rf
    integer :: i, n
    character(120)    :: filename = ''
  
    allocate(dcount(npa_chords%nchan))
    do i=1,npa_chords%nchan
        dcount(i) = count(npa%part%detector.eq.i)
    enddo
  
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa.h5"
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    d(1) = 1
    dim2 = [fbm%nenergy, npa%nchan]
    call h5ltmake_dataset_int_f(fid,"/nenergy", 0, d, [fbm%nenergy], error)
    call h5ltmake_dataset_int_f(fid,"/nchan", 0, d, [npa%nchan], error)
    call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim2(1:1),&
         npa%energy, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim2(2:2),&
         npa_chords%radius, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/flux",2,dim2,npa%flux, error)
    call h5ltmake_compressed_dataset_int_f(fid,"/count",1,dim2(2:2), dcount, error)
  
    !Add attributes
    call h5ltset_attribute_string_f(fid,"/","description", &
         "NPA flux calculated by FIDASIM",error)
    call h5ltset_attribute_string_f(fid,"/nenergy","description",&
         "Number of energy values",error)
    call h5ltset_attribute_string_f(fid,"/nchan","description",&
         "Number of channels",error)
    call h5ltset_attribute_string_f(fid,"/energy","description", &
         "Energy array", error)
    call h5ltset_attribute_string_f(fid,"/energy","units","keV", error)
    call h5ltset_attribute_string_f(fid,"/radius","description", &
         "Detector line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius","units","cm",error)
    call h5ltset_attribute_string_f(fid,"/flux", "description", &
         "Neutral flux: flux(energy,chan)", error)
    call h5ltset_attribute_string_f(fid,"/flux", "units","neutrals/(s*dE)", error)
    call h5ltset_attribute_string_f(fid,"/count","description", &
         "Number of particles that hit the detector: count(chan)", error)
  
    deallocate(dcount)
  
    if((npa%npart.ne.0).and.(inputs%calc_npa.ge.2)) then
        n = npa%npart
        allocate(ri(3,n),rf(3,n))
        ri(1,:) = npa%part(1:n)%xi
        ri(2,:) = npa%part(1:n)%yi
        ri(3,:) = npa%part(1:n)%zi
        rf(1,:) = npa%part(1:n)%xf
        rf(2,:) = npa%part(1:n)%yf
        rf(3,:) = npa%part(1:n)%zf
  
        !Create Group
        call h5gcreate_f(fid,"/particles",gid, error)
        call h5ltmake_dataset_int_f(gid, "nparticle", 0, d, [npa%npart], error)
        d(1) = npa%npart
        dim2 = [3, n]
        call h5ltmake_compressed_dataset_double_f(gid,"ri",2,dim2, ri, error)
        call h5ltmake_compressed_dataset_double_f(gid,"rf",2,dim2, rf, error)
        call h5ltmake_compressed_dataset_double_f(gid,"pitch",1,d, &
             npa%part(1:n)%pitch, error)
        call h5ltmake_compressed_dataset_double_f(gid,"energy",1,d,&
             npa%part(1:n)%energy, error)
        call h5ltmake_compressed_dataset_double_f(gid,"weight",1,d,&
             npa%part(1:n)%weight, error)
        call h5ltmake_compressed_dataset_int_f(gid,"detector",1,d,&
             npa%part(1:n)%detector, error)
  
        !Add attributes
        call h5ltset_attribute_string_f(gid,"nparticle","description", &
             "Number of particles that hit a detector", error)
        call h5ltset_attribute_string_f(gid,"ri","description", &
             "Neutral particle's birth position in machine coordinates: ri([x,y,z],particle)", error)
        call h5ltset_attribute_string_f(gid,"ri","units", "cm", error)
        call h5ltset_attribute_string_f(gid,"rf","description", &
             "Neutral particle's hit position in machine coordinates: rf([x,y,z],particle)", error)
        call h5ltset_attribute_string_f(gid,"rf","units", "cm", error)
        call h5ltset_attribute_string_f(gid,"pitch","description", &
             "Pitch value of the neutral particle: p = v_parallel/v  w.r.t. the magnetic field", error)
        call h5ltset_attribute_string_f(gid,"energy","description", &
             "Energy value of the neutral particle", error)
        call h5ltset_attribute_string_f(gid,"energy","units","keV",error)
        call h5ltset_attribute_string_f(gid,"weight","description", &
             "Neutral particle's contribution to the flux", error)
        call h5ltset_attribute_string_f(gid,"weight","units","neutrals/s",error)
        call h5ltset_attribute_string_f(gid,"detector","description", &
             "Detector that the neutral particle hit", error)
  
        call h5ltset_attribute_string_f(fid,"/particles","coordinate_system", &
             "Right-handed cartesian",error)
        call h5ltset_attribute_string_f(fid,"/particles","description", &
             "Monte Carlo particles",error)
  
        !Close group
        call h5gclose_f(gid, error)
        deallocate(ri,rf)
    endif
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
  
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'NPA data written to: ',trim(filename)
    endif

end subroutine write_npa

subroutine write_spectra
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(3) :: dims
    integer(HSIZE_T), dimension(1) :: d
    integer :: error
  
    character(120)    :: filename
    integer         :: i
    real(double), dimension(:)  , allocatable :: lambda_arr
  
    allocate(lambda_arr(inputs%nlambda))
    do i=1,inputs%nlambda
        lambda_arr(i) = (i-0.5)*inputs%dlambda*0.1d0 &
                      + inputs%lambdamin*0.1d0
    enddo
  
    !! convert [Ph/(s*wavel_bin*cm^2*all_directions)] to [Ph/(s*nm*sr*m^2)]!
    spec%brems=spec%brems/(0.1d0*inputs%dlambda)/(4.d0*pi)*1.d4
    spec%bes=spec%bes/(0.1d0*inputs%dlambda)/(4.d0*pi)*1.d4
    spec%fida=spec%fida/(0.1d0*inputs%dlambda)/(4.d0*pi)*1.d4
  
    !! write to file
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_spectra.h5"
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    d(1) = 1
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [spec_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nlambda", 0, d, [inputs%nlambda], error)
    dims(1) = inputs%nlambda
    dims(2) = spec_chords%nchan
    dims(3) = particles%nclass
    call h5ltmake_compressed_dataset_double_f(fid, "/lambda", 1, dims(1:1), &
         lambda_arr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, dims(2:2), &
         spec_chords%radius, error)
  
    !Add attributes
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nlambda", "description", &
         "Number of wavelengths", error)
    call h5ltset_attribute_string_f(fid,"/lambda","description", &
         "Wavelength array", error) 
    call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius","units","cm", error)
  
    if(inputs%calc_brems.ge.1) then 
        !Write variables
        call h5ltmake_compressed_dataset_double_f(fid, "/brems", 2, &
             dims(1:2), spec%brems, error)
        !Add attributes
        call h5ltset_attribute_string_f(fid,"/brems","description", &
             "Visible Bremsstrahlung: brems(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/brems","units",&
             "Ph/(s*nm*sr*m^2)",error )
    endif
  
    if(inputs%calc_bes.ge.1) then
        !Write variables
        call h5ltmake_compressed_dataset_double_f(fid, "/full", 2, dims(1:2), &
             spec%bes(:,:,nbif_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/half", 2, dims(1:2), &
             spec%bes(:,:,nbih_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/third", 2, dims(1:2),&
             spec%bes(:,:,nbit_type), error)
        call h5ltmake_compressed_dataset_double_f(fid, "/halo", 2, dims(1:2), &
             spec%bes(:,:,halo_type), error)
        !Add attributes
        call h5ltset_attribute_string_f(fid,"/full","description", &
             "Full energy component of the beam emmision: full(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/full","units","Ph/(s*nm*sr*m^2)",error )
        call h5ltset_attribute_string_f(fid,"/half","description", &
             "Half energy component of the beam emmision: half(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/half","units","Ph/(s*nm*sr*m^2)",error )
        call h5ltset_attribute_string_f(fid,"/third","description", &
             "Third energy component of the beam emmision: third(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/third","units","Ph/(s*nm*sr*m^2)",error )
        call h5ltset_attribute_string_f(fid,"/halo","description", &
             "Halo component of the beam emmision (includes dcx): halo(lambda,chan)", error)
        call h5ltset_attribute_string_f(fid,"/halo","units","Ph/(s*nm*sr*m^2)",error )
    endif
  
    if(inputs%calc_fida.ge.1) then
        !Write variables
        if(particles%nclass.le.1) then
            call h5ltmake_compressed_dataset_double_f(fid, "/fida", 2, &
                 dims(1:2), spec%fida(:,:,1), error)
            !Add attributes
            call h5ltset_attribute_string_f(fid,"/fida","description", &
                 "Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan)", error)
        else
            call h5ltmake_dataset_int_f(fid,"/nclass", 0, d, [particles%nclass], error)
            call h5ltmake_compressed_dataset_double_f(fid, "/fida", 2, &
                 dims, spec%fida, error)
            !Add attributes
            call h5ltset_attribute_string_f(fid,"/fida","description", &
                 "Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan,class)", error)
       endif
        call h5ltset_attribute_string_f(fid,"/fida","units","Ph/(s*nm*sr*m^2)",error )
    endif
  
    call h5ltset_attribute_string_f(fid,"/","description",&
         "Spectra calculated by FIDASIM", error)
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
  
    if(inputs%verbose.ge.1) then 
        write(*,'(T4,a,a)') 'Spectra written to: ', trim(filename)
    endif

end subroutine write_spectra

subroutine write_fida_weights
    !! HDF5 variables
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: dim1
    integer :: error
  
    character(120) :: filename
    integer :: i,ie,ip,ic
    real(double), dimension(:),   allocatable :: lambda_arr
    real(double), dimension(:),   allocatable :: ebarr,ptcharr
    real(double), dimension(:,:), allocatable :: jacobian,e_grid,p_grid
    real(double), dimension(:,:), allocatable :: vpa_grid,vpe_grid
    real(double) :: dlambda, wtot, dE, dP
  
    dlambda=(inputs%lambdamax_wght-inputs%lambdamin_wght)/inputs%nlambda_wght
    allocate(lambda_arr(inputs%nlambda_wght))
    do i=1,inputs%nlambda_wght
        lambda_arr(i)=(i-0.5)*dlambda + inputs%lambdamin_wght
    enddo
  
    !! define arrays
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))
  
    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))
  
    !! define 2d grids
    !! define energy grid
    allocate(e_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%ne_wght
        e_grid(i,:) = ebarr(i)
    enddo
  
    !! define pitch grid
    allocate(p_grid(inputs%ne_wght,inputs%np_wght))
    do i=1,inputs%np_wght
        p_grid(:,i) = ptcharr(i)
    enddo
  
    !! define velocity space grid
    allocate(vpe_grid(inputs%ne_wght,inputs%np_wght)) !! V perpendicular
    allocate(vpa_grid(inputs%ne_wght,inputs%np_wght)) !! V parallel
    vpa_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid)*p_grid ! [cm/s]
    vpe_grid = 100*sqrt((((2.0d3)*e0)/(mass_u*inputs%ab))*e_grid*(1.0-p_grid**2.0)) ![cm/s]
  
    !! define jacobian to convert between E-p to velocity
    allocate(jacobian(inputs%ne_wght,inputs%np_wght))
    jacobian = ((inputs%ab*mass_u)/(e0*1.0d3)) *vpe_grid/sqrt(vpa_grid**2.0 + vpe_grid**2.0)
  
    !! normalize mean_f
    do ic=1,spec_chords%nchan
        do ip=1,inputs%np_wght
            do ie=1,inputs%ne_wght
                wtot = sum(fweight%weight(:,ie,ip,ic))
                if((wtot.gt.0.d0)) then
                    fweight%mean_f(ie,ip,ic) = fweight%mean_f(ie,ip,ic)/wtot
                endif
            enddo
        enddo
    enddo
    if(inputs%calc_fida_wght.eq.1) then
        fweight%mean_f = fweight%mean_f/(dE*dP)
    endif
  
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_fida_weights.h5"
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    dim1(1) = 1
    dim2 = [inputs%nlambda_wght, spec_chords%nchan]
    dim4 = [inputs%nlambda_wght, inputs%ne_wght, inputs%np_wght, spec_chords%nchan]
    call h5ltmake_dataset_int_f(fid,"/nenergy",0,dim1,[inputs%ne_wght], error)
    call h5ltmake_dataset_int_f(fid,"/npitch",0,dim1,[inputs%np_wght], error)
    call h5ltmake_dataset_int_f(fid,"/nchan",0,dim1,[spec_chords%nchan], error)
    call h5ltmake_compressed_dataset_double_f(fid,"/weight",4,dim4,fweight%weight,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/fida",2,dim2,fweight%fida,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/mean_f",3,dim4(2:4),fweight%mean_f,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/lambda",1,dim4(1:1),lambda_arr,error)
    call h5ltmake_compressed_dataset_double_f(fid,"/energy",1,dim4(2:2),ebarr, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/pitch",1,dim4(3:3),ptcharr, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/radius",1,dim4(4:4),spec_chords%radius, error)
    dim2 = [inputs%ne_wght, inputs%np_wght]
    call h5ltmake_compressed_dataset_double_f(fid,"/jacobian",2,dim2, jacobian, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vpe_grid",2,dim2,vpe_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/vpa_grid",2,dim2,vpa_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/e_grid",2,dim2,e_grid, error)
    call h5ltmake_compressed_dataset_double_f(fid,"/p_grid",2,dim2,p_grid, error)
  
    !Add attributes
    if(inputs%calc_fida_wght.eq.1) then 
        call h5ltset_attribute_string_f(fid,"/", "description", &
             "Line of Sight averaged FIDA E-p space sensitivity/weights " // &
             "and spectra calculated by FIDASIM", error)
    else
        call h5ltset_attribute_string_f(fid,"/", "description", &
             "Full FIDA E-p space sensitivity/weights and spectra calculated " // &
             "by FIDASIM via Monte Carlo method", error)
    endif
    call h5ltset_attribute_string_f(fid,"/weight","description", &
         "E-p space sensivity/weight of FIDA diagnostic: weight(lambda,energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/weight","units", &
         "(Ph*cm)/(s*nm*sr*fast-ion*dE*dP)",error)
    call h5ltset_attribute_string_f(fid,"/fida","units", &
         "Ph/(s*nm*sr*m^2)",error )
    call h5ltset_attribute_string_f(fid,"/fida","description", &
         "Fast-ion D-alpha (FIDA) emmision: fida(lambda,chan)", error)
    call h5ltset_attribute_string_f(fid,"/mean_f","description", &
         "Mean fast-ion distribution function seen by los: mean_f(energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/mean_f","units", &
         "fast-ion/(dE*dP*cm^3)", error)
    call h5ltset_attribute_string_f(fid,"/lambda","description", &
         "Wavelength array", error) 
    call h5ltset_attribute_string_f(fid,"/lambda","units","nm", error)
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(fid,"/npitch", "description", &
         "Number of pitch value", error)
    call h5ltset_attribute_string_f(fid,"/energy","description", &
         "Energy array", error) 
    call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
    call h5ltset_attribute_string_f(fid,"/pitch", "description", &
         "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error) 
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius", "units","cm", error)
    call h5ltset_attribute_string_f(fid,"/jacobian","description", &
         "Jacobian used to convert from E-p space to velocity space", error)
    call h5ltset_attribute_string_f(fid,"/jacobian","units", &
         "(dE*dP)/(dvpa*dvpe)", error)
    call h5ltset_attribute_string_f(fid,"/e_grid","description", &
         "2D energy grid", error)
    call h5ltset_attribute_string_f(fid,"/e_grid","units","keV", error)
    call h5ltset_attribute_string_f(fid,"/p_grid","description", &
         "2D pitch grid", error)
    call h5ltset_attribute_string_f(fid,"/vpe_grid","description", &
         "2D perpendicular velocity grid", error)
    call h5ltset_attribute_string_f(fid,"/vpe_grid","units","cm/s", error)
    call h5ltset_attribute_string_f(fid,"/vpa_grid","description", &
         "2D parallel velocity grid", error)
    call h5ltset_attribute_string_f(fid,"/vpa_grid","units","cm/s", error)
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
  
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'FIDA weights written to: ', trim(filename)
    endif

end subroutine write_fida_weights

subroutine write_npa_weights
    character(120) :: filename
    integer :: i
    real(double), dimension(:), allocatable :: ebarr,ptcharr
   
    !! HDF5 variables
    integer(HID_T) :: fid
    integer(HSIZE_T), dimension(5) :: dim5
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(1) :: d
    integer :: error
  
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
  
    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
  
    filename=trim(adjustl(inputs%result_dir))//"/"//trim(adjustl(inputs%runid))//"_npa_weights.h5"
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fid, error)
  
    !Write variables
    d(1) = 1
    dim2 = [inputs%ne_wght, npa_chords%nchan]
    dim3 = [inputs%ne_wght, inputs%np_wght, npa_chords%nchan]
    dim5 = [inputs%ne_wght, beam_grid%nx, beam_grid%ny, beam_grid%nz, npa_chords%nchan]
    call h5ltmake_dataset_int_f(fid, "/nchan", 0, d, [npa_chords%nchan], error)
    call h5ltmake_dataset_int_f(fid, "/nenergy", 0, d, [inputs%ne_wght], error)
    call h5ltmake_dataset_int_f(fid, "/npitch", 0, d, [inputs%np_wght], error)
    call h5ltmake_compressed_dataset_double_f(fid, "/radius", 1, &
         dim2(2:2), npa_chords%radius, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/energy", 1, &
         dim2(1:1), ebarr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/pitch", 1, &
         dim3(2:2), ptcharr, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/flux", 2, &
         dim2, nweight%flux, error)
    call h5ltmake_compressed_dataset_double_f(fid, "/weight", 3, &
         dim3, nweight%weight, error)
  
    !Add attributes
    call h5ltset_attribute_string_f(fid,"/", "description", &
         "NPA E-p space sensitivity/weights and Flux calculated by FIDASIM", error)
    call h5ltset_attribute_string_f(fid,"/nchan", "description", &
         "Number of channels", error)
    call h5ltset_attribute_string_f(fid,"/nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(fid,"/npitch", "description", &
         "Number of pitch value", error)
    call h5ltset_attribute_string_f(fid,"/energy","description", &
         "Energy array", error) 
    call h5ltset_attribute_string_f(fid,"/energy", "units","keV", error)
    call h5ltset_attribute_string_f(fid,"/pitch", "description", &
         "Pitch array: p = v_parallel/v  w.r.t. the magnetic field", error) 
    call h5ltset_attribute_string_f(fid,"/radius", "description", &
         "Line of sight radius at midplane or tangency point", error)
    call h5ltset_attribute_string_f(fid,"/radius", "units","cm", error)
    call h5ltset_attribute_string_f(fid,"/flux", "description", &
         "Neutral flux: flux(energy,chan)", error)
    call h5ltset_attribute_string_f(fid,"/flux", "units", &
         "neutrals/(s*dE)", error)
    call h5ltset_attribute_string_f(fid,"/weight", "description", &
         "E-p space sensivity/weight of NPA diagnostics: weight(energy,pitch,chan)", error)
    call h5ltset_attribute_string_f(fid,"/weight","units", &
         "neutrals/(s*fast-ion*dE*dP)",error)
  
    if(inputs%calc_npa_wght.ge.2) then !Write diagnostic variables
        call write_beam_grid(fid, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/emissivity", 4, &
             dim5(2:5), nweight%emissivity, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/attenuation", 5, &
             dim5, nweight%attenuation, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/cx", 5, &
             dim5, nweight%cx, error)
        call h5ltmake_compressed_dataset_double_f(fid, "/phit", 4, &
             dim5(2:5), npa_chords%phit%p, error)
  
        call h5ltset_attribute_string_f(fid,"/emissivity", "description", &
             "Neutral emissivity: emissivity(x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/emissivity", "units", &
             "neutrals/(s*dV)", error)
        call h5ltset_attribute_string_f(fid,"/cx", "description", &
             "Charge-exchange rate: cx(energy,x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/cx", "units", "s^(-1)", error)
        call h5ltset_attribute_string_f(fid,"/attenuation","description", &
             "Attenuation factor i.e. survival probability: attenuation(energy,x,y,z,chan)", error)
        call h5ltset_attribute_string_f(fid,"/phit","description", &
             "Probability of hitting the detector given an isotropic source: phit(x,y,z,chan)", error)
    endif
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
    if(inputs%verbose.ge.1) then
        write(*,'(T4,a,a)') 'NPA weights written to: ',trim(filename)
    endif
  
end subroutine write_npa_weights

subroutine read_neutrals
    integer(HID_T) :: fid, gid
    integer(HSIZE_T), dimension(4) :: dims
    integer :: error,nx,ny,nz
    logical :: exis
  
    if(inputs%verbose.ge.1) then 
        write(*,'(a)') '---- loading neutrals ----'
    endif
  
    inquire(file=inputs%neutrals_file,exist=exis)
    if(exis) then
        write(*,'(T2,"Neutrals file: ",a)') trim(inputs%neutrals_file)
        write(*,*) ''
    else
        write(*,'(a,a)') 'READ_NEUTRALS: Neutrals file does not exist: ',inputs%neutrals_file
        stop
    endif
  
    !Open HDF5 interface
    call h5open_f(error)
    
    !Create file overwriting any existing file
    call h5fopen_f(inputs%neutrals_file, H5F_ACC_RDONLY_F, fid, error)
  
    call h5gopen_f(fid, "/grid", gid, error)
    call h5ltread_dataset_int_scalar_f(gid,"nx", nx, error) 
    call h5ltread_dataset_int_scalar_f(gid,"ny", ny, error) 
    call h5ltread_dataset_int_scalar_f(gid,"nz", nz, error)
    call h5gclose_f(gid, error)
  
    if((nx.ne.beam_grid%nx).or. &
       (ny.ne.beam_grid%ny).or. &
       (nz.ne.beam_grid%nz)) then
        write(*,'(a)') 'READ_NEUTRALS: Neutrals file has incompatable grid dimensions'
        stop
    endif
  
    dims = [nlevs, nx, ny, nz]
    call h5ltread_dataset_double_f(fid,"/fdens", &
         neut%dens(:,nbif_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/hdens", &
         neut%dens(:,nbih_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/tdens", &
         neut%dens(:,nbit_type,:,:,:), dims, error)
    call h5ltread_dataset_double_f(fid,"/halodens", &
         neut%dens(:,halo_type,:,:,:), dims, error)
  
    !Close file
    call h5fclose_f(fid, error)
  
    !Close HDF5 interface
    call h5close_f(error)
    
end subroutine read_neutrals

!=============================================================================
!-----------------------------Geometry Routines-------------------------------
!=============================================================================
function cross_product(u, v) result(s)
    real(double), dimension(3), intent(in) :: u
    real(double), dimension(3), intent(in) :: v
    real(double), dimension(3)             :: s
  
    s(1) = u(2)*v(3) - u(3)*v(2)
    s(2) = u(3)*v(1) - u(1)*v(3)
    s(3) = u(1)*v(2) - u(2)*v(1)

end function cross_product

function normp(u, p_in) result(n)
    real(double), dimension(:), intent(in) :: u
    integer, intent(in), optional          :: p_in
    real(double)                           :: n
  
    integer :: p
    
    IF(present(p_in)) THEN
        p = p_in
    ELSE
        p = 2
    ENDIF
  
    SELECT CASE (p)
        CASE (1)
            n = sum(abs(u))
        CASE (2)
            n = sqrt(dot_product(u,u))
        CASE DEFAULT
            write(*,'("NORMP: Unknown p value: ",i2)'),p
            stop
    END SELECT

end function normp

subroutine tb_zyx(alpha, beta, gamma, basis, inv_basis)
    !! Creates active rotation matrix for z-y'-x" rotation given Tait-bryan angles
    real(double), intent(in)                            :: alpha
    real(double), intent(in)                            :: beta
    real(double), intent(in)                            :: gamma
    real(double), dimension(3,3), intent(out)           :: basis
    real(double), dimension(3,3), intent(out), optional :: inv_basis
  
    real(double) :: sa, sb, sg, ca, cb, cg
  
    sa = sin(alpha) ; sb = sin(beta) ; sg = sin(gamma)
    ca = cos(alpha) ; cb = cos(beta) ; cg = cos(gamma)
  
    basis(1,1) = ca*cb ; basis(1,2) = ca*sb*sg - cg*sa ; basis(1,3) = sa*sg + ca*cg*sb
    basis(2,1) = cb*sa ; basis(2,2) = ca*cg + sa*sb*sg ; basis(2,3) = cg*sa*sb - ca*sg
    basis(3,1) = -sb   ; basis(3,2) = cb*sg            ; basis(3,3) = cb*cg
  
    if(present(inv_basis)) inv_basis = transpose(basis)
  
end subroutine tb_zyx

subroutine line_basis(r0, v0, basis, inv_basis)
    !calculates basis from a line with +x in the direction of line
    real(double), dimension(3), intent(in)              :: r0
    real(double), dimension(3), intent(in)              :: v0
    real(double), dimension(3,3), intent(out)           :: basis
    real(double), dimension(3,3), intent(out), optional :: inv_basis
  
    real(double), dimension(3) :: rf
    real(double) :: alpha, beta, dis
  
    rf = r0 + v0
    dis = sqrt(sum((rf - r0)**2.0))
    beta = asin((r0(3) - rf(3))/dis)
    alpha = atan2(rf(2)-r0(2),rf(1)-r0(1))
  
    call tb_zyx(alpha,beta,0.d0,basis)
  
    if(present(inv_basis)) inv_basis = transpose(basis)
  
end subroutine line_basis

subroutine plane_basis(center, redge, tedge, basis, inv_basis)
    !calculates basis from 3 points on a plane with +z being the plane normal
    real(double), dimension(3), intent(in)              :: center
    real(double), dimension(3), intent(in)              :: redge
    real(double), dimension(3), intent(in)              :: tedge
    real(double), dimension(3,3), intent(out)           :: basis
    real(double), dimension(3,3), intent(out), optional :: inv_basis
  
    real(double), dimension(3) :: u1,u2,u3
  
    u1 = (redge - center)
    u1 = u1/normp(u1)
    u2 = (tedge - center)
    u2 = u2/normp(u2)
    u3 = cross_product(u1,u2)
    u3 = u3/normp(u3)
  
    basis(:,1) = u1
    basis(:,2) = u2
    basis(:,3) = u3
  
    if(present(inv_basis)) inv_basis = transpose(basis)
  
end subroutine plane_basis

subroutine plane_intercept(l0, l, p0, n, p, t)
    real(double), dimension(3), intent(in)  :: l0 ! Point on line
    real(double), dimension(3), intent(in)  :: l  ! Ray of line
    real(double), dimension(3), intent(in)  :: p0 ! Point on plane
    real(double), dimension(3), intent(in)  :: n  ! Normal vector of plane
    real(double), dimension(3), intent(out) :: p  ! Line-plane intercept point
    real(double), intent(out)               :: t ! "time" to intercept
  
    t = dot_product(p0 - l0, n)/dot_product(l, n)
  
    p = l0 + t*l
  
end subroutine plane_intercept

function in_boundary(bplane, p) result(in_b)
    type(BoundedPlane), intent(in)         :: bplane
    real(double), dimension(3), intent(in) :: p
  
    real(double), dimension(3) :: pp
    real(double) :: hh, hw
    logical :: in_b
  
    hh = bplane%hh
    hw = bplane%hw
    pp = matmul(bplane%inv_basis, p - bplane%origin)
    in_b = .False.
    SELECT CASE (bplane%shape)
        CASE (1) !Rectangular boundary
            if((abs(pp(1)).le.hw).and. &
               (abs(pp(2)).le.hh)) then
                in_b = .True.
            endif
        CASE (2) !Circular/Ellipsoidal boundary
            if(((hh*pp(1))**2.0 + (hw*pp(2))**2.0).le.((hh*hw)**2.0)) then
                in_b = .True.
            endif
        CASE DEFAULT
            write(*,'("IN_BOUNDARY: Unknown boundary shape: ",i2)'),bplane%shape
            stop
    END SELECT
  
end function in_boundary

subroutine hit_npa_detector(r0, v0, d_index, rd)
    real(double), dimension(3), intent(in)            :: r0
    real(double), dimension(3), intent(in)            :: v0
    integer, intent(out)                              :: d_index
    real(double), dimension(3), intent(out), optional :: rd
  
    real(double), dimension(3) :: d, a
    real(double) :: t_a,t_d
    integer :: i, ndet
  
    ndet = npa_chords%nchan
  
    d_index = 0
    detector_loop: do i=1,ndet
        !! Find where trajectory crosses detector plane
        call plane_intercept(r0,v0,npa_chords%det(i)%detector%origin, &
             npa_chords%det(i)%detector%basis(:,3),d,t_d)
  
        !! Find where trajectory crosses aperture plane
        call plane_intercept(r0,v0,npa_chords%det(i)%aperture%origin, &
             npa_chords%det(i)%aperture%basis(:,3),a,t_a)
  
        !! If both points are in plane boundaries and the 
        !! particle is heading toward the detector then its a hit
        if(in_boundary(npa_chords%det(i)%aperture,a) .and. &
           in_boundary(npa_chords%det(i)%detector,d) .and. &
           (t_d.gt.0.0) ) then
            d_index = i
            if(present(rd)) rd = d
            exit detector_loop
        endif
    enddo detector_loop

end subroutine hit_npa_detector

subroutine xyz_to_uvw(xyz, uvw)
    real(double), dimension(3), intent(in)  :: xyz
    real(double), dimension(3), intent(out) :: uvw
  
    real(double), dimension(3) :: origin
    real(double), dimension(3,3) :: basis
  
    origin = beam_grid%origin
    basis = beam_grid%basis
    uvw = matmul(basis,xyz)
    uvw = uvw + origin
  
end subroutine xyz_to_uvw

subroutine uvw_to_xyz(uvw,xyz)
    real(double), dimension(3), intent(in)  :: uvw
    real(double), dimension(3), intent(out) :: xyz
  
    real(double), dimension(3) :: origin, uvw_p
    real(double), dimension(3,3) :: basis
  
    origin = beam_grid%origin
    basis = beam_grid%inv_basis
    uvw_p = uvw - origin
    xyz = matmul(basis,uvw_p)
  
end subroutine uvw_to_xyz

subroutine grid_intersect(r0, v0, length, r_enter, r_exit, center_in, lwh_in)
    real(double), dimension(3), intent(in)           :: r0
    real(double), dimension(3), intent(in)           :: v0
    real(double), intent(out)                        :: length
    real(double), dimension(3), intent(out)          :: r_enter
    real(double), dimension(3), intent(out)          :: r_exit
    real(double), dimension(3), intent(in), optional :: center_in
    real(double), dimension(3), intent(in), optional :: lwh_in
  
    real(double), dimension(3,6) :: ipnts
    real(double), dimension(3) :: vi
    real(double), dimension(3) :: center
    real(double), dimension(3) :: lwh
    integer, dimension(6) :: side_inter
    integer, dimension(2) :: ind
    integer :: i, j, nunique, ind1, ind2
  
    if (present(center_in)) then
        center = center_in
    else
        center = beam_grid%center
    endif
  
    if (present(lwh_in)) then
        lwh = lwh_in
    else
        lwh = beam_grid%lwh
    endif
  
    side_inter = 0
    ipnts = 0.d0
    do i=1,6
        j = int(ceiling(i/2.0))
        if (j.eq.1) ind = [2,3]
        if (j.eq.2) ind = [1,3]
        if (j.eq.3) ind = [1,2]
        if (abs(v0(j)).gt.0.d0) then
            ipnts(:,i) = r0 + v0*( ( (center(j) + &
                         (mod(i,2) - 0.5)*lwh(j)) - r0(j))/v0(j) )
            if ((abs(ipnts(ind(1),i) - center(ind(1))).le.(0.5*lwh(ind(1)))).and. &
                (abs(ipnts(ind(2),i) - center(ind(2))).le.(0.5*lwh(ind(2))))) then
                side_inter(i) = 1
            endif
        endif
    enddo
  
    length = 0.d0
    r_enter = r0
    r_exit  = r0
  
    if (sum(side_inter).ge.2) then
        ! Find first intersection side
        i=1
        do while (i.le.6)
            if(side_inter(i).eq.1) exit
            i=i+1
        enddo
        ind1=i
        !Find number of unique points
        nunique = 0
        do i=ind1,6
            if (side_inter(i).ne.1) cycle
            if (sqrt( sum( ( ipnts(:,i)-ipnts(:,ind1) )**2.0 ) ).gt.0.001) then
                ind2=i
                nunique = 2
                exit
            endif
        enddo
  
        if(nunique.eq.2) then
            vi = ipnts(:,ind2) - ipnts(:,ind1)
            if (dot_product(v0,vi).gt.0.0) then
                r_enter = ipnts(:,ind1)
                r_exit  = ipnts(:,ind2)
            else
                r_enter = ipnts(:,ind2)
                r_exit  = ipnts(:,ind1)
            endif
            length = sqrt(sum((r_exit - r_enter)**2.0))
        endif
    endif
  
end subroutine grid_intersect

subroutine get_indices(pos, ind)
    real(double),  dimension(3), intent(in)  :: pos
    integer(long), dimension(3), intent(out) :: ind
  
    real(double),  dimension(3) :: mini
    integer(long), dimension(3) :: maxind
    integer :: i
  
    maxind(1) = beam_grid%nx
    maxind(2) = beam_grid%ny
    maxind(3) = beam_grid%nz
  
    mini(1) = minval(beam_grid%xc) - 0.5*beam_grid%dr(1)
    mini(2) = minval(beam_grid%yc) - 0.5*beam_grid%dr(2)
    mini(3) = minval(beam_grid%zc) - 0.5*beam_grid%dr(3)
  
    do i=1,3
        ind(i) = floor((pos(i)-mini(i))/beam_grid%dr(i)) + 1
        if (ind(i).gt.maxind(i)) ind(i)=maxind(i)
        if (ind(i).lt.1) ind(i)=1
    enddo
  
end subroutine get_indices

subroutine get_position(ind, pos)
    integer(long), dimension(3), intent(in)  :: ind
    real(double), dimension(3), intent(out)  :: pos
  
    pos(1) = beam_grid%xc(ind(1))
    pos(2) = beam_grid%yc(ind(2))
    pos(3) = beam_grid%zc(ind(3))

end subroutine get_position

subroutine track(rin, vin, tracks, ncell, los_intersect)
    !!track computes the path of a neutral through the beam grid
    real(double), dimension(3), intent(in)           :: rin  ! initial position
    real(double), dimension(3), intent(in)           :: vin  ! velocitiy
    type(ParticleTrack), dimension(:), intent(inout) :: tracks
    integer(long), intent(out)                       :: ncell
    logical, intent(out), optional                   :: los_intersect
  
    integer :: cc, i, ii, mind
    integer, dimension(3) :: ind
    logical :: in_plasma1, in_plasma2, in_plasma_tmp, los_inter
    real(double) :: dT, dt1, inv_50
    real(double), dimension(3) :: dt_arr, dr
    real(double), dimension(3) :: vn, inv_vn
    real(double), dimension(3) :: ri, ri_tmp, ri_cell 
    integer, dimension(3) :: sgn
    integer, dimension(3) :: gdims
    integer, dimension(1) :: minpos
  
    vn = vin ;  ri = rin ; sgn = 0 ; ncell = 0
  
    if(dot_product(vin,vin).eq.0.0) then
        return
    endif
    
    gdims(1) = beam_grid%nx
    gdims(2) = beam_grid%ny
    gdims(3) = beam_grid%nz
  
    !! define actual cell
    call get_indices(ri,ind)
    ri_cell = [beam_grid%xc(ind(1)), &
               beam_grid%yc(ind(2)), &
               beam_grid%zc(ind(3))]
  
    do i=1,3
        if (vn(i).gt.0.0) sgn(i) = 1
        if (vn(i).lt.0.0) sgn(i) =-1
        if (vn(i).eq.0.0) vn(i)  = 1.0d-3
    enddo
  
    dr = beam_grid%dr*sgn
    inv_vn = 1/vn
    inv_50 = 1.0/50.0
    cc=1
    los_inter = .False.
    tracks%time = 0.d0
    tracks%flux = 0.d0
    call in_plasma(ri,in_plasma1)
    track_loop: do i=1,beam_grid%ntrack
        if(cc.gt.beam_grid%ntrack) exit track_loop
  
        if(spec_chords%los_inter(ind(1),ind(2),ind(3)).and.(.not.los_inter))then
            los_inter = .True.
        endif
        dt_arr = abs(( (ri_cell + 0.5*dr) - ri)*inv_vn)
        minpos = minloc(dt_arr)
        mind = minpos(1)
        dT = dt_arr(mind)
        ri_tmp = ri + dT*vn
        call in_plasma(ri_tmp,in_plasma2)
        if(in_plasma1.neqv.in_plasma2) then
            dt1 = 0.0
            track_fine: do ii=1,50
                dt1 = dt1 + dT*inv_50
                ri_tmp = ri + vn*dt1
                call in_plasma(ri_tmp,in_plasma_tmp)
                if(in_plasma2.eqv.in_plasma_tmp) exit track_fine
            enddo track_fine
            tracks(cc)%pos = ri + 0.5*dt1*vn
            tracks(cc+1)%pos = ri + 0.5*(dt1 + dT)*vn
            tracks(cc)%time = dt1
            tracks(cc+1)%time = dT - dt1
            tracks(cc)%ind = ind
            tracks(cc+1)%ind = ind
            cc = cc + 2
        else
            tracks(cc)%pos = ri + 0.5*dT*vn
            tracks(cc)%time = dT
            tracks(cc)%ind = ind
            cc = cc + 1
        endif
        in_plasma1 = in_plasma2
  
        ri = ri + dT*vn
        ind(mind) = ind(mind) + sgn(mind)
        ri_cell(mind) = ri_cell(mind) + dr(mind)
  
        if (ind(mind).gt.gdims(mind)) exit track_loop
        if (ind(mind).lt.1) exit track_loop
    enddo track_loop
    ncell = cc-1
    if(present(los_intersect)) then
        los_intersect = los_inter
    endif

end subroutine track

!============================================================================
!---------------------------Interpolation Routines---------------------------
!============================================================================
subroutine interpol1D_coeff(xmin,dx,nx,xout,i,b1,b2,err)
    !!Linear interpolation coefficients and index for a 1D grid y(x)
    real(double), intent(in)       :: xmin
    real(double), intent(in)       :: dx
    integer, intent(in)            :: nx
    real(double), intent(in)       :: xout
    integer, intent(out)           :: i
    real(double), intent(out)      :: b1
    real(double), intent(out)      :: b2
    integer, intent(out), optional :: err
  
    real(double) :: x1, xp
    integer :: err_status
  
    err_status = 1
    xp = max(xout,xmin)
    i = floor((xp - xmin)/dx)+1
  
    if ((i.gt.0).and.(i.le.(nx-1))) then
        x1 = xmin + (i-1)*dx
  
        b2 = (xp - x1)/dx
        b1 = (1.0 - b2)
  
        err_status = 0
    else
        b1 = 0.d0 ; b2 = 0.d0
    endif
  
    if(present(err)) err = err_status

end subroutine interpol1D_coeff

subroutine interpol1D_coeff_arr(x,xout,i,b1,b2,err)
    !!Linear interpolation coefficients and index for a 1D grid y(x)
    real(double), dimension(:), intent(in) :: x
    real(double), intent(in)               :: xout
    integer, intent(out)                   :: i
    real(double), intent(out)              :: b1
    real(double), intent(out)              :: b2
    integer, intent(out), optional         :: err
  
    real(double) :: xmin, dx
    integer :: sx,err_status
  
    err_status = 1
    sx = size(x)
    xmin = minval(x)
    dx = abs(x(2)-x(1))
  
    call interpol1D_coeff(xmin, dx, sx, xout, &
                          i, b1, b2, err_status)
  
    if(present(err)) err = err_status

end subroutine interpol1D_coeff_arr

subroutine interpol2D_coeff(xmin,dx,nx,ymin,dy,ny,xout,yout,i,j,b11,b12,b21,b22,err)
    !!Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(double), intent(in)       :: xmin
    real(double), intent(in)       :: dx
    integer, intent(in)            :: nx
    real(double), intent(in)       :: ymin
    real(double), intent(in)       :: dy
    integer, intent(in)            :: ny
    real(double), intent(in)       :: xout
    real(double), intent(in)       :: yout
    integer, intent(out)           :: i
    integer, intent(out)           :: j
    real(double), intent(out)      :: b11
    real(double), intent(out)      :: b12
    real(double), intent(out)      :: b21
    real(double), intent(out)      :: b22
    integer, intent(out), optional :: err
  
    real(double) :: x1, x2, y1, y2, xp, yp
    integer :: err_status
  
    err_status = 1
    xp = max(xout,xmin)
    yp = max(yout,ymin)
    i = floor((xp-xmin)/dx)+1
    j = floor((yp-ymin)/dy)+1
  
    if (((i.gt.0).and.(i.le.(nx-1))).and.((j.gt.0).and.(j.le.(ny-1)))) then
        x1 = xmin + (i-1)*dx
        x2 = x1 + dx
        y1 = ymin + (j-1)*dy
        y2 = y1 + dy
  
        b11 = ((x2 - xp) * (y2 - yp))/(dx*dy)
        b21 = ((xp - x1) * (y2 - yp))/(dx*dy)
        b12 = ((x2 - xp) * (yp - y1))/(dx*dy)
        b22 = ((xp - x1) * (yp - y1))/(dx*dy)
  
        err_status = 0
    else
        b11 = 0.d0 ; b12 = 0.d0 ; b21 = 0.d0 ; b22 = 0.d0
    endif
  
    if(present(err)) err = err_status

end subroutine interpol2D_coeff

subroutine interpol2D_coeff_arr(x,y,xout,yout,i,j,b11,b12,b21,b22,err)
    !!Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(double), dimension(:), intent(in) :: x
    real(double), dimension(:), intent(in) :: y
    real(double), intent(in)               :: xout
    real(double), intent(in)               :: yout
    integer, intent(out)                   :: i
    integer, intent(out)                   :: j
    real(double), intent(out)              :: b11
    real(double), intent(out)              :: b12
    real(double), intent(out)              :: b21
    real(double), intent(out)              :: b22
    integer, intent(out), optional         :: err
  
    real(double) :: xmin, ymin, dx, dy
    integer :: sx, sy, err_status
  
    err_status = 1
    sx = size(x)
    sy = size(y)
    xmin = minval(x)
    ymin = minval(y)
    dx = abs(x(2)-x(1))
    dy = abs(y(2)-y(1))
  
    call interpol2D_coeff(xmin, dx, sx, ymin, dy, sy, xout, yout, &
                          i, j, b11, b12, b21, b22, err_status)
  
    if(present(err)) err = err_status

end subroutine interpol2D_coeff_arr

subroutine interpol1D_arr(x, y, xout, yout, err)
    !!Interpolate on a uniform 1D grid y(x)
    real(double), dimension(:), intent(in) :: x
    real(double), dimension(:), intent(in) :: y
    real(double), intent(in)               :: xout
    real(double), intent(out)              :: yout
    integer, intent(out), optional         :: err
  
    real(double) :: b1
    real(double) :: b2
    integer :: i, err_status
  
    err_status = 1
    call interpol_coeff(x,xout,i,b1,b2,err_status)
    if(err_status.eq.0) then
        yout = b1*y(i) + b2*y(i+1)
    else
        yout = 0.d0
    endif
  
    if(present(err)) err = err_status
  
end subroutine interpol1D_arr

subroutine interpol2D_arr(x, y, z, xout, yout, zout, err)
    !!Interpolate on a 2D grid z(x,y)
    real(double), dimension(:), intent(in)   :: x
    real(double), dimension(:), intent(in)   :: y
    real(double), dimension(:,:), intent(in) :: z
    real(double), intent(in)                 :: xout
    real(double), intent(in)                 :: yout
    real(double), intent(out)                :: zout
    integer, intent(out), optional           :: err
  
    real(double) :: b11
    real(double) :: b12
    real(double) :: b21
    real(double) :: b22
    integer :: i, j, err_status
  
    err_status = 1
    call interpol_coeff(x,y,xout,yout,i,j,b11,b12,b21,b22,err_status)
    if(err_status.eq.0) then
        zout = b11*z(i,j) + b12*z(i,j+1) + b21*z(i+1,j) + b22*z(i+1,j+1)
    else
        zout = 0.d0
    endif
  
    if(present(err)) err = err_status

end subroutine interpol2D_arr

subroutine interpol2D_2D_arr(x, y, z, xout, yout, zout, err)
    !!Interpolate on a 2D grid of 2D arrays z(:,:,x,y)
    real(double), dimension(:), intent(in)       :: x
    real(double), dimension(:), intent(in)       :: y
    real(double), dimension(:,:,:,:), intent(in) :: z
    real(double), intent(in)                     :: xout
    real(double), intent(in)                     :: yout
    real(double), dimension(:,:), intent(out)    :: zout
    integer, intent(out), optional               :: err
  
    real(double) :: b11
    real(double) :: b12
    real(double) :: b21
    real(double) :: b22
    integer :: i, j, err_status
  
    err_status = 1
    call interpol_coeff(x,y,xout,yout,i,j,b11,b12,b21,b22,err_status)
    if(err_status.eq.0) then
        zout = b11*z(:,:,i,j) + b12*z(:,:,i,j+1) + b21*z(:,:,i+1,j) + b22*z(:,:,i+1,j+1)
    else
        zout = 0.0
    endif
  
    if(present(err)) err = err_status
  
end subroutine interpol2D_2D_arr

!=============================================================================
!-------------------------Profiles and Fields Routines------------------------
!=============================================================================
subroutine in_plasma(xyz, inp)
    real(double), dimension(3), intent(in) :: xyz
    logical, intent(out)                   :: inp
  
    real(double), dimension(3)             :: uvw
    real(double)                           :: R, W, phi, mask
    integer                                :: err
  
    err = 1
    !! Convert to machine coordinates
    call xyz_to_uvw(xyz,uvw)
    R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    W = uvw(3)
    phi = atan2(uvw(2),uvw(1))
  
    !! Interpolate mask value
    call interpol(inter_grid%r, inter_grid%z, equil%mask, R, W, mask, err)
  
    if((mask.ge.0.5).and.(err.eq.0)) then
        inp = .True.
    else
        inp = .False.
    endif

end subroutine in_plasma

subroutine get_plasma(plasma, pos, ind)
    type(LocalProfiles), intent(out)                  :: plasma
    real(double), dimension(3), intent(in), optional  :: pos
    integer(long), dimension(3), intent(in), optional :: ind
  
    logical :: inp
    real(double), dimension(3) :: xyz, uvw, vrot_uvw
    real(double) :: R, W, phi, s, c
    real(double) :: b11, b12, b21, b22
    integer :: i, j
  
    plasma%in_plasma = .False.
  
    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos
  
    call in_plasma(xyz,inp)
    if(inp) then
        !! Convert to machine coordinates
        call xyz_to_uvw(xyz,uvw)
        R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        W = uvw(3)
        phi = atan2(uvw(2),uvw(1))
  
        call interpol_coeff(inter_grid%r, inter_grid%z, R, W, & 
                              i, j, b11, b12, b21, b22)
  
        plasma = b11*equil%plasma(i,j)   + b12*equil%plasma(i,j+1) + &
                 b21*equil%plasma(i+1,j) + b22*equil%plasma(i+1,j+1)
  
        s = sin(phi) ; c = cos(phi)
        vrot_uvw(1) = plasma%vr*c - plasma%vt*s 
        vrot_uvw(2) = plasma%vr*s + plasma%vt*c
        vrot_uvw(3) = plasma%vz
        plasma%vrot = matmul(beam_grid%inv_basis,vrot_uvw)
        plasma%pos = xyz
        plasma%in_plasma = .True.
    endif

end subroutine get_plasma

subroutine calc_perp_vectors(b, a, c)
  !!Returns normalized vectors that are perpendicular to b
  !!such that a x c = bnorm
  real(double), dimension(3), intent(in)  :: b
  real(double), dimension(3), intent(out) :: a
  real(double), dimension(3), intent(out) :: c

  real(double), dimension(3) :: bnorm

  bnorm=b/normp(b)

  if (abs(bnorm(3)).eq.1) then
      a=[1.d0,0.d0,0.d0]
      c=[0.d0,1.d0,0.d0]
  else
      if (bnorm(3).eq.0.) then
          a=[0.d0,0.d0,1.d0]
          c=[bnorm(2),-bnorm(1), 0.d0]/sqrt(bnorm(1)**2+bnorm(2)**2)
      else
          a=[bnorm(2),-bnorm(1),0.d0]/sqrt(bnorm(1)**2+bnorm(2)**2)
          c=-[ a(2) , -a(1) , (a(1)*bnorm(2)-a(2)*bnorm(1))/bnorm(3) ]
          c=c/normp(c)
          if(bnorm(3).lt.0.0) then
              c=-c
          endif
      endif
  endif

end subroutine calc_perp_vectors

subroutine get_fields(fields,pos,ind)
    type(LocalEMFields),intent(out)                   :: fields
    real(double), dimension(3), intent(in), optional  :: pos
    integer(long), dimension(3), intent(in), optional :: ind
  
    logical :: inp
    real(double), dimension(3) :: xyz, uvw
    real(double), dimension(3) :: uvw_bfield, uvw_efield
    real(double), dimension(3) :: xyz_bfield, xyz_efield
    real(double) :: R, W, phi, s, c
    real(double) :: b11, b12, b21, b22
    integer :: i, j
  
    fields%in_plasma = .False.
  
    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos
    
    call in_plasma(xyz,inp)
    if(inp) then
        !! Convert to machine coordinates
        call xyz_to_uvw(xyz,uvw)
        R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
        W = uvw(3)
        phi = atan2(uvw(2),uvw(1))
        s = sin(phi) ; c = cos(phi)
  
        call interpol_coeff(inter_grid%r, inter_grid%z, R, W, & 
                              i, j, b11, b12, b21, b22)
  
        fields = b11*equil%fields(i,j) + b12*equil%fields(i,j+1) + &
                 b21*equil%fields(i+1,j) + b22*equil%fields(i+1,j+1)
  
        !Convert cylindrical coordinates to uvw
        uvw_bfield(1) = c*fields%br - s*fields%bt
        uvw_bfield(2) = s*fields%br + c*fields%bt
        uvw_bfield(3) = fields%bz
        uvw_efield(1) = c*fields%er - s*fields%et
        uvw_efield(2) = s*fields%er + c*fields%et
        uvw_efield(3) = fields%ez
  
        !Represent fields in beam grid coordinates
        xyz_bfield = matmul(beam_grid%inv_basis,uvw_bfield)
        xyz_efield = matmul(beam_grid%inv_basis,uvw_efield)
  
        !Calculate field directions and magnitudes
        fields%b_abs = normp(xyz_bfield)
        fields%e_abs = normp(xyz_efield)
        if(fields%b_abs.ne.0) fields%b_norm = xyz_bfield/fields%b_abs
        if(fields%e_abs.ne.0) fields%e_norm = xyz_efield/fields%e_abs
  
        call calc_perp_vectors(fields%b_norm,fields%a_norm,fields%c_norm)
  
        fields%pos = xyz
        fields%in_plasma = .True.
    endif

end subroutine get_fields

subroutine get_distribution(fbeam, pos, ind)
    real(double), dimension(:,:), intent(out)        :: fbeam
    real(double), dimension(3), intent(in), optional :: pos
    integer(long), dimension(3),intent(in), optional :: ind

    real(double), dimension(3) :: xyz, uvw
    real(double) :: R, Z
    logical :: in_plasma1
    integer :: err
  
    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos
  
    !! Convert to machine coordinates
    call xyz_to_uvw(xyz,uvw)
    R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    Z = uvw(3)
  
    in_plasma1= .False.
    call in_plasma(xyz,in_plasma1)
    if(in_plasma1) then
        call interpol(inter_grid%r, inter_grid%z, fbm%f, R, Z, fbeam,err)
    else
        fbeam = 0.0
    endif

end subroutine get_distribution

subroutine get_ep_denf(energy, pitch, denf, pos,ind)
    real(double), intent(in)                          :: energy
    real(double), intent(in)                          :: pitch
    real(double), intent(out)                         :: denf
    real(double), dimension(3), intent(in), optional  :: pos
    integer(long), dimension(3), intent(in), optional :: ind
  
    real(double), dimension(3) :: xyz, uvw
    real(double), dimension(fbm%nenergy,fbm%npitch)  :: fbeam 
    integer(long), dimension(2) :: epi
    integer(long), dimension(1) :: dummy
    real(double) :: R, Z
    real(double) :: dE, dp
    logical :: in_plasma1
  
    if(present(ind)) call get_position(ind,xyz)
    if(present(pos)) xyz = pos
  
    !! Convert to machine coordinates
    call xyz_to_uvw(xyz,uvw)
    R = sqrt(uvw(1)*uvw(1) + uvw(2)*uvw(2))
    Z = uvw(3)
  
    dummy = minloc(abs(fbm%energy - energy))
    epi(1) = dummy(1)
    dummy = minloc(abs(fbm%pitch - pitch))
    epi(2) = dummy(1)
    dE = abs(fbm%energy(epi(1)) - energy)
    dp = abs(fbm%pitch(epi(2)) - pitch)
  
    in_plasma1= .False.
    call in_plasma(xyz,in_plasma1)
    if(in_plasma1.and.(dE.le.fbm%dE).and.(dp.le.fbm%dp)) then
        call interpol(inter_grid%r, inter_grid%z, fbm%f, R, Z, fbeam)
        denf = fbeam(epi(1),epi(2))
    else
        denf = 0.0
    endif

end subroutine get_ep_denf

!=============================================================================
!--------------------------Result Storage Routines----------------------------
!=============================================================================
subroutine store_neutrals(ind, neut_type, dens, store_iter)
    integer(long), dimension(3), intent(in) :: ind
    integer, intent(in)                     :: neut_type
    real(double), dimension(:), intent(in)  :: dens
    logical, intent(in), optional           :: store_iter
  
    logical :: iter
  
    if(present(store_iter)) then
        iter = store_iter
    else
        iter = .False.
    endif
  
    !$OMP CRITICAL(store_neutrals_1)
    if(iter) halo_iter_dens(neut_type) = halo_iter_dens(neut_type) + sum(dens)
      neut%dens(:,neut_type,ind(1),ind(2),ind(3)) = &
      neut%dens(:,neut_type,ind(1),ind(2),ind(3))+dens ![neutrals/cm^3]
    !$OMP END CRITICAL(store_neutrals_1)
end subroutine store_neutrals

subroutine store_births(ind, neut_type, dflux)
    integer(long), dimension(3), intent(in) :: ind
    integer(long), intent(in)               :: neut_type
    real(double), intent(in)                :: dflux
  
    !$OMP CRITICAL(store_births_1)
    birth%dens( neut_type,ind(1),ind(2),ind(3))= &
     birth%dens(neut_type,ind(1),ind(2),ind(3)) + dflux
    !$OMP END CRITICAL(store_births_1)
end subroutine store_births

subroutine store_npa(det, ri, rf, vn, flux)
    integer, intent(in)                    :: det !! Detector/Channel Number
    real(double), dimension(3), intent(in) :: ri !! birth position
    real(double), dimension(3), intent(in) :: rf !! detector position
    real(double), dimension(3), intent(in) :: vn !! particle velocity
    real(double), intent(in)               :: flux !sum(states)/(nlaunch*nloop)*beam_grid%dv [neutrals/s]a
  
    type(LocalEMFields) :: fields
    real(double), dimension(3) :: uvw_ri, uvw_rf,vn_norm
    real(double) :: energy, pitch
    integer(long), dimension(1) :: ienergy
    type(NPAParticle), dimension(:), allocatable :: parts
  
    ! Convert to machine coordinates
    call xyz_to_uvw(ri,uvw_ri)
    call xyz_to_uvw(rf,uvw_rf)
  
    ! Calculate energy 
    energy = inputs%ab*v2_to_E_per_amu*dot_product(vn,vn)
  
    ! Calculate pitch if distribution actually uses pitch
    if(inputs%dist_type.le.2) then
        call get_fields(fields, pos = ri)
        vn_norm = vn/normp(vn)
        pitch = dot_product(fields%b_norm,vn_norm)
    else
        pitch = 0.d0
    endif
  
    !$OMP CRITICAL(store_npa_1)
    npa%npart = npa%npart + 1
    if(npa%npart.gt.npa%nmax) then
        allocate(parts(npa%npart-1))
        parts = npa%part
        deallocate(npa%part)
        npa%nmax = int(npa%nmax*2)
        allocate(npa%part(npa%nmax))
        npa%part(1:(npa%npart-1)) = parts
        deallocate(parts)
    endif
    npa%part(npa%npart)%detector = det
    npa%part(npa%npart)%xi = uvw_ri(1)
    npa%part(npa%npart)%yi = uvw_ri(2)
    npa%part(npa%npart)%zi = uvw_ri(3)
    npa%part(npa%npart)%xf = uvw_rf(1)
    npa%part(npa%npart)%yf = uvw_rf(2)
    npa%part(npa%npart)%zf = uvw_rf(3)
    npa%part(npa%npart)%energy = energy
    npa%part(npa%npart)%pitch = pitch
    npa%part(npa%npart)%weight = flux
    ienergy = minloc(abs(npa%energy - energy))
    npa%flux(ienergy(1),det) = npa%flux(ienergy(1),det) + flux/fbm%dE
    !$OMP END CRITICAL(store_npa_1)

end subroutine store_npa

!=============================================================================
!--------------------------Atomic Physics Routines----------------------------
!=============================================================================
subroutine neut_rates(denn, vi, vn, rates)
    real(double), dimension(nlevs), intent(in)  :: denn  !!density of neutrals cm-3
    real(double), dimension(3),     intent(in)  :: vi,vn !!of neutrals/ions (cm/s)
    real(double), dimension(nlevs), intent(out) :: rates !! rates
  
    real(double), dimension(nlevs,nlevs) :: neut  !!rate coeff
    real(double) :: eb !! relative Energy
    real(double) :: b1, b2, dlogE, logEmin, logeb
    real(double) :: vrel !! relative velocity
    integer :: ebi, neb, err
  
    !Eeff
    vrel=normp(vi-vn)
    eb=v2_to_E_per_amu*vrel**2  ! [kev/amu]
    logeb = log10(eb)
    logEmin = tables%H_H_cx%logemin
    dlogE = tables%H_H_cx%dlogE
    neb = tables%H_H_cx%nenergy
    call interpol_coeff(logEmin,dlogE,neb,logeb,ebi,b1,b2,err)
    if(err.eq.1) then
        write(*,'(a)') "NEUT_RATES: Eb out of range of H_H_cx table"
        stop
    endif
  
    neut(:,:) = (b1*tables%H_H_cx%log_cross(:,:,ebi) + &
                 b2*tables%H_H_cx%log_cross(:,:,ebi+1))
  
    where (neut.lt.tables%H_H_cx%minlog_cross)
        neut = 0.d0
    elsewhere
        neut = 10.d0**neut
    end where
  
    rates=matmul(neut,denn)*vrel

end subroutine neut_rates

subroutine get_beam_cx_prob(ind, pos, v_ion, types, prob)
    integer(long), dimension(3), intent(in)     :: ind
    real(double), dimension(3), intent(in)      :: pos
    real(double), dimension(3), intent(in)      :: v_ion
    integer(long), dimension(:),intent(in)      :: types
    real(double), dimension(nlevs), intent(out) :: prob
  
    integer :: ntypes, i, ii
    real(double), dimension(nlevs) :: rates, denn
    real(double), dimension(3) :: vhalo, vnbi ,vn
  
    vnbi = pos - nbi%src
    vnbi = vnbi/normp(vnbi)*nbi%vinj
  
    ntypes = size(types)
    prob = 0
  
    do i=1,ntypes
        if((types(i).le.3).and.(types(i).ne.0)) then
            ! CX with full type'th energy NBI neutrals
            denn = neut%dens(:,types(i),ind(1),ind(2),ind(3))
            vn = vnbi/sqrt(real(types(i)))
            call neut_rates(denn,v_ion,vn,rates)
            prob = prob+rates
        else
            denn = neut%dens(:,types(i),ind(1),ind(2),ind(3))
            do ii=1,int(n_halo_neutrate)
                call mc_halo(ind, vhalo)
                call neut_rates(denn,v_ion,vhalo,rates)
                prob=prob+rates/n_halo_neutrate
            enddo
        endif
    enddo

end subroutine get_beam_cx_prob

subroutine get_rate_matrix(plasma, i_type, eb, rmat)
    type(LocalProfiles), intent(in)                   :: plasma
    integer, intent(in)                               :: i_type
    real(double), intent(in)                          :: eb
    real(double), dimension(nlevs,nlevs), intent(out) :: rmat
  
    real(double) :: logEmin, dlogE, logeb
    real(double) :: logTmin, dlogT, logti, logte
    integer :: neb, nt
    real(double) :: b11, b12, b21, b22, dene, denp, denimp
    real(double), dimension(nlevs,nlevs) :: H_H_pop, H_e_pop, H_Aq_pop
    real(double), dimension(nlevs) :: H_H_depop, H_e_depop, H_Aq_depop
    integer :: ebi, tii, tei, n, err_status
  
    H_H_pop = 0.d0
    H_e_pop = 0.d0
    H_Aq_pop = 0.d0
    H_H_depop = 0.d0
    H_e_depop = 0.d0
    H_Aq_depop = 0.d0
    denp = plasma%denp
    dene = plasma%dene
    denimp = plasma%denimp
    logeb = log10(eb)
    logti = log10(plasma%ti)
    logte = log10(plasma%te)
    !!H_H
    err_status = 1
    logEmin = tables%H_H%logemin
    logTmin = tables%H_H%logtmin
    dlogE = tables%H_H%dlogE
    dlogT = tables%H_H%dlogT
    neb = tables%H_H%nenergy
    nt = tables%H_H%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logti, &
                        ebi, tii, b11, b12, b21, b22, err_status)
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_H table. Setting H_H rates to zero"
        write(*,'("eb = ",f6.3," [keV]")') eb
        write(*,'("ti = ",f6.3," [keV]")') plasma%ti
        denp = 0.d0
    endif
  
    H_H_pop = (b11*tables%H_H%log_pop(:,:,ebi,tii,i_type)   + &
               b12*tables%H_H%log_pop(:,:,ebi,tii+1,i_type) + &
               b21*tables%H_H%log_pop(:,:,ebi+1,tii,i_type) + &
               b22*tables%H_H%log_pop(:,:,ebi+1,tii+1,i_type))
    where (H_H_pop.lt.tables%H_H%minlog_pop)
        H_H_pop = 0.d0
    elsewhere
        H_H_pop = denp * 10.d0**H_H_pop
    end where
  
    H_H_depop = (b11*tables%H_H%log_depop(:,ebi,tii,i_type)   + &
                 b12*tables%H_H%log_depop(:,ebi,tii+1,i_type) + &
                 b21*tables%H_H%log_depop(:,ebi+1,tii,i_type) + &
                 b22*tables%H_H%log_depop(:,ebi+1,tii+1,i_type))
    where (H_H_depop.lt.tables%H_H%minlog_depop)
        H_H_depop = 0.d0
    elsewhere
        H_H_depop = denp * 10.d0**H_H_depop
    end where
  
    !!H_e
    err_status = 1
    logEmin = tables%H_e%logemin
    logTmin = tables%H_e%logtmin
    dlogE = tables%H_e%dlogE
    dlogT = tables%H_e%dlogT
    neb = tables%H_e%nenergy
    nt = tables%H_e%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logte, &
                        ebi, tei, b11, b12, b21, b22, err_status)
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Te out of range of H_e table. Setting H_e rates to zero"
        write(*,'("eb = ",f6.3," [keV]")') eb
        write(*,'("te = ",f6.3," [keV]")') plasma%te
        dene = 0.d0
    endif
  
    H_e_pop = (b11*tables%H_e%log_pop(:,:,ebi,tei,i_type)   + &
               b12*tables%H_e%log_pop(:,:,ebi,tei+1,i_type) + &
               b21*tables%H_e%log_pop(:,:,ebi+1,tei,i_type) + &
               b22*tables%H_e%log_pop(:,:,ebi+1,tei+1,i_type))
    where (H_e_pop.lt.tables%H_e%minlog_pop)
        H_e_pop = 0.d0
    elsewhere
        H_e_pop = dene * 10.d0**H_e_pop
    end where
  
    H_e_depop = (b11*tables%H_e%log_depop(:,ebi,tei,i_type)   + &
                 b12*tables%H_e%log_depop(:,ebi,tei+1,i_type) + &
                 b21*tables%H_e%log_depop(:,ebi+1,tei,i_type) + &
                 b22*tables%H_e%log_depop(:,ebi+1,tei+1,i_type))
                         
    where (H_e_depop.lt.tables%H_e%minlog_depop)
        H_e_depop = 0.d0
    elsewhere
        H_e_depop = dene * 10.d0**H_e_depop
    end where
  
    !!H_Aq
    err_status = 1
    logEmin = tables%H_Aq%logemin
    logTmin = tables%H_Aq%logtmin
    dlogE = tables%H_Aq%dlogE
    dlogT = tables%H_Aq%dlogT
    neb = tables%H_Aq%nenergy
    nt = tables%H_Aq%ntemp
    call interpol_coeff(logEmin, dlogE, neb, logTmin, dlogT, nt, &
                        logeb, logti, &
                        ebi, tii, b11, b12, b21, b22, err_status)
    if(err_status.eq.1) then
        write(*,'(a)') "GET_RATE_MATRIX: Eb or Ti out of range of H_Aq table. Setting H_Aq rates to zero"
        write(*,'("eb = ",f6.3," [keV]")') eb
        write(*,'("ti = ",f6.3," [keV]")') plasma%ti
        denimp = 0.d0
    endif
  
    H_Aq_pop = (b11*tables%H_Aq%log_pop(:,:,ebi,tii,i_type)   + &
                b12*tables%H_Aq%log_pop(:,:,ebi,tii+1,i_type) + &
                b21*tables%H_Aq%log_pop(:,:,ebi+1,tii,i_type) + &
                b22*tables%H_Aq%log_pop(:,:,ebi+1,tii+1,i_type))
    where (H_Aq_pop.lt.tables%H_Aq%minlog_pop)
        H_Aq_pop = 0.d0
    elsewhere
        H_Aq_pop = denimp * 10.d0**H_Aq_pop
    end where
    H_Aq_depop = (b11*tables%H_Aq%log_depop(:,ebi,tii,i_type)   + &
                  b12*tables%H_Aq%log_depop(:,ebi,tii+1,i_type) + &
                  b21*tables%H_Aq%log_depop(:,ebi+1,tii,i_type) + &
                  b22*tables%H_Aq%log_depop(:,ebi+1,tii+1,i_type))
  
    where (H_Aq_depop.lt.tables%H_Aq%minlog_depop)
        H_Aq_depop = 0.d0
    elsewhere
        H_Aq_depop = denimp * 10.d0**H_Aq_depop
    end where
  
    rmat = tables%einstein + H_H_pop + H_e_pop + H_Aq_pop
    do n=1,nlevs
        rmat(n,n) = -sum(tables%einstein(:,n)) - H_H_depop(n) - H_e_depop(n) - H_Aq_depop(n)
    enddo

end subroutine get_rate_matrix

subroutine colrad(plasma,i_type,vn,dt,states,dens,photons)
    type(LocalProfiles), intent(in)             :: plasma !! Local plasma parameters
    integer, intent(in)                         :: i_type !! Ion type (beam,thermal)
    real(double), dimension(:), intent(in)      :: vn  !!velocitiy (cm/s)
    real(double), intent(in)                    :: dt  !!time interval in cell
    real(double), dimension(:), intent(inout)   :: states  !! Density of states
    real(double), dimension(nlevs), intent(out) :: dens !! Density of neutrals
    real(double), intent(out)                   :: photons !! Emitted photons(3->2)
  
    real(double), dimension(nlevs,nlevs) :: matrix  !! Matrix
    real(double) :: b_amu
    real(double) :: vnet_square    !! netto velocity of neutrals
    real(double) :: eb             !! Energy of the fast neutral

    real(double), dimension(nlevs,nlevs) :: eigvec, eigvec_inv
    real(double), dimension(nlevs) :: eigval, coef
    real(double), dimension(nlevs) :: exp_eigval_dt
    real(double) :: iflux !!Initial total flux
    integer :: n

    photons=0.d0
    dens=0.d0
  
    iflux=sum(states)
    if(iflux.lt.colrad_threshold .and. inputs%calc_npa.eq.0)then
        return
    endif
  
    if(.not.plasma%in_plasma) then
        dens = states*dt
        return
    endif
  
    if(i_type.eq.beam_ion) then
        b_amu = inputs%ab
    else
        b_amu = inputs%ai
    endif
    vnet_square=dot_product(vn-plasma%vrot,vn-plasma%vrot)          ![cm/s]
    eb = v2_to_E_per_amu*b_amu*vnet_square ![kev]
    call get_rate_matrix(plasma, i_type, eb, matrix)
  
    call eigen(nlevs,matrix, eigvec, eigval)
    call matinv(eigvec, eigvec_inv)
    coef = matmul(eigvec_inv, states)!coeffs determined from states at t=0
    exp_eigval_dt = exp(eigval*dt)   ! to improve speed (used twice)
    do n=1,nlevs
        if(eigval(n).eq.0.0) eigval(n)=eigval(n)+1 !protect against dividing by zero
    enddo
  
    states = matmul(eigvec, coef * exp_eigval_dt)  ![neutrals/cm^3/s]!
    dens   = matmul(eigvec,coef*(exp_eigval_dt-1.d0)/eigval)
  
    if ((minval(states).lt.0).or.(minval(dens).lt.0)) then
        do n=1,nlevs
            if(states(n).lt.0) states(n)=0.d0
            if(dens(n).lt.0) dens(n)=0.d0
        enddo
    endif
  
    photons=dens(3)*tables%einstein(2,3) !! - [Ph/(s*cm^3)] - !!
  
end subroutine colrad

subroutine attenuate(ri, rf, vi, states, dstep_in)
    real(double), dimension(3), intent(in)        :: ri
    real(double), dimension(3), intent(in)        :: rf
    real(double), dimension(3), intent(in)        :: vi
    real(double), dimension(nlevs), intent(inout) :: states
    real(double), intent(in), optional            :: dstep_in
  
    type(LocalProfiles) :: plasma
    real(double) :: photons, vabs, dt, dstep, dis,max_dis
    real(double), dimension(3) :: r0
    real(double), dimension(nlevs) :: dens
  
    if(present(dstep_in)) then
        dstep=dstep_in
    else
        dstep = sqrt(inter_grid%da) !cm
    endif
  
    max_dis = normp(rf-ri)
  
    vabs = normp(vi)
    dt = dstep/vabs
  
    call get_plasma(plasma,pos=ri)
    r0 = ri
    dis = 0.d0
    do while (dis.le.max_dis)
        call colrad(plasma,beam_ion,vi,dt,states,dens,photons)
        r0 = r0 + vi*dt
        dis = dis + dstep
        call get_plasma(plasma,pos=r0)
    enddo

end subroutine attenuate

subroutine spectrum(vecp, vi, fields, sigma_pi, photons, dlength, lambda, intensity)
    real(double), dimension(3), intent(in)        :: vecp
    real(double), dimension(3), intent(in)        :: vi
    type(LocalEMFields), intent(in)               :: fields
    real(double), intent(in)                      :: sigma_pi
    real(double), intent(in)                      :: photons
    real(double), intent(in)                      :: dlength
    real(double), dimension(n_stark), intent(out) :: lambda
    real(double), dimension(n_stark), intent(out) :: intensity
  
    real(double), dimension(3) :: vp, vn
    real(double), dimension(3) :: bfield, efield
    real(double) :: E, cos_los_Efield, lambda_shifted
    integer, parameter, dimension(n_stark) ::stark_sign= +1*stark_sigma -1*stark_pi
  
    !! vector directing towards the optical head
    vp=vecp/normp(vecp)
  
    ! Calculate Doppler shift
    vn=vi*0.01d0 ! [m/s]
    lambda_shifted = lambda0*(1.d0 + dot_product(vn,vp)/c0)
  
    !! Calculate Stark Splitting
    ! Calculate E-field
    bfield = fields%b_norm*fields%b_abs
    efield = fields%e_norm*fields%e_abs
    efield(1) = efield(1) +  vn(2)*bfield(3) - vn(3)*bfield(2)
    efield(2) = efield(2) - (vn(1)*bfield(3) - vn(3)*bfield(1))
    efield(3) = efield(3) +  vn(1)*bfield(2) - vn(2)*bfield(1)
    E = normp(efield)
  
    !Stark Splitting
    lambda =  lambda_shifted + E * stark_wavel ![A]
  
    !Intensities of stark components
    if (E .eq. 0.d0) then
        cos_los_Efield = 0.d0
    else
        cos_los_Efield = dot_product(vp,efield) / E
    endif
  
    intensity = stark_intens*(1.d0+ stark_sign* cos_los_Efield**2.d0)
    !! --- E.g. mirrors may change the pi to sigma intensity ratio  --- !!
    where (stark_sigma .eq. 1)
        intensity = intensity * sigma_pi
    endwhere
  
    !! --- normalize and multiply with photon density from colrad --- !!
    intensity = intensity/sum(intensity)*photons*dlength

endsubroutine spectrum

subroutine store_bes_photons(pos, vi, photons, neut_type)
    real(double), dimension(3), intent(in) :: pos !! Position of neutral
    real(double), dimension(3), intent(in) :: vi !!velocitiy of neutral [cm/s]
    real(double), intent(in)               :: photons !! photons from colrad [Ph/(s*cm^3)]
    integer,intent(in)                     :: neut_type
  
    real(double), dimension(n_stark) :: lambda, intensity
    real(double) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(long), dimension(3) :: ind
    real(double), dimension(3) :: vp
    integer :: ichan,i,bin
  
    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)
  
    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call spectrum(vp,vi,fields,sigma_pi,photons, &
                      dlength,lambda,intensity)
  
        loop_over_stark: do i=1,n_stark
            bin=floor((lambda(i)-inputs%lambdamin)/inputs%dlambda) + 1
            if (bin.lt.1) cycle loop_over_stark
            if (bin.gt.inputs%nlambda) cycle loop_over_stark
            !$OMP CRITICAL(bes_spectrum)
            spec%bes(bin,ichan,neut_type)= &
              spec%bes(bin,ichan,neut_type) + intensity(i)
            !$OMP END CRITICAL(bes_spectrum)
        enddo loop_over_stark
    enddo loop_over_channels
  
end subroutine store_bes_photons

subroutine store_fida_photons(pos, vi, photons, orbit_class)
    real(double), dimension(3), intent(in) :: pos !! Position of neutral
    real(double), dimension(3), intent(in) :: vi !!velocitiy of fast neutral [cm/s]
    real(double), intent(in)               :: photons !! photons from colrad [Ph/(s*cm^3)]
    integer, intent(in), optional          :: orbit_class
  
    real(double), dimension(n_stark) :: lambda, intensity
    real(double) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(long), dimension(3) :: ind
    real(double), dimension(3) :: vp
    integer :: ichan, i, bin, iclass
  
    if(present(orbit_class)) then
        iclass = orbit_class
    else
        iclass = 1
    endif
  
    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)
  
    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call spectrum(vp,vi,fields,sigma_pi,photons, &
                      dlength,lambda,intensity)
  
        loop_over_stark: do i=1,n_stark
            bin=floor((lambda(i)-inputs%lambdamin)/inputs%dlambda) + 1
            if (bin.lt.1) cycle loop_over_stark
            if (bin.gt.inputs%nlambda) cycle loop_over_stark
            !$OMP CRITICAL(fida_spectrum)
            spec%fida(bin,ichan,iclass)= &
              spec%fida(bin,ichan,iclass) + intensity(i)
            !$OMP END CRITICAL(fida_spectrum)
        enddo loop_over_stark
    enddo loop_over_channels
  
end subroutine store_fida_photons

subroutine store_fw_photons_at_chan(ichan,eind,pind,vp,vi,fields,dlength,sigma_pi,denf,photons)
    integer, intent(in)                    :: ichan
    integer, intent(in)                    :: eind
    integer, intent(in)                    :: pind
    real(double), dimension(3), intent(in) :: vp !! Photon direction
    real(double), dimension(3), intent(in) :: vi !!velocity of neutral [cm/s]
    type(LocalEMFields), intent(in)        :: fields
    real(double), intent(in)               :: dlength
    real(double), intent(in)               :: sigma_pi
    real(double), intent(in)               :: denf
    real(double), intent(in)               :: photons !! photons from colrad [Ph/(s*cm^3)]
  
    real(double), dimension(n_stark) :: lambda,intensity
    real(double) :: dlambda,intens_fac
    integer :: i,bin
  
  
    dlambda=(inputs%lambdamax_wght-inputs%lambdamin_wght)/inputs%nlambda_wght
    intens_fac = (1.d0)/(4.d0*pi*dlambda)
    call spectrum(vp,vi,fields,sigma_pi,photons, &
                  dlength,lambda,intensity)
  
    !$OMP CRITICAL(fida_wght)
    loop_over_stark: do i=1,n_stark
        bin=floor((lambda(i)*0.1 - inputs%lambdamin_wght)/dlambda) + 1
        if (bin.lt.1)                   cycle loop_over_stark
        if (bin.gt.inputs%nlambda_wght) cycle loop_over_stark
        fweight%fida(bin,ichan)= fweight%fida(bin,ichan) + & 
          (denf*intens_fac*1.d4)*intensity(i) !ph/(s*nm*sr*m^2)
        fweight%weight(bin,eind,pind,ichan) = & 
          fweight%weight(bin,eind,pind,ichan) + intensity(i)*intens_fac !(ph*cm)/(s*nm*sr*fast-ion*dE*dp)
    enddo loop_over_stark
    if(denf.gt.0.d0) then
        fweight%mean_f(eind,pind,ichan) = fweight%mean_f(eind,pind,ichan) + &
                                          (denf*intens_fac)*sum(intensity)
    endif
    !$OMP END CRITICAL(fida_wght)
  
end subroutine store_fw_photons_at_chan

subroutine store_fw_photons(eind, pind, pos, vi, denf, photons)
    integer, intent(in)                    :: eind
    integer, intent(in)                    :: pind
    real(double), dimension(3), intent(in) :: pos !! Position of neutral
    real(double), dimension(3), intent(in) :: vi !!velocity of neutral [cm/s]
    real(double), intent(in)               :: denf
    real(double), intent(in)               :: photons !! photons from colrad [Ph/(s*cm^3)]
  
    real(double) :: dlength, sigma_pi
    type(LocalEMFields) :: fields
    integer(long), dimension(3):: ind
    real(double), dimension(3) :: vp
    integer :: ichan
  
    call get_indices(pos,ind)
    call get_fields(fields,pos=pos)
  
    loop_over_channels: do ichan=1,spec_chords%nchan
        dlength = spec_chords%dlength(ichan,ind(1),ind(2),ind(3))
        if(dlength.le.0.0) cycle loop_over_channels
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        vp = pos - spec_chords%los(ichan)%lens
        call store_fw_photons_at_chan(ichan, eind, pind, &
             vp, vi, fields, dlength, sigma_pi, denf, photons)
    enddo loop_over_channels
  
end subroutine store_fw_photons

!=============================================================================
!---------------------------Monte Carlo Routines------------------------------
!=============================================================================
subroutine get_nlaunch(nr_markers,papprox,papprox_tot,nlaunch)
    !! routine to define the number of MC particles started in one cell
    integer(long), intent(in)                   :: nr_markers
    real(double), dimension(:,:,:), intent(in)  :: papprox
    real(double)                  , intent(in)  :: papprox_tot
    real(double), dimension(:,:,:), intent(out) :: nlaunch

    integer  :: i, j, k, cc
    real(double), dimension(:), allocatable :: randomu

    do i=1,1000
       nlaunch(:,:,:)=papprox(:,:,:)/papprox_tot*nr_markers*(1.+i*0.01)
       if(sum(nlaunch).gt.nr_markers) then 
          exit
       endif
    enddo

    allocate(randomu(count(nlaunch.gt.0)))
    call randu(randomu)
    cc=1
    do k = 1, beam_grid%nz
        do j = 1, beam_grid%ny
            do i = 1, beam_grid%nx
                if(nlaunch(i,j,k).gt.0.)then
                    if(mod(nlaunch(i,j,k),1.).gt.randomu(cc))then
                        nlaunch(i,j,k)=nlaunch(i,j,k)+1.
                    endif
                    cc=cc+1
                endif
            enddo
        enddo
    enddo

    do k = 1, beam_grid%nz
        do j = 1, beam_grid%ny
            do i = 1, beam_grid%nx
                nlaunch(i,j,k)=floor(nlaunch(i,j,k))
            enddo
        enddo
    enddo
    deallocate(randomu)

end subroutine get_nlaunch

subroutine pitch_to_vec(pitch, gyroangle, fields, vi_norm)
    real(double), intent(in)                :: pitch
    real(double), intent(in)                :: gyroangle
    type(LocalEMFields), intent(in)         :: fields
    real(double), dimension(3), intent(out) :: vi_norm
  
    real(double) :: sinus
  
    sinus = sqrt(1.d0-pitch**2)
    vi_norm = (sinus*cos(gyroangle)*fields%a_norm + &
               pitch*fields%b_norm + &
               sinus*sin(gyroangle)*fields%c_norm)

end subroutine pitch_to_vec

subroutine gyro_correction(vi, fields, r_gyro)
    real(double), dimension(3), intent(in)  :: vi
    type(LocalEMFields), intent(in)         :: fields
    real(double), dimension(3), intent(out) :: r_gyro

    real(double), dimension(3) :: vxB
    real(double) :: one_over_omega
  
    one_over_omega=inputs%ab*mass_u/(fields%b_abs*e0)
    vxB = cross_product(vi,fields%b_norm)
    r_gyro = vxB*one_over_omega

end subroutine gyro_correction

subroutine mc_fastion(ind,ri,vi,at_guiding_center)
    integer, dimension(3), intent(in)         :: ind !ind of actual cell
    real(double), dimension(3), intent(out)   :: ri !starting position
    real(double), dimension(3), intent(out)   :: vi !velocity [cm/s]
    logical, intent(in), optional             :: at_guiding_center !indicates that ri is at gyrocenter
                                                                   !Defaults to true
    type(LocalEMFields) :: fields
    real(double), dimension(fbm%nenergy,fbm%npitch) :: fbeam
    real(double), dimension(3) :: r_gyro, rp
    real(double) :: eb ,ptch
    integer :: ii, ienergy, ipitch
    real(double) :: vabs, phi
    real(double), dimension(3) :: randomu3
    real(double), dimension(4) :: randomu4
    integer, dimension(1) :: minpos
    integer, dimension(2,1) :: ep_ind
    real(double), dimension(1) :: max_fbm
    logical :: use_inverse_sampler
   
    if(present(at_guiding_center)) then
        use_inverse_sampler = at_guiding_center
    else
        use_inverse_sampler = .True.
    endif
  
    call randu(randomu3)
    ri(1) = beam_grid%xc(ind(1)) + beam_grid%dr(1)*(randomu3(1) - 0.5)
    ri(2) = beam_grid%yc(ind(2)) + beam_grid%dr(2)*(randomu3(2) - 0.5)
    ri(3) = beam_grid%zc(ind(3)) + beam_grid%dr(3)*(randomu3(3) - 0.5)
    vi=0.d0
  
    call get_fields(fields,pos=ri)
    if(.not.fields%in_plasma) return 
  
    if(use_inverse_sampler) then
        call get_distribution(fbeam,pos=ri)
        call randind(fbeam,ep_ind)
        call randu(randomu3)
        eb = fbm%energy(ep_ind(1,1)) + fbm%dE*(randomu3(1)-0.5)
        ptch = fbm%pitch(ep_ind(2,1)) + fbm%dp*(randomu3(2)-0.5)
        phi = 2.d0*pi*randomu3(3)
  
        !! Calculate gyroradius
        vabs  = sqrt(eb/(v2_to_E_per_amu*inputs%ab))
        call pitch_to_vec(ptch,phi,fields,vi)
        vi = vabs*vi
        call gyro_correction(vi,fields,r_gyro)
  
        !! Move a gyro-orbit away
        ri=ri-r_gyro
    else
        !! Use rejection method to determine velocity vector
        rejection_loop: do ii=1,10000
            call randu(randomu4)
            !! Pick a random energy, pitch, and gyro angle
            eb   = fbm%emin + fbm%e_range * randomu4(1)
            ptch = fbm%pmin + fbm%p_range * randomu4(2)
            phi  = 2.d0*pi*randomu4(3)
  
            !! Calculate gyroradius
            vabs  = sqrt(eb/(v2_to_E_per_amu*inputs%ab))
            call pitch_to_vec(ptch,phi,fields,vi)
            vi = vabs*vi
            call gyro_correction(vi,fields,r_gyro)
  
            !! Move a gyro-orbit away and sample distribution there
            rp=ri+r_gyro
  
            call get_distribution(fbeam,pos=rp)
            max_fbm = maxval(fbeam)
            !! Find new cell
            if(max_fbm(1).gt.0.0) then
                !! take point in FBM distribution closest to eb, ptch.
                minpos=minloc(abs(eb - fbm%energy))
                ienergy= minpos(1)
                minpos=minloc(abs(ptch - fbm%pitch ))
                ipitch = minpos(1)
                if((fbeam(ienergy,ipitch)).gt.(randomu4(4)*max_fbm(1)))then
                    return
                endif
            endif
            vi=0.d0
        enddo rejection_loop
  
        write(*,'(a)') 'MC_FASTION: Rejection method found no solution!'
    endif

end subroutine mc_fastion

subroutine mc_halo(ind,vhalo,ri,plasma_in)
    integer, dimension(3), intent(in)                 :: ind    !! index of actual cell
    real(double), dimension(3), intent(out)           :: vhalo !! velocity [cm/s]
    real(double), dimension(3), intent(out), optional :: ri
    type(LocalProfiles), intent(in), optional         :: plasma_in

    type(LocalProfiles) :: plasma
    real(double), dimension(3) :: random3
  
    if(.not.present(plasma_in)) then
        if(present(ri)) then
            call randu(random3)
            ri(1) = beam_grid%xc(ind(1)) + beam_grid%dr(1)*(random3(1) - 0.5)
            ri(2) = beam_grid%yc(ind(2)) + beam_grid%dr(2)*(random3(2) - 0.5)
            ri(3) = beam_grid%zc(ind(3)) + beam_grid%dr(3)*(random3(3) - 0.5)
            call get_plasma(plasma,pos=ri)
        else 
            call get_plasma(plasma,ind=ind)
        endif
    else
        plasma=plasma_in
    endif
  
    call randn(random3)
  
    vhalo = plasma%vrot + sqrt(plasma%ti*0.5/(v2_to_E_per_amu*inputs%ai))*random3 !![cm/s]

end subroutine mc_halo

subroutine mc_nbi(vnbi,efrac,rnbi)
    integer, intent(in)                               :: efrac !! energy fraction
    real(double), dimension(3), intent(out)           :: vnbi  !! velocity [cm/s]
    real(double), dimension(3), intent(out)           :: rnbi  !! position
  
    real(double), dimension(3) :: r_exit
    real(double), dimension(3) :: uvw_src    !! Start position on ion source
    real(double), dimension(3) :: xyz_src    !! Start position on ion source
    real(double), dimension(3) :: uvw_ray    !! NBI velocity in uvw coords
    real(double), dimension(3) :: xyz_ray    !! NBI velocity in xyz coords
    real(double), dimension(2) :: randomu    !! uniform random numbers
    real(double), dimension(2) :: randomn    !! normal random numbers
    real(double) :: length, sqrt_rho, theta
  
    call randu(randomu)
    select case (nbi%shape)
        case (1)
            ! Uniformally sample in rectangle
            xyz_src(1) =  0.d0
            xyz_src(2) =  nbi%widy * 2.d0*(randomu(1)-0.5d0)
            xyz_src(3) =  nbi%widz * 2.d0*(randomu(2)-0.5d0)
        case (2)
            ! Uniformally sample in ellipse
            sqrt_rho = sqrt(randomu(1))
            theta = 2*pi*randomu(2)
            xyz_src(1) = 0.d0
            xyz_src(2) = nbi%widy*sqrt_rho*cos(theta)
            xyz_src(3) = nbi%widz*sqrt_rho*sin(theta)
    end select
  
    !! Create random velocity vector
    call randn(randomn)
    xyz_ray(1)= 1.d0
    xyz_ray(2)=(xyz_src(2)/nbi%focy + tan(nbi%divy(efrac)*randomn(1)))
    xyz_ray(3)=(xyz_src(3)/nbi%focz + tan(nbi%divz(efrac)*randomn(2)))
  
    !! Convert to beam sightline coordinates to beam grid coordinates
    uvw_src = matmul(nbi%basis,xyz_src) + nbi%src
    uvw_ray = matmul(nbi%basis,xyz_ray)
    vnbi = uvw_ray/normp(uvw_ray)
  
    !! ----------- Determine start postition on FIDASIM beam grid --------- !!
    call grid_intersect(uvw_src,vnbi,length,rnbi,r_exit)
    if(length.le.0.0)then
        rnbi=[-1.d0,0.d0,0.d0]
        nbi_outside = nbi_outside + 1
    endif
  
    !! ---- Determine velocity of neutrals corrected by efrac ---- !!
    vnbi = vnbi*nbi%vinj/sqrt(real(efrac))
    
end subroutine mc_nbi

!=============================================================================
!------------------------Primary Simulation Routines--------------------------
!=============================================================================
subroutine ndmc
    integer :: neut_type !! full half third energy
    real(double)  :: nlaunch   !! nr. of markers
    real(double)  :: nneutrals !! # NBI particles
    real(double), dimension(3) :: vnbi      !! velocities(full..)
    real(double), dimension(3) :: rnbi      !! initial position
  
    integer :: jj, ii, kk
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    integer, dimension(3) :: nl_birth
    type(LocalProfiles) :: plasma
    real(double), dimension(nlevs) :: states, dens
    real(double) :: photons, iflux
    integer(long), dimension(3) :: ind
    real(double), dimension(1) :: randomu
    integer, dimension(1) :: randi
  
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') inputs%n_nbi
    endif
  
    !! ------------- calculate nr. of injected neutrals ---------------- !!
    !! # of injected neutrals = NBI power/energy_per_particle
    nneutrals=1.d6*nbi%pinj/ (1.d3*nbi%einj*e0 &
         *( nbi%species_mix(1)      &
         +  nbi%species_mix(2)/2.d0 &
         +  nbi%species_mix(3)/3.d0 ) )
  
    nlaunch=real(inputs%n_nbi)

    !$OMP PARALLEL DO schedule(guided) & 
    !$OMP& private(vnbi,rnbi,tracks,ncell,plasma,nl_birth,randi, &
    !$OMP& states,dens,iflux,photons,neut_type,jj,ii,kk,ind)
    loop_over_markers: do ii=1,inputs%n_nbi
        !! (type = 1: full energy, =2: half energy, =3: third energy
        if(inputs%calc_birth.ge.1) then
            nl_birth = 0
            do kk=1,inputs%n_birth
                call randind(nbi%species_mix,randi)
                nl_birth(randi(1)) = nl_birth(randi(1)) + 1
            enddo
        endif
        energy_fractions: do neut_type=1,3
            call mc_nbi(vnbi,neut_type,rnbi)
            if(rnbi(1).eq.-1)cycle loop_over_markers
  
            call track(rnbi,vnbi,tracks,ncell)
            if(ncell.eq.0) cycle loop_over_markers
  
            !! --------- solve collisional radiative model along track ----- !!
            states=0.d0
            states(1)=nneutrals*nbi%species_mix(neut_type)/beam_grid%dv
            loop_along_track: do jj=1,ncell
                iflux = sum(states)
                ind = tracks(jj)%ind
                call get_plasma(plasma,pos=tracks(jj)%pos)
  
                call colrad(plasma,beam_ion,vnbi,tracks(jj)%time,states,dens,photons)
                call store_neutrals(ind,neut_type,dens/nlaunch)
                tracks(jj)%flux = (iflux - sum(states))*beam_grid%dv/nlaunch
  
                if(inputs%calc_birth.ge.1) then
                    call store_births(ind,neut_type,tracks(jj)%flux)
                endif
  
                if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                    call store_bes_photons(tracks(jj)%pos,vnbi,photons/nlaunch,neut_type)
                endif
            enddo loop_along_track
            if(inputs%calc_birth.ge.1) then
                do kk=1,nl_birth(neut_type)
                    call randind(tracks(1:ncell)%flux,randi)
                    call randu(randomu)
                    birth%vi(:,birth%ind) = vnbi
                    birth%ri(:,birth%ind) = tracks(randi(1))%pos + &
                                            vnbi*(tracks(randi(1))%time*(randomu(1)-0.5))
                    birth%ind = birth%ind+1
                enddo
            endif
        enddo energy_fractions
    enddo loop_over_markers
    !$OMP END PARALLEL DO
  
    if(nbi_outside.gt.0)then
         write(*,'(a, f6.2)') 'Percent of markers outside the grid: ', &
                              100.*nbi_outside/(3.*inputs%n_nbi)
         if(sum(neut%dens).eq.0) stop 'Beam does not intersect the grid!'
    endif

end subroutine ndmc

subroutine bremsstrahlung
    type(LocalProfiles) :: plasma
    integer :: i, ichan, nc, ic
    real(double) :: dlength, gaunt, max_length
    real(double) :: spot_size, theta, sqrt_rho
    real(double), dimension(2) :: randomu
    real(double), dimension(3) :: vi, xyz, r0
    real(double), dimension(3,3) :: basis
    real(double), dimension(:), allocatable :: lambda_arr,brems
  
    allocate(lambda_arr(inputs%nlambda))
    allocate(brems(inputs%nlambda))
  
    do i=1,inputs%nlambda
        lambda_arr(i)=(i-0.5)*inputs%dlambda+inputs%lambdamin ! [A]
    enddo
  
    dlength = 0.3 !cm
    !! $OMP PARALLEL DO schedule(guided) private(ichan,xyz,vi,basis,spot_size, &
    !! $OMP& ic, nc,randomu,sqrt_rho,theta,r0,plasma,gaunt,brems)
    loop_over_channels: do ichan=1,spec_chords%nchan
        xyz = spec_chords%los(ichan)%lens
        vi = spec_chords%los(ichan)%axis
        vi = vi/normp(vi)
        spot_size = spec_chords%los(ichan)%spot_size
        call line_basis(xyz,vi,basis)
  
        if(spot_size.le.0.d0) then 
            nc = 1
        else
            nc = 100
        endif
        
        do ic=1,nc
            call randu(randomu)
            sqrt_rho = sqrt(randomu(1))
            theta = 2*pi*randomu(2)
            r0(1) = 0.d0
            r0(2) = spot_size*sqrt_rho*cos(theta)
            r0(3) = spot_size*sqrt_rho*sin(theta)
            r0 = matmul(basis,r0) + xyz
  
            ! Find edge of plasma
            call get_plasma(plasma,pos=r0)
            max_length=0.0
            do while (.not.plasma%in_plasma)
                r0 = r0 + vi*dlength ! move dlength
                call get_plasma(plasma,pos=r0)
                max_length = max_length + dlength
                if(max_length.gt.300) cycle loop_over_channels
            enddo
  
            ! Calculate bremsstrahlung along los
            do while (plasma%in_plasma)
                if(plasma%te.gt.0.0) then
                    gaunt = 5.542-(3.108-log(plasma%te))*(0.6905-0.1323/plasma%zeff)
                    brems = 7.57d-9*gaunt*plasma%dene**2*plasma%zeff/(lambda_arr &
                            *sqrt(plasma%te*1000.0))*exp(-h_planck*c0/(lambda_arr*plasma%te*1000.0)) &
                            *inputs%dlambda*(4.d0*pi)*1.d-4
  
                    spec%brems(:,ichan)= spec%brems(:,ichan) + (brems*dlength*1.d-2)/nc
                endif
                ! Take a step
                r0 = r0 + vi*dlength
                call get_plasma(plasma,pos=r0)
            enddo
        enddo
  
        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') 100*ichan/real(spec_chords%nchan),char(13)
        endif
    enddo loop_over_channels
    !! $OMP END PARALLEL DO

    deallocate(lambda_arr,brems)

end subroutine bremsstrahlung

subroutine dcx
    integer :: i,j,k !indices of cells
    integer :: idcx !! counter
    real(double), dimension(3) :: ri    !! start position
    real(double), dimension(3) :: vihalo
    integer,dimension(3) :: ind    !! actual cell
    integer,dimension(3) :: neut_types = [1,2,3]
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(double), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs) :: prob    !!  Prob. for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle tracks
    integer :: jj       !! counter along track
    real(double):: tot_denn, photons  !! photon flux
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
    real(double) :: papprox_tot, ccnt, inv_ng
  
    halo_iter_dens(halo_type) = 0.d0
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
    tot_denn=0.d0
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind = [i,j,k]
                call get_plasma(plasma,ind=ind)
                tot_denn = sum(neut%dens(:,nbif_type,i,j,k)) + &
                           sum(neut%dens(:,nbih_type,i,j,k)) + &
                           sum(neut%dens(:,nbit_type,i,j,k))
                papprox(i,j,k)= tot_denn*(plasma%denp-plasma%denf)
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
  
    call get_nlaunch(inputs%n_dcx,papprox,papprox_tot,nlaunch)
  
    if(inputs%verbose.ge.1) then
       write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif
  
    ccnt=0.d0
    inv_ng = 100.0/real(beam_grid%ngrid)
    !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,idcx,ind,vihalo, &
    !$OMP& ri,tracks,ncell,prob,denn,states,jj,photons,plasma)
    loop_along_z: do k = 1, beam_grid%nz
        loop_along_y: do j = 1, beam_grid%ny
            loop_along_x: do i = 1, beam_grid%nx
                !! ------------- loop over the markers ---------------------- !!
                loop_over_dcx: do idcx=1,int(nlaunch(i,j,k))
                    !! ---------------- calculate ri,vhalo and track ----------!!
                    ind = [i, j, k]
                    call mc_halo(ind,vihalo,ri)
                    call track(ri,vihalo,tracks,ncell)
                    if(ncell.eq.0) cycle loop_over_dcx
  
                    !! ---------------- calculate CX probability --------------!!
                    call get_beam_cx_prob(tracks(1)%ind,ri,vihalo,neut_types,prob)
                    if(sum(prob).le.0.) cycle loop_over_dcx
  
                    !! --------- solve collisional radiative model along track-!!
                    call get_plasma(plasma,pos=tracks(1)%pos)
  
                    states = prob*(plasma%denp - plasma%denf)
  
                    loop_along_track: do jj=1,ncell
                        call get_plasma(plasma,pos=tracks(jj)%pos)
    
                        call colrad(plasma,thermal_ion,vihalo,tracks(jj)%time,states,denn,photons)
                        call store_neutrals(tracks(jj)%ind,halo_type,denn/nlaunch(i,j,k),plasma%in_plasma)
    
                        if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                          call store_bes_photons(tracks(jj)%pos,vihalo,photons/nlaunch(i,j,k),halo_type)
                        endif
    
                    enddo loop_along_track
                enddo loop_over_dcx
                ccnt=ccnt+1
                if (inputs%verbose.eq.2)then
                    WRITE(*,'(f7.2,"% completed",a,$)') ccnt*inv_ng,char(13)
                endif
            enddo loop_along_x
        enddo loop_along_y
    enddo loop_along_z
    !$OMP END PARALLEL DO

end subroutine dcx

subroutine halo
    integer :: i,j,k !indices of cells
    integer :: ihalo !! counter
    real(double), dimension(3) :: ri    !! start position
    real(double), dimension(3) :: vihalo!! velocity bulk plasma ion
    integer,dimension(3) :: ind    !! actual cell
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(double), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs) :: prob    !!  Prob. for CX
    !! Collisiional radiative model along track
    real(double), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks  !! Particle Tracks
    integer :: jj       !! counter along track
    real(double) :: tot_denn, photons  !! photon flux
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
    real(double) :: papprox_tot, ccnt, inv_ng
    !! Halo iteration
    integer :: hh !! counters
    real(double) :: dcx_dens, halo_iteration_dens
    integer :: s1type  ! halo iteration
    integer :: s2type  ! halo iteration
  
    s1type = fida_type
    s2type = brems_type
  
    dcx_dens = halo_iter_dens(halo_type)
    if(dcx_dens.eq.0) then
        write(*,'(a)') 'HALO: Density of DCX-neutrals is zero'
        stop
    endif
    inv_ng = 100.0/real(beam_grid%ngrid)
    neut%dens(:,s1type,:,:,:) = neut%dens(:,halo_type,:,:,:)
    iterations: do hh=1,200
        !! ------------- calculate papprox needed for guess of nlaunch --------!!
        papprox=0.d0
        papprox_tot=0.d0
        tot_denn=0.d0
        halo_iter_dens(s2type) = 0.d0
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    ind = [i,j,k]
                    call get_plasma(plasma,ind=ind)
                    tot_denn = sum(neut%dens(:,s1type,i,j,k))
                    papprox(i,j,k)= tot_denn*(plasma%denp-plasma%denf)
  
                    if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
                enddo
            enddo
        enddo
  
        call get_nlaunch(inputs%n_halo,papprox,papprox_tot,nlaunch)
  
        if(inputs%verbose.ge.1) then
            write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
        endif
        ccnt=0.d0
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(i,j,k,ihalo,ind,vihalo, &
        !$OMP& ri,tracks,ncell,prob,denn,states,jj,photons,plasma)
        loop_along_z: do k = 1, beam_grid%nz
            loop_along_y: do j = 1, beam_grid%ny
                loop_along_x: do i = 1, beam_grid%nx
                    !! ------------- loop over the markers ---------------------- !!
                    loop_over_halos: do ihalo=1,int(nlaunch(i,j,k))
                        !! ---------------- calculate ri,vhalo and track ----------!!
                        ind = [i, j, k]
                        call mc_halo(ind,vihalo,ri)
                        call track(ri,vihalo,tracks,ncell)
                        if(ncell.eq.0)cycle loop_over_halos
  
                        !! ---------------- calculate CX probability --------------!!
                        call get_beam_cx_prob(tracks(1)%ind,ri,vihalo,[s1type],prob)
                        if(sum(prob).le.0.)cycle loop_over_halos
  
                        !! --------- solve collisional radiative model along track-!!
                        call get_plasma(plasma,pos=tracks(1)%pos)
  
                        states = prob*plasma%denp
  
                        loop_along_track: do jj=1,ncell
                            call get_plasma(plasma,pos=tracks(jj)%pos)
    
                            call colrad(plasma,thermal_ion,vihalo,tracks(jj)%time,states,denn,photons)
                            call store_neutrals(tracks(jj)%ind,s2type, &
                                 denn/nlaunch(i,j,k),plasma%in_plasma)
    
                            if((photons.gt.0.d0).and.(inputs%calc_bes.ge.1)) then
                              call store_bes_photons(tracks(jj)%pos,vihalo,photons/nlaunch(i,j,k),halo_type)
                            endif
    
                        enddo loop_along_track
                    enddo loop_over_halos
                    ccnt=ccnt+1
                    if (inputs%verbose.eq.2)then
                        WRITE(*,'(f7.2,"% completed",a,$)') ccnt*inv_ng,char(13)
                    endif
                enddo loop_along_x
            enddo loop_along_y
        enddo loop_along_z
        !$OMP END PARALLEL DO
       
        halo_iteration_dens = halo_iter_dens(s2type)
        neut%dens(:,halo_type,:,:,:)= neut%dens(:,halo_type,:,:,:) &
                                           + neut%dens(:,s2type,:,:,:)
        neut%dens(:,s1type,:,:,:)= neut%dens(:,s2type,:,:,:)
        neut%dens(:,s2type,:,:,:)= 0.
  
        if(halo_iteration_dens/dcx_dens.gt.1)exit iterations
  
        inputs%n_halo=int(inputs%n_dcx*halo_iteration_dens/dcx_dens)
  
        if(inputs%n_halo.lt.inputs%n_dcx*0.01)exit iterations
    enddo iterations
    !! set the neutral density in s1type(fida_type) and s2type (brems) to 0!
    neut%dens(:,s1type,:,:,:) = 0.d0
    neut%dens(:,s2type,:,:,:) = 0.d0

end subroutine halo

subroutine fida_f
    integer :: i,j,k  !! indices  x,y,z of cells
    integer(kind=8) :: iion,ip
    real(double), dimension(3) :: ri      !! start position
    real(double), dimension(3) :: vi      !! velocity of fast ions
    integer, dimension(3) :: ind      !! new actual cell
    integer, dimension(4) :: neut_types=[1,2,3,4]
    logical :: los_intersect
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    real(double), dimension(nlevs) :: prob !! Prob. for CX
  
    !! Collisiional radiative model along track
    integer :: ncell
    integer :: jj      !! counter along track
    type(ParticleTrack),dimension(beam_grid%ntrack) :: tracks
  
    real(double) :: photons !! photon flux
    real(double), dimension(nlevs) :: states  !! Density of n-states
    real(double), dimension(nlevs) :: denn
  
    !! Number of particles to launch
    integer(kind=8) :: pcnt
    real(double) :: papprox_tot, inv_maxcnt, cnt
    integer, dimension(3,beam_grid%ngrid) :: pcell
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
  
    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    papprox_tot=0.d0
    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))* &
                                plasma%denf
                if(papprox(i,j,k).gt.0) then
                    pcell(:,pcnt)= ind
                    pcnt=pcnt+1
                endif
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt=pcnt-1
    inv_maxcnt=100.0/real(pcnt)
    call get_nlaunch(inputs%n_fida,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif
  
    !! Loop over all cells that have neutrals
    cnt=0.d0
    !$OMP PARALLEL DO schedule(guided) private(ip,i,j,k,iion,ind,vi,ri, &
    !$OMP tracks,ncell,jj,plasma,prob,denn,states,photons)
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        loop_over_fast_ions: do iion=1,int8(nlaunch(i, j, k))
            !! Sample fast ion distribution for velocity and position
            call mc_fastion(ind, ri, vi)
            if(sum(vi).eq.0) cycle loop_over_fast_ions
  
            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ncell,los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ncell.eq.0) cycle loop_over_fast_ions
  
            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(tracks(1)%ind, ri, vi, neut_types, prob)
            if(sum(prob).le.0.) cycle loop_over_fast_ions
  
            !! Calculate initial states of particle
            call get_plasma(plasma,pos=tracks(1)%pos)
            states=prob*plasma%denf
  
            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ncell
                call get_plasma(plasma,pos=tracks(jj)%pos)
  
                call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)
  
                call store_fida_photons(tracks(jj)%pos, vi, photons/nlaunch(i,j,k))
            enddo loop_along_track
        enddo loop_over_fast_ions
        cnt=cnt+1
  
        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells
    !$OMP END PARALLEL DO

end subroutine fida_f

subroutine fida_mc      
    integer :: iion,iphi
    type(FastIon) :: fast_ion
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(double) :: phi,theta
    real(double), dimension(3) :: ri      !! start position
    real(double), dimension(3) :: vi,vi_norm      !! velocity of fast ions
    !! Determination of the CX probability
    real(double), dimension(nlevs) :: denn    !!  neutral dens (n=1-4)
    real(double), dimension(nlevs) :: prob    !! Prob. for CX 
    !! Collisiional radiative model along track
    real(double), dimension(nlevs) :: states  ! Density of n-states
    integer :: ncell
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    logical :: los_intersect
    integer :: jj      !! counter along track
    real(double) :: photons !! photon flux 
    integer, dimension(4) :: neut_types=[1,2,3,4]
    real(double), dimension(3) :: xyz, uvw, r_gyro, uvw_vi
    real(double)  :: s, c
    real(double)  :: maxcnt, inv_maxcnt, cnt
    real(double), dimension(2) :: randomu
    integer(long) :: nlaunch
  
    maxcnt=particles%nparticle
    inv_maxcnt = 100.d0/maxcnt
  
    nlaunch = ceiling(dble(inputs%n_fida)/particles%nparticle)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers per mc particle: ",i7)') int(nlaunch)
    endif

    cnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(iion,iphi,fast_ion,vi,vi_norm,ri,phi,fields,tracks,s,c, &
    !$OMP& plasma,theta,randomu,xyz,uvw,uvw_vi,r_gyro,ncell,jj,prob,denn,los_intersect,states,photons)
    loop_over_fast_ions: do iion=1,particles%nparticle
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0) cycle loop_over_fast_ions
        if(fast_ion%cross_grid) then
            phi_loop: do iphi=1,nlaunch
                !! Pick random torodial angle
                call randu(randomu)
                phi = fast_ion%phi_enter + fast_ion%delta_phi*randomu(1)
                s = sin(phi) ; c = cos(phi)
                uvw(1) = fast_ion%r*c
                uvw(2) = fast_ion%r*s
                uvw(3) = fast_ion%z
                call uvw_to_xyz(uvw,xyz)
  
                if(particles%guiding_center) then          
                    !! Pick random gyroangle
                    theta = 2*pi*randomu(2)
                    call get_fields(fields,pos=xyz)
   
                    call pitch_to_vec(fast_ion%pitch,theta,fields,vi_norm)
                    vi = fast_ion%vabs*vi_norm
   
                    call gyro_correction(vi,fields,r_gyro)
                    ri = xyz - r_gyro
                else
                    uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
                    uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
                    uvw_vi(3) = fast_ion%vz
  
                    vi = matmul(beam_grid%inv_basis,uvw_vi)
                    ri = xyz
                endif
   
                call track(ri, vi, tracks, ncell, los_intersect)
                if(.not.los_intersect) cycle phi_loop
                if(ncell.eq.0)cycle phi_loop
  
                !! ---------------- calculate CX probability --------------!!
                call get_beam_cx_prob(tracks(1)%ind,ri,vi,neut_types,prob)
                if(sum(prob).le.0.)cycle phi_loop
  
                !! Calculate the spectra produced in each cell along the path
                states=prob*fast_ion%weight
                loop_along_track: do jj=1,ncell
                    call get_plasma(plasma,pos=tracks(jj)%pos)
  
                    call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)
  
                    call store_fida_photons(tracks(jj)%pos, vi, photons/dble(nlaunch), fast_ion%class)
                enddo loop_along_track
            enddo phi_loop
        endif
        cnt=cnt+1
        if (inputs%verbose.ge.2)then
          WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

end subroutine fida_mc

subroutine npa_f
    integer :: i,j,k  !! indices  x,y,z  of cells
    integer :: iion, det, ip
    real(double), dimension(3) :: ri      !! start position
    real(double), dimension(3) :: rf      !! end position
    real(double), dimension(3) :: vi      !! velocity of fast ions
    integer, dimension(3) :: ind      !! new actual cell
    integer, dimension(3,beam_grid%ngrid) :: pcell
  
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    integer, dimension(4) :: neut_types=[1,2,3,4]
    real(double), dimension(nlevs) :: prob    !! Prob. for CX
  
    !! Collisiional radiative model along track
    real(double), dimension(nlevs) :: states  !! Density of n-states
    real(double) :: flux !! flux
  
    integer :: inpa,pcnt
    real(double) :: papprox_tot, maxcnt, cnt, inv_maxcnt
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox, nlaunch !! approx. density
  
    !! ------------- calculate papprox needed for guess of nlaunch --------!!
    papprox=0.d0
    papprox_tot=0.d0
  
    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))* &
                                plasma%denf
  
                if(papprox(i,j,k).gt.0.and.(npa_chords%hit(i,j,k))) then
                    pcell(:,pcnt)= ind
                    pcnt = pcnt + 1
                endif
  
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt = pcnt - 1
    maxcnt=real(pcnt)
    inv_maxcnt = 100.0/maxcnt
  
    call get_nlaunch(inputs%n_npa,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif
  
    !! Loop over all cells that can contribute to NPA signal
    cnt=0.d0
    !$OMP PARALLEL DO schedule(guided) private(ip,i,j,k,inpa,iion,ind, &
    !$OMP& vi,ri,rf,det,plasma,prob,states,flux)
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        npa_loop: do inpa=1,npa%nloop
            loop_over_fast_ions: do iion=1,int(nlaunch(i, j, k))
                !! Sample fast ion distribution for velocity and position
                call mc_fastion(ind, ri, vi, at_guiding_center = .False.)
                if(sum(vi).eq.0)cycle loop_over_fast_ions
  
                !! Check if particle hits a NPA detector
                call hit_npa_detector(ri, vi ,det, rf)
                if(det.eq.0) cycle loop_over_fast_ions
                
                !! Calculate CX probability with beam and halo neutrals
                call get_beam_cx_prob(ind,ri,vi,neut_types,prob)
                if(sum(prob).le.0.)cycle loop_over_fast_ions
  
                !! Attenuate states as the particle move through plasma
                call get_plasma(plasma,pos=ri)
                states=prob*plasma%denf
                call attenuate(ri,rf,vi,states)
  
                !! Store NPA Flux
                flux = sum(states)*beam_grid%dv/(nlaunch(i,j,k)*real(npa%nloop))
                call store_npa(det,ri,rf,vi,flux)
            enddo loop_over_fast_ions
        enddo npa_loop
        cnt=cnt+1
        if (inputs%verbose.eq.2)then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells
    !$OMP END PARALLEL DO
    write(*,'("Number of NPA particles that hit a detector: ",i8)') npa%npart

end subroutine npa_f

subroutine npa_mc      
    integer :: iion,iphi
    type(FastIon) :: fast_ion
    type(LocalEMFields) :: fields
    type(LocalProfiles) :: plasma
    real(double) :: phi, theta
    real(double), dimension(3) :: ri, rf      !! positions
    real(double), dimension(3) :: vi, vi_norm !! velocity of fast ions
    integer :: det !! detector
    real(double), dimension(nlevs) :: prob    !! Prob. for CX 
    real(double), dimension(nlevs) :: states  ! Density of n-states
    real(double) :: flux
    integer, dimension(4) :: neut_types=[1,2,3,4]
    integer, dimension(3) :: ind
    real(double), dimension(3) :: xyz, uvw, r_gyro, uvw_vi
    real(double) :: s,c
    real(double) :: maxcnt, inv_maxcnt, cnt
    real(double), dimension(2) :: randomu
    integer(long) :: nlaunch
  
    maxcnt=particles%nparticle
    inv_maxcnt = 100.d0/maxcnt
  
    nlaunch = ceiling(dble(inputs%n_npa)/particles%nparticle)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers per mc particle: ",i7)') int(nlaunch)
    endif
    cnt=0.0
    !$OMP PARALLEL DO schedule(guided) private(iion,ind,iphi,fast_ion,vi,vi_norm,ri,rf,phi,fields,s,c, &
    !$OMP& plasma,theta,randomu,xyz,uvw,uvw_vi,r_gyro,prob,states,flux,det)
    loop_over_fast_ions: do iion=1,particles%nparticle
        fast_ion = particles%fast_ion(iion)
        if(fast_ion%vabs.eq.0)cycle loop_over_fast_ions
        if(fast_ion%cross_grid) then
            phi_loop: do iphi=1,nlaunch
                !! Pick random torodial angle
                call randu(randomu)
                phi = fast_ion%phi_enter + fast_ion%delta_phi*randomu(1)
                s = sin(phi) ; c = cos(phi)
                uvw(1) = fast_ion%r*c
                uvw(2) = fast_ion%r*s
                uvw(3) = fast_ion%z
                call uvw_to_xyz(uvw,xyz)
  
                if(particles%guiding_center) then          
                    !! Pick random gyroangle
                    theta = 2*pi*randomu(2)
                    call get_fields(fields,pos=xyz)
   
                    call pitch_to_vec(fast_ion%pitch,theta,fields,vi_norm)
                    vi = fast_ion%vabs*vi_norm
   
                    call gyro_correction(vi,fields,r_gyro)
                    ri = xyz - r_gyro
                else
                    uvw_vi(1) = c*fast_ion%vr - s*fast_ion%vt
                    uvw_vi(2) = s*fast_ion%vr + c*fast_ion%vt
                    uvw_vi(3) = fast_ion%vz
  
                    vi = matmul(beam_grid%inv_basis,uvw_vi)
                    ri = xyz
                endif
  
                !! Get beam grid indices 
                call get_indices(ri,ind)
  
                !! Check if particle hits a NPA detector
                call hit_npa_detector(ri, vi ,det, rf)
                if(det.eq.0) cycle phi_loop
                
                !! Calculate CX probability with beam and halo neutrals
                call get_beam_cx_prob(ind,ri,vi,neut_types,prob)
                if(sum(prob).le.0.)cycle phi_loop
  
                !! Attenuate states as the particle moves though plasma
                call get_plasma(plasma,pos=ri)
                states=prob*plasma%denf
                call attenuate(ri,rf,vi,states)
  
                !! Store NPA Flux
                flux = sum(states)*beam_grid%dv/dble(nlaunch)
                call store_npa(det,ri,rf,vi,flux)
            enddo phi_loop
        endif
        cnt=cnt+1
        if (inputs%verbose.ge.2)then
          WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_fast_ions
    !$OMP END PARALLEL DO

end subroutine npa_mc

subroutine fida_weights_mc
    integer :: i,j,k !! indices  x,y,z of cells
    integer(kind=8) :: iion,ip
    real(double), dimension(3) :: ri,r_gyro      !! start position
    real(double), dimension(3) :: vi,vi_norm    !! velocity of fast ions
    integer,dimension(3) :: ind      !! new actual cell
    integer,dimension(4) :: neut_types=[1,2,3,4]
    logical :: los_intersect
  
    !! Determination of the CX probability
    type(LocalProfiles) :: plasma
    type(LocalEMFields) :: fields
    real(double), dimension(nlevs) :: prob !! Prob. for CX
  
    !! Collisiional radiative model along track
    integer :: ncell
    integer :: jj,kk      !! counter along track
    type(ParticleTrack),dimension(beam_grid%ntrack) :: tracks
  
    real(double) :: photons !! photon flux
    real(double), dimension(nlevs) :: states  !! Density of n-states
    real(double), dimension(nlevs) :: denn
  
    integer :: nwav
    real(double) :: etov2, energy, pitch, phi
    real(double) :: dE, dP, dEdP
    real(double), dimension(:), allocatable :: ebarr, ptcharr
    integer, dimension(1) :: ienergy, ipitch
    real(double), dimension(3) :: randomu3
  
    !! Number of particles to launch
    integer(kind=8) :: pcnt
    real(double) :: papprox_tot,inv_maxcnt,cnt,fbm_denf,phase_area
    integer,dimension(3,beam_grid%ngrid) :: pcell
    real(double), dimension(beam_grid%nx,beam_grid%ny,beam_grid%nz) :: papprox,nlaunch !! approx. density
  
    nwav = inputs%nlambda_wght
  
    !! define arrays
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))
  
    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))
    dEdP = dE*dP
    phase_area = dEdP*real(inputs%np_wght)*real(inputs%ne_wght)
  
    !! allocate storage arrays
    allocate(fweight%weight(nwav,inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
    allocate(fweight%fida(nwav,spec_chords%nchan))
    allocate(fweight%mean_f(inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Number of Channels: ",i3)') spec_chords%nchan
        write(*,'(T2,"Nlambda: ",i4)') nwav
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T2,"LOS averaged: ",a)') "False"
    endif
  
    !! zero out arrays
    fweight%weight = 0.d0
    fweight%fida = 0.d0
    fweight%mean_f = 0.d0
  
    etov2 = 1.d0/(v2_to_E_per_amu*inputs%ab)
  
    !! Estimate how many particles to launch in each cell
    papprox=0.d0
    papprox_tot=0.d0
    pcnt=1
    do k=1,beam_grid%nz
        do j=1,beam_grid%ny
            do i=1,beam_grid%nx
                ind =[i,j,k]
                call get_plasma(plasma,ind=ind)
                papprox(i,j,k)=(sum(neut%dens(:,nbif_type,i,j,k)) + &
                                sum(neut%dens(:,nbih_type,i,j,k)) + &
                                sum(neut%dens(:,nbit_type,i,j,k)) + &
                                sum(neut%dens(:,halo_type,i,j,k)))
                if(papprox(i,j,k).gt.0) then
                    pcell(:,pcnt)= ind
                    pcnt=pcnt+1
                endif
                if(plasma%in_plasma) papprox_tot=papprox_tot+papprox(i,j,k)
            enddo
        enddo
    enddo
    pcnt=pcnt-1
    inv_maxcnt=100.0/real(pcnt)
    call get_nlaunch(10*inputs%n_fida,papprox,papprox_tot,nlaunch)
    if(inputs%verbose.ge.1) then
        write(*,'(T6,"# of markers: ",i9)') int(sum(nlaunch))
    endif
  
    !! Loop over all cells that have neutrals
    cnt=0.d0
    !$OMP PARALLEL DO schedule(guided) private(ip,i,j,k,iion,ind,vi,vi_norm,ri,ienergy,ipitch, &
    !$OMP tracks,ncell,jj,kk,plasma,fields,prob,denn,states,photons,energy,pitch,phi, &
    !$OMP r_gyro,los_intersect,randomu3,fbm_denf)
    loop_over_cells: do ip = 1, int(pcnt)
        i = pcell(1,ip)
        j = pcell(2,ip)
        k = pcell(3,ip)
        ind = [i, j, k]
        loop_over_fast_ions: do iion=1,int8(nlaunch(i, j, k))
            !! Sample fast ion distribution uniformally
            call randind(inputs%ne_wght, ienergy)
            call randind(inputs%np_wght, ipitch)
            call randu(randomu3)
            phi = 2*pi*randomu3(1)
            energy = ebarr(ienergy(1)) + dE*(randomu3(2)-0.5)
            pitch = ptcharr(ipitch(1)) + dP*(randomu3(3)-0.5)
  
            call randu(randomu3)
            ri = [beam_grid%xc(i),beam_grid%yc(j),beam_grid%zc(k)] + beam_grid%dr*(randomu3-0.5)
  
            fbm_denf = 0.0
            call get_ep_denf(energy,pitch,fbm_denf,pos=ri)
  
            !! Get velocity
            call get_fields(fields,pos=ri)
            call pitch_to_vec(pitch,phi,fields,vi_norm)
            vi = sqrt(energy*etov2)*vi_norm 
            if(energy.eq.0) cycle loop_over_fast_ions
  
            !! Correct for gyro motion
            call gyro_correction(vi,fields,r_gyro)
            ri = ri - r_gyro
  
            !! Find the particles path through the beam grid
            call track(ri, vi, tracks, ncell, los_intersect)
            if(.not.los_intersect) cycle loop_over_fast_ions
            if(ncell.eq.0) cycle loop_over_fast_ions
  
            !! Calculate CX probability with beam and halo neutrals
            call get_beam_cx_prob(tracks(1)%ind, ri, vi, neut_types, prob)
            if(sum(prob).le.0.) cycle loop_over_fast_ions
            states=prob*1.d20
  
            !! Calculate the spectra produced in each cell along the path
            loop_along_track: do jj=1,ncell
                call get_plasma(plasma,pos=tracks(jj)%pos)
  
                call colrad(plasma,beam_ion, vi, tracks(jj)%time, states, denn, photons)
  
                call store_fw_photons(ienergy(1), ipitch(1), &
                     tracks(jj)%pos, vi, fbm_denf, photons/nlaunch(i,j,k))
            enddo loop_along_track
        enddo loop_over_fast_ions
        cnt=cnt+1
  
        if(inputs%verbose.eq.2) then
            WRITE(*,'(f7.2,"% completed",a,$)') cnt*inv_maxcnt,char(13)
        endif
    enddo loop_over_cells
    !$OMP END PARALLEL DO
  
    fweight%weight = ((1.d-20)*phase_area/dEdP)*fweight%weight
    fweight%fida = ((1.d-20)*phase_area)*fweight%fida
    fweight%mean_f = ((1.d-20)*phase_area/dEdP)*fweight%mean_f
    call write_fida_weights()
  
end subroutine fida_weights_mc

subroutine fida_weights_los
    type(LocalProfiles) :: plasma, plasma_cell
    type(LocalEMFields) :: fields, fields_cell
    real(double) :: denf
    real(double) :: wght, wght_tot
    real(double) :: photons !! photon flux
    real(double) :: length
    type(ParticleTrack), dimension(beam_grid%ntrack) :: tracks
    integer :: nwav
    integer(long) :: i, j, k, ienergy
    integer(long) :: ipitch, igyro, icell, ichan
    real(double), dimension(:), allocatable :: ebarr,ptcharr,phiarr
    real(double), dimension(:,:), allocatable :: mean_f
    real(double), dimension(3) :: vi, vi_norm, vp
    real(double), dimension(3) :: vnbi_f, vnbi_h, vnbi_t, vhalo
    real(double), dimension(3) :: r_enter, r_exit
    real(double) :: vabs, dE, dP
   !! Determination of the CX probability
    real(double), dimension(nlevs) :: fdens,hdens,tdens,halodens
    real(double), dimension(nlevs) :: rates
    real(double), dimension(nlevs) :: states ! Density of n-states
    real(double), dimension(nlevs) :: denn  ! Density of n-states
    !! COLRAD
    real(double) :: dt, max_dens, dlength, sigma_pi
    real(double) :: eb, ptch, phi
    !! ---- Solution of differential equation  ---- !
    integer, dimension(3) :: ind  !!actual cell
    real(double), dimension(3) :: ri
    integer(long) :: ncell
    
    real(double):: etov2, dEdP
  
    nwav = inputs%nlambda_wght
  
    !! Define energy array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE = abs(ebarr(2)-ebarr(1))
  
    !! Define pitch array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP = abs(ptcharr(2)-ptcharr(1))
    dEdP = dE*dP
  
    !! define gyro - array
    allocate(phiarr(inputs%nphi_wght))
    do i=1,inputs%nphi_wght
        phiarr(i)=real(i-0.5)*2.d0*pi/real(inputs%nphi_wght)
    enddo
  
    !! allocate storage arrays
    allocate(fweight%fida(nwav,spec_chords%nchan))
    allocate(fweight%mean_f(inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
    allocate(fweight%weight(nwav,inputs%ne_wght,inputs%np_wght,spec_chords%nchan))
  
    allocate(mean_f(inputs%ne_wght,inputs%np_wght))
    !! zero out arrays
    fweight%weight = 0.d0
    fweight%fida = 0.d0
    fweight%mean_f = 0.d0
    mean_f = 0.d0
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Number of Channels: ",i3)') spec_chords%nchan
        write(*,'(T2,"Nlambda: ",i4)') nwav
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T2,"Ngyro: ", i3)') inputs%nphi_wght
        write(*,'(T2,"Maximum Energy: ",f7.2)') inputs%emax_wght
        write(*,'(T2,"LOS averaged: ",a)') "True"
        write(*,*) ''
    endif
  
    etov2 = 1.0/(v2_to_E_per_amu*inputs%ab)
  
    chan_loop: do ichan=1,spec_chords%nchan 
        fdens = 0.d0 ; hdens = 0.d0 ; tdens = 0.d0 ; halodens = 0.d0
        plasma = plasma*0.d0
        fields = fields*0.d0
        wght_tot = 0.d0
        mean_f = 0.d0
        do k=1,beam_grid%nz
            do j=1,beam_grid%ny
                do i=1,beam_grid%nx
                    if(spec_chords%los_inter(i,j,k)) then
                        ind = [i,j,k]
                        dlength = spec_chords%dlength(ichan,i,j,k)
                        fdens = fdens + neut%dens(:,nbif_type,i,j,k)*dlength
                        hdens = hdens + neut%dens(:,nbih_type,i,j,k)*dlength
                        tdens = tdens + neut%dens(:,nbit_type,i,j,k)*dlength
                        halodens = halodens + neut%dens(:,halo_type,i,j,k)*dlength
                        wght = sum(neut%dens(3,1:4,i,j,k))*dlength
                        
                        call get_plasma(plasma_cell,ind=ind)
                        call get_fields(fields_cell,ind=ind)
                        plasma = plasma + wght*plasma_cell
                        fields = fields + wght*fields_cell
                        do ipitch=1,inputs%np_wght
                            do ienergy=1,inputs%ne_wght
                                call get_ep_denf(ebarr(ienergy),ptcharr(ipitch),denf,ind=ind)
                                mean_f(ienergy,ipitch) = mean_f(ienergy,ipitch) + wght*denf
                            enddo
                        enddo
                        wght_tot = wght_tot + wght
                    endif
                enddo
            enddo
        enddo
  
        if(wght_tot.le.0) then 
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Skipping channel ",i3,": Neutral density is zero")') ichan
            endif
            cycle chan_loop
        else
            plasma = plasma/wght_tot
            plasma%in_plasma = .True.
            fields = fields/wght_tot
            fields%in_plasma= .True.
            call calc_perp_vectors(fields%b_norm,fields%a_norm,fields%c_norm)
            mean_f = mean_f/wght_tot  
            if(inputs%verbose.ge.1) then
                write(*,'(T4,"Channel: ",i3)') ichan
                write(*,'(T4,"Radius: ",f7.2)') spec_chords%radius(ichan)
                write(*,'(T4,"Mean Fast-ion Density: ",ES14.5)') sum(mean_f)*dEdP 
                write(*,*)''
            endif
        endif
  
        ri = plasma%pos
        vp = ri - spec_chords%los(ichan)%lens
        vnbi_f = ri - nbi%src
        vnbi_f = vnbi_f/normp(vnbi_f)*nbi%vinj
        vnbi_h = vnbi_f/sqrt(2.d0)
        vnbi_t = vnbi_f/sqrt(3.d0)
        sigma_pi = spec_chords%los(ichan)%sigma_pi
        dlength = 1.d0
  
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(eb,vabs,ptch,phi,vi,vi_norm, &
        !$OMP& ri,r_enter,r_exit,length,max_dens,ind,tracks,ncell,dt,icell,states,rates, &
        !$OMP& vhalo,denn,denf,photons,ienergy,ipitch,igyro)
        do ienergy=1,inputs%ne_wght
            do ipitch=1,inputs%np_wght
                do igyro=1,inputs%nphi_wght
                    eb = ebarr(ienergy)
                    vabs = sqrt(eb*etov2)
                    ptch = ptcharr(ipitch)
                    phi = phiarr(igyro)
                    call pitch_to_vec(ptch,phi,fields,vi_norm)
                    vi = vabs*vi_norm
  
                    call grid_intersect(ri,vi,length,r_enter,r_exit)
                    call track(r_enter, vi, tracks, ncell)
                    max_dens = 0.d0
                    do icell=1,ncell
                        ind = tracks(icell)%ind
                        tracks(icell)%flux = sum(neut%dens(3,1:4,ind(1),ind(2),ind(3)))
                        if(tracks(icell)%flux.gt.max_dens) max_dens=tracks(icell)%flux
                    enddo
                    dt = 0.d0
                    do icell=1,ncell
                        if(tracks(icell)%flux.gt.(0.5*max_dens)) then
                            dt = dt + tracks(icell)%time
                        endif
                    enddo
  
                    states=0.d0
                    call neut_rates(fdens,vi,vnbi_f,rates)
                    states = states + rates
                    call neut_rates(hdens,vi,vnbi_h,rates)
                    states = states + rates
                    call neut_rates(tdens,vi,vnbi_t,rates)
                    states = states + rates
                    do i=1,int(n_halo_neutrate)
                        call mc_halo(ind,vhalo,plasma_in=plasma)
                        call neut_rates(halodens,vi,vhalo,rates)
                        states = states + rates/real(n_halo_neutrate)
                    enddo
  
                    call colrad(plasma,beam_ion,vi,dt,states,denn,photons)
                    denf = mean_f(ienergy,ipitch)*dEdP
                    photons = photons/real(inputs%nphi_wght)
                    call store_fw_photons_at_chan(ichan, ienergy, ipitch, &
                         vp, vi, fields, dlength, sigma_pi, denf, photons)
                         
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
    enddo chan_loop
  
    call write_fida_weights()
  
end subroutine fida_weights_los

subroutine npa_weights
    type(LocalEMFields) :: fields
    real(double) :: pitch
    real(double) :: pcxa
    integer(long) :: det
    integer(long) :: ii, jj, kk, i, ic   !!indices
    integer,dimension(1) :: ipitch
    real(double), dimension(3) :: vi,vi_norm
    real(double) :: vabs, fbm_denf, dE, dP, ccnt
    real(double), dimension(nlevs) :: pcx   !! Rate coefficiants for CX
    real(double), dimension(nlevs) :: states, states_i  ! Density of n-states
    integer, dimension(4) :: neut_types=[1,2,3,4]
    real(double), dimension(3) :: pos,dpos,r_gyro
    integer(long) :: ichan
    real(double), dimension(:), allocatable :: ebarr, ptcharr
   
    !! define energy - array
    allocate(ebarr(inputs%ne_wght))
    do i=1,inputs%ne_wght
        ebarr(i)=real(i-0.5)*inputs%emax_wght/real(inputs%ne_wght)
    enddo
    dE=abs(ebarr(2)-ebarr(1))
  
    !! define pitch - array
    allocate(ptcharr(inputs%np_wght))
    do i=1,inputs%np_wght
        ptcharr(i)=real(i-0.5)*2./real(inputs%np_wght)-1.
    enddo
    dP=abs(ptcharr(2)-ptcharr(1))
  
    if(inputs%verbose.ge.1) then
        write(*,'(T2,"Number of Channels: ",i3)') npa_chords%nchan
        write(*,'(T2,"Nenergy: ",i3)') inputs%ne_wght
        write(*,'(T2,"Npitch: ",i3)') inputs%np_wght
        write(*,'(T2,"Maximum energy: ",f7.2)') inputs%emax_wght
        write(*,*) ''
    endif
  
    !! define storage arrays
    allocate(nweight%emissivity(beam_grid%nx, &
                                beam_grid%ny, &
                                beam_grid%nz, &
                                npa_chords%nchan))
  
    allocate(nweight%attenuation(inputs%ne_wght, &
                                 beam_grid%nx, &
                                 beam_grid%ny, &
                                 beam_grid%nz, &
                                 npa_chords%nchan))
  
    allocate(nweight%cx(inputs%ne_wght, &
                        beam_grid%nx, &
                        beam_grid%ny, &
                        beam_grid%nz, & 
                        npa_chords%nchan))
  
    allocate(nweight%weight(inputs%ne_wght, &
                            inputs%np_wght, &
                            npa_chords%nchan))
  
    allocate(nweight%flux(inputs%ne_wght, npa_chords%nchan))
  
    nweight%emissivity = 0.d0
    nweight%attenuation = 0.d0
    nweight%cx = 0.d0
    nweight%weight = 0.d0
    nweight%flux = 0.d0
  
    loop_over_channels: do ichan=1,npa_chords%nchan
        if(inputs%verbose.ge.1) then
            write(*,'(T4,"Channel: ",i3)') ichan
            write(*,'(T4,"Radius: ",f10.3)') npa_chords%radius(ichan)
        endif
  
        ccnt=0.d0
        !$OMP PARALLEL DO schedule(guided) collapse(3) private(ii,jj,kk,fields, &
        !$OMP& ic,det,pos,dpos,r_gyro,pitch,ipitch,vabs,vi,pcx,pcxa,states,states_i,vi_norm,fbm_denf)  
        loop_along_z: do kk=1,beam_grid%nz
            loop_along_y: do jj=1,beam_grid%ny
                loop_along_x: do ii=1,beam_grid%nx
                    if(npa_chords%phit(ii,jj,kk,ichan)%p.gt.0.d0) then
                        pos = [beam_grid%xc(ii), beam_grid%yc(jj), beam_grid%zc(kk)]
                        call get_fields(fields,pos=pos)
                        if(.not.fields%in_plasma) cycle loop_along_x
  
                        !!Determine velocity vector
                        dpos = npa_chords%phit(ii,jj,kk,ichan)%eff_rd
                        vi_norm = (dpos - pos)/normp(dpos - pos)
  
                        !!Check if it hits a detector just to make sure
                        call hit_npa_detector(pos,vi_norm,det)
                        if (det.ne.ichan) then
                            write(*,'(a)') 'NPA_WEIGHTS: Missed detector'
                            cycle loop_along_x
                        endif
  
                        !! Determine the angle between the B-field and the Line of Sight
                        pitch = dot_product(fields%b_norm,vi_norm)
                        ipitch=minloc(abs(ptcharr - pitch))
                        loop_over_energy: do ic = 1, inputs%ne_wght !! energy loop
                            vabs = sqrt(ebarr(ic)/(v2_to_E_per_amu*inputs%ab))
                            vi = vi_norm*vabs
                            !!Correct for gyro orbit
                            call gyro_correction(vi,fields,r_gyro)
  
                            fbm_denf=0
                            if (inputs%dist_type.eq.1) then
                                call get_ep_denf(ebarr(ic),pitch,fbm_denf,pos=(pos+r_gyro))
                            endif
                            if (fbm_denf.ne.fbm_denf) cycle loop_over_energy
  
                            !! -------------- calculate CX probability -------!!
                            call get_beam_cx_prob([ii,jj,kk],pos,vi,neut_types,pcx)
                            if(sum(pcx).le.0) cycle loop_over_energy
  
                            !!Calculate attenuation
                            states = pcx*1.0d14 !!needs to be large aribitrary number so colrad works
                            states_i=states
                            call attenuate(pos,dpos,vi,states)
                            pcxa=sum(states)/sum(states_i)
  
                            !$OMP CRITICAL(npa_wght)
                            nweight%attenuation(ic,ii,jj,kk,ichan) = pcxa
                            nweight%cx(ic,ii,jj,kk,ichan) = sum(pcx)
                            nweight%weight(ic,ipitch(1),ichan) = nweight%weight(ic,ipitch(1),ichan) + &
                                      2*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p*beam_grid%dv/dP
        
                            nweight%flux(ic,ichan) = nweight%flux(ic,ichan) + &
                                      2*beam_grid%dv*fbm_denf*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p
                                      !Factor of 2 above is to convert fbm to ions/(cm^3 dE (domega/4pi))
                            nweight%emissivity(ii,jj,kk,ichan)=nweight%emissivity(ii,jj,kk,ichan)+ &
                                      2*fbm_denf*sum(pcx)*pcxa*npa_chords%phit(ii,jj,kk,ichan)%p*dE
                            !$OMP END CRITICAL(npa_wght)
                        enddo loop_over_energy
                    endif
                    ccnt=ccnt+1
                    if (inputs%verbose.eq.2)then
                        WRITE(*,'(f7.2,"% completed",a,$)') ccnt/real(beam_grid%ngrid)*100,char(13)
                    endif
                enddo loop_along_x
            enddo loop_along_y
        enddo loop_along_z
        !$OMP END PARALLEL DO
  
       if(inputs%verbose.ge.1) then
           write(*,'(T4,A,ES14.5)'),'Flux:   ',sum(nweight%flux(:,ichan))*dE
           write(*,'(T4,A,ES14.5)'),'Weight: ',sum(nweight%weight(:,:,ichan))*dE*dP
           write(*,*) ''
       endif
    enddo loop_over_channels
  
    call write_npa_weights()
  
end subroutine npa_weights

end module simulation

!=============================================================================
!-------------------------------Main Program----------------------------------
!=============================================================================
program fidasim
    use simulation
    use hdf5_extra
#ifdef _OMP
    use omp_lib
#endif
    implicit none
    character(3)          :: arg = ''
    integer, dimension(8) :: time_arr,time_start,time_end !Time array
    integer               :: i,narg,nthreads,max_threads
    integer               :: hour,minu,sec
  
    call print_banner()
  
    narg = command_argument_count()
    if(narg.eq.0) then
        write(*,'(a)') "usage: ./fidasim namelist_file [num_threads]"
        stop
    else 
        call get_command_argument(1,namelist_file)
    endif
  
    !! Check if compression is possible
    call check_compression_availability()
  
    !! measure time
    call date_and_time (values=time_start)
  
    call read_inputs()
  
#ifdef _OMP
    max_threads = OMP_get_num_procs()
    if(narg.ge.2) then
        call get_command_argument(2,arg)
        read(arg,'(i3)') nthreads
    else
        nthreads = max_threads
    endif
    max_threads = min(nthreads,max_threads)
    if(inputs%verbose.ge.1) then
        write(*,'(a)') "---- OpenMP settings ----"
        write(*,'(T2,"Number of threads: ",i2)') max_threads
        write(*,*) ''
    endif
    call OMP_set_num_threads(max_threads)
#else
    max_threads = 1
#endif

    !! ----------------------------------------------------------
    !! ------ INITIALIZE THE RANDOM NUMBER GENERATOR  -----------
    !! ----------------------------------------------------------
    allocate(rng(max_threads))
    do i=1,max_threads
        call rng_init(rng(i),932117 + i)
    enddo
  
    !! ----------------------------------------------------------
    !! ------- READ GRIDS, PROFILES, LOS, TABLES, & FBM --------
    !! ----------------------------------------------------------
    call make_beam_grid()
    call read_equilibrium()
    call read_beam()
    call read_tables()
    call read_distribution()
    
    if((inputs%calc_spec.ge.1).or.(inputs%calc_fida_wght.ge.1)) then
        call read_chords()
    endif
  
    if((inputs%calc_npa.ge.1).or.(inputs%calc_npa_wght.ge.1)) then
        call read_npa()
    endif
  
    !! ----------------------------------------------------------
    !! --------------- ALLOCATE THE RESULT ARRAYS ---------------
    !! ----------------------------------------------------------
    !! neutral density array!
    allocate(neut%dens(nlevs,ntypes,beam_grid%nx,beam_grid%ny,beam_grid%nz))
    neut%dens = 0.d0
  
    !! birth profile
    if(inputs%calc_birth.ge.1) then
        allocate(birth%dens(3, & 
                            beam_grid%nx, &
                            beam_grid%ny, &
                            beam_grid%nz))
        allocate(birth%ri(3,int(inputs%n_birth*inputs%n_nbi)))
        allocate(birth%vi(3,int(inputs%n_birth*inputs%n_nbi)))
        birth%dens = 0.d0
        birth%ri = 0.d0
        birth%vi = 0.d0
    endif
  
    if(inputs%calc_spec.ge.1) then
        allocate(spec%brems(inputs%nlambda,spec_chords%nchan))
        allocate(spec%bes(inputs%nlambda,spec_chords%nchan,4))
        allocate(spec%fida(inputs%nlambda,spec_chords%nchan,particles%nclass)) 
        spec%brems = 0.d0
        spec%bes = 0.d0
        spec%fida = 0.d0
    endif
  
    if(inputs%calc_npa.ge.1)then
        npa%nchan = npa_chords%nchan
        allocate(npa%part(npa%nmax))
        allocate(npa%energy(fbm%nenergy))
        allocate(npa%flux(fbm%nenergy,npa%nchan))
        npa%energy = fbm%energy
        npa%flux = 0.0
    endif
  
    !! -----------------------------------------------------------------------
    !! --------------- CALCULATE/LOAD the BEAM and HALO DENSITY---------------
    !! -----------------------------------------------------------------------
    if(inputs%load_neutrals.eq.1) then
        call read_neutrals()
    else
        !! ----------- BEAM NEUTRALS ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'ndmc:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call ndmc
        if(inputs%calc_birth.eq.1)then
            call write_birth_profile()
        endif
        write(*,'(30X,a)') ''
  
        !! ---------- DCX (Direct charge exchange) ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'dcx:     ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call dcx()
        if(inputs%dump_dcx.eq.1) call write_dcx()
        write(*,'(30X,a)') ''
  
        !! ---------- HALO ---------- !!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'halo:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call halo()
        !! ---------- WRITE NEUTRALS ---------- !!
        call write_neutrals()
        write(*,'(30X,a)') ''
    endif
  
    !! -----------------------------------------------------------------------
    !!----------------------------- BREMSSTRAHLUNG ---------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_brems.ge.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'bremsstrahlung:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call bremsstrahlung()
        write(*,'(30X,a)') ''
    endif
  
    !! -----------------------------------------------------------------------
    !! --------------------- CALCULATE the FIDA RADIATION --------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_fida.ge.1)then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'fida:    ' , &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%dist_type.eq.1) then
            call fida_f()
        else
            call fida_mc()
        endif
        write(*,'(30X,a)') ''
    endif
  
    if(inputs%calc_spec.ge.1) call write_spectra()
  
    !! -----------------------------------------------------------------------
    !! ----------------------- CALCULATE the NPA FLUX ------------------------
    !! -----------------------------------------------------------------------
    if(inputs%calc_npa.ge.1)then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'npa:    ' , &
            time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%dist_type.eq.1) then
            call npa_f()
        else
            call npa_mc()
        endif
        write(*,'(30X,a)') ''
    endif
  
    if(inputs%calc_npa.ge.1) call write_npa()
  
    !! -------------------------------------------------------------------
    !! ----------- Calculation of weight functions -----------------------
    !! -------------------------------------------------------------------
    if(inputs%calc_fida_wght.ge.1) then
        colrad_threshold=0. !! to speed up simulation!
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'fida weight function:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        if(inputs%calc_fida_wght.eq.1) then
            call fida_weights_los()
        else
            call fida_weights_mc()
        endif
        write(*,'(30X,a)') ''
    endif
  
    if(inputs%calc_npa_wght.eq.1) then
        call date_and_time (values=time_arr)
        if(inputs%verbose.ge.1) then
            write(*,'(A,I2,":",I2.2,":",I2.2)') 'npa weight function:    ',  &
                  time_arr(5),time_arr(6),time_arr(7)
        endif
        call npa_weights()
        write(*,'(30X,a)') ''
    endif
  
    call date_and_time (values=time_arr)
    if(inputs%verbose.ge.1) then
        write(*,'(A,I2,":",I2.2,":",I2.2)') 'END: hour, minute, second: ',&
              time_arr(5),time_arr(6),time_arr(7)
    endif

    call date_and_time (values=time_end)
    hour = time_end(5) - time_start(5)
    minu = time_end(6) - time_start(6)
    sec  = time_end(7) - time_start(7)
    if (minu.lt.0.) then
        minu = minu +60
        hour = hour -1
    endif
    if (sec.lt.0.) then
        sec  = sec +60
        minu = minu -1
    endif
  
    if(inputs%verbose.ge.1) then
        write(*,'(A,18X,I2,":",I2.2,":",I2.2)') 'duration:',hour,minu,sec
    endif

end program fidasim
