module atom_class
    !
    ! Purpose:
    !           This module contains definitions of user types and all constants
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 18.02.2014    	Svenja M. Janke		Original
    !			Sascha Kandratsenka
    !			Dan J. Auerbach

    implicit none
    public         ! public for performance in accessing components
    save

    ! Various useful constants
    real(8), parameter          :: sqrt2    = 1.41421356237d0
    real(8), parameter          :: isqrt2   = 0.7071067811865d0
    real(8), parameter          :: sqrt3    = 1.73205080757d0
    real(8), parameter          :: pi       = 3.14159265359d0
    real(8), parameter          :: kB       = 8.61733238496d-5
    real(8), parameter          :: hbar     = 0.6582119280967d0
    real(8), parameter          :: twelfth  = 1.0d0/12.0d0
    real(8), parameter          :: tolerance = 1.0e-9
    integer, parameter          :: randseed(13)   = (/8,6,7,5,3,11,9,1,17,2,9,6,4/)
    real(8), parameter          :: e_max    = 20.0d0
    real(8), parameter          :: e_max_aimd= 100.0d0 ! cut-off for fitting aimd-input-data-set
    real(8), parameter          :: cutoff   = 1.50d0 ! Cut-off for energy calculations
    real(8)                     :: de_aimd_max = 0.0d0 ! maximal energy difference between AIMD points

    character(len=80), parameter:: fit_dir    = 'fitdata/'
    character(len=80), parameter:: fit_eq     = 'energy_eq.dat'
    character(len=80), parameter:: aimd_pos   = 'XDATCAR.dat'
    character(len=80), parameter:: aimd_e     = 'analyse.dat'
    character(len=80), parameter:: parsname(7)= (/'eta2=  ','n0=    ','E0=    ',&
                                                'lambda=','V0=    ','kappa= ','s0=    '/)

    ! geometrical factor for fcc metals
    ! beta = (16 Pi / 3)^(1/3)/Sqrt(2)
    real(8), parameter          :: beta     = 1.8093997906d0

    ! Number of nearest, next-nearest and next-next-nearest neighbours
    ! fcc only. MODIFY FOR NON-FCC-METALS
    integer, dimension(3), parameter   :: nneighs = (/12, 6, 24/)

    ! highest permitted projectile position
    real(8), parameter          :: proj_upgone = 6.1

    ! Criteria for number of bounces
    real(8), parameter          :: crit_imp = cos(pi/36)
    integer, parameter          :: step_imp = 10

    ! Conversion constants to program units
    !
    ! Program basic units
    !           Length : Ang
    !           Time   : fs
    !           Energy : eV
    ! Program derived units
    !           Mass   : eV fs^2 / A^2 = 1/103.6382 amu
    !           Angle  : radian = 180 deg
    !           bohr   : bohr = 0.5291772 Angstroem
    real(8), parameter          :: amu2mass = 103.6382d0
    real(8), parameter          :: deg2rad  = pi/180.0d0
    real(8), parameter          :: bohr2ang = 0.529177211d0
    real(8), parameter          :: p2GPa    = 160.2176565d0

    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3,n_atom)
    type atoms
        integer                                :: n_atoms  ! number of atoms
        integer                                :: nofix    ! number of non-fixed atoms
        real(8), dimension(:),   allocatable   :: dens     ! embedded electron density
        real(8), dimension(:,:), allocatable   :: r        ! positions
        real(8), dimension(:,:), allocatable   :: v        ! velocities
        real(8), dimension(:,:), allocatable   :: f        ! forces
        real(8), dimension(:,:), allocatable   :: a        ! accelerations
        real(8), dimension(:,:), allocatable   :: vp       ! predicted velocities
        real(8), dimension(:,:), allocatable   :: vc       ! corrected velocities
        real(8), dimension(:,:), allocatable   :: ao       ! old accelerations
        real(8), dimension(:,:), allocatable   :: au       ! very old accelerations
    end type atoms

    !-------------------------------------------------------------------------|
    !                  Constructor for type atoms                             |
    !-------------------------------------------------------------------------|
    !    input:  n_atoms                                                      |
    !    allocates arrays as (3,n_atom)                                       |
    !                                                                         |
    !-------------------------------------------------------------------------|

    interface atoms
        module procedure new_atoms
    end interface

contains

    function new_atoms(n_atoms)
        integer, intent(in) :: n_atoms
        type(atoms) new_atoms

        allocate(new_atoms%dens(n_atoms))
        allocate(new_atoms%r(3,n_atoms))    !   allocate
        allocate(new_atoms%v(3,n_atoms))
        allocate(new_atoms%f(3,n_atoms))
        allocate(new_atoms%a(3,n_atoms))

        allocate(new_atoms%vp(3,n_atoms))
        allocate(new_atoms%vc(3,n_atoms))
        allocate(new_atoms%ao(3,n_atoms))
        allocate(new_atoms%au(3,n_atoms))

        new_atoms%n_atoms = n_atoms         !   initialize
        new_atoms%nofix   = n_atoms
        new_atoms%dens = 0.0001d0
        new_atoms%r    = 0.0d0
        new_atoms%v    = 0.0d0
        new_atoms%f    = 0.0d0
        new_atoms%a    = 0.0d0

        new_atoms%vp = 0.0d0
        new_atoms%vc = 0.0d0
        new_atoms%ao = 0.0d0
        new_atoms%au = 0.0d0

    end function


end module
