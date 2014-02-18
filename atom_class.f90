module atom_class
    !
    ! Purpose:
    !           This module containes the definitions of all the types and all
    !           constants that are used in the program.
    !
    ! Date          Author          History of Revison
    ! ====          ======          ==================
    ! 07.01.2014    Sascha&Svenja   Dan's data structure
    ! 09.10.2013    Sascha&Svenja   Original
    !

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
    real(8), parameter          :: twelfth  = 1./12.
    integer, parameter          :: randseed(13) = (/8,6,7,5,3,11,9,1,17,2,9,6,4/)

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
