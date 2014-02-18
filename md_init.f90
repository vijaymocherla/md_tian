module md_init
    !
    ! Purpose:
    !    prepare (almost) everything for md simulations.
    !    Should be able to do the following things:
    !      1. Read in an initial configuration
    !      2. Construct simulations cell from the read-in configuration
    !      3. Assign initial thermal velocities to atoms
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 18.02.2014    	Svenja M. Janke		Original
    !			Sascha Kandratsenka	
    !			Dan J. Auerbach

    use atom_class
    use open_file
    use useful_things
    implicit none
    save

    integer :: start_tr     = 1     ! a trajectory to start with
    integer :: ntrajs       = 10    ! number of trajectories
    real(8) :: einc         = 5     ! incidence energy (eV)
    real(8) :: inclination  = 0     ! incidence polar angle (degree)
    real(8) :: azimuth      = 0     ! incidence azimuthal angle (degree)
                                    ! [10-1] = 60 deg
                                    ! [11-2] = 90 deg
    real(8) :: Tsurf        = 300   ! surface temperature (Kelvin)
    real(8) :: step         =-1.0d0 ! time step in fs
    integer :: nsteps       = 100   ! number of steps
    integer :: md_algo_l    = 1     ! md propagation algorithm: 1 - verlet
                                    !                           2 - beeman
                                    !                           3 - langevin
                                    !                           4 - langevin (series)
    integer :: md_algo_p    = 0     !  0 means no projectile
    integer :: pip_sign     =-1     ! -1 : read in from configuration file. Default.
                                    !  0 : assign random positions
                                    !  1 : coordinates for projectile via key word: top, fcc, hcp
                                    !  2 : take positions from file
    real(8) :: height = 6.0d0       ! Read-in height projectile
    integer :: n_confs = 1          ! Number of configurations to read in
    integer, dimension(2) :: wstep   = (/-1,1/)   ! way and interval to save data
    real(8) :: a_lat                ! lattice constant
    character(len=80) :: name_p = 'Elerium'
    character(len=80) :: key_p_pos = 'top'
    character(len=80) :: pot_p = 'emt'
    character(len=80) :: key_p = 'empty'
    real(8) :: mass_p = 1.0d0
    integer :: npars_p = 0
    character(len=80) :: name_l, pot_l, key_l
    character(len= 7) :: confname
    integer :: npars_l
    real(8) :: mass_l
    real(8) :: Epot = 0.0d0

    real(8),dimension(3,3) :: cell_mat, cell_imat ! simulation cell matrix and its inverse
    real(8), dimension(:), allocatable :: pars_l, pars_p ! potential parameters
    character(len=80) :: confname_file

    logical :: Tsurf_key = .false., md_algo_l_key = .false., md_algo_p_key = .false.


contains

subroutine simbox_init(slab, teil)
!
! Purpose:
!           Initialise the entire system:
!               1. Geometry
!               2. Interaction Potential
!               3. Velocities
!
    implicit none

    type(atoms), intent(out) :: slab, teil   ! hold r, v and f for atoms in the box

    character(len=80) :: pos_init_file
    character(len=80) :: buffer, biffer, label, mdpa_name_p,mdpa_name_l
    character(len= 1) :: coord_sys

    integer :: pos1, ios = 0, line = 0
    integer :: rep = 2  ! number of repetition layers arounf an original cell
    integer :: n_l0, n_l=1, n_p, n_p0=0, itemp
    integer :: i, j, k, l, s, r
    integer :: randk = 13
    integer :: nlnofix, nlno = -1

    integer, dimension(3) :: celldim=(/2,2,4/)  ! input cell structure

    real(8) :: cellscale, v_pdof, rtemp

    real(8), dimension(3,3) :: c_matrix, d_matrix

    real(8), dimension(:,:), allocatable :: start_l, start_p, d_l, d_p
    real(8), dimension(:,:), allocatable :: pos_l, vel_l, pos_p
!    real(8), dimension(:,:), allocatable :: d_ref

    logical :: exists

    integer :: itraj, q, nspec, npnofix

    if (iargc() == 0) stop " I need an input file"
    call getarg(1, pos_init_file)


! size of seed for random number generator
call random_seed(size=randk)

!------------------------------------------------------------------------------
!                       READ IN INPUT FILE
!                       ==================
!------------------------------------------------------------------------------

    call open_for_read(38, pos_init_file)

    ! ios < 0: end of record condition encountered or endfile condition detected
    ! ios > 0: an error is detected
    ! ios = 0  otherwise

    do while (ios == 0)
        read(38, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)
            case('start')
                read(buffer,*,iostat=ios) start_tr
            case('ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('Einc')
                read(buffer,*,iostat=ios) einc
            case('inclination')
                read(buffer,*,iostat=ios) inclination
                inclination = inclination*deg2rad
            case('azimuth')
                read(buffer,*,iostat=ios) azimuth
                azimuth = azimuth*deg2rad
            case('Tsurf')
                read(buffer,*,iostat=ios) Tsurf
                Tsurf_key = .true.
            case('step')
                read(buffer,*,iostat=ios) step
            case('nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('wstep')
                read(buffer,*,iostat=ios) wstep
                if (wstep(1)==-1) wstep(2) = nsteps + 1
            case ('projectile')
                read(buffer, *, iostat=ios) name_p, mass_p, pot_p, npars_p, &
                                            key_p, mdpa_name_p, n_p0
                md_algo_p_key = .true.
                mass_p=mass_p*amu2mass
                call lower_case(mdpa_name_p)
                select case (mdpa_name_p(1:3))
                    case ('ver')
                        md_algo_p = 1
                    case ('bee')
                        md_algo_p = 2
                    case ('lan')
                        md_algo_p = 3
                    case ('sla')
                        md_algo_p = 4
                   case default
                        print *, 'algorithm ', trim(mdpa_name_p), ' unknown'
                        stop
                end select
                allocate(pars_p(npars_p))
                call open_for_read(23,key_p)
                read(23,'(/)')
                do i = 1, npars_p
                    read(23,*) biffer, pars_p(i)
                end do
                close(23)

            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, npars_l, &
                                            key_l, mdpa_name_l, nlno
                md_algo_l_key = .true.
                mass_l=mass_l*amu2mass
                call lower_case(mdpa_name_l)
                select case (mdpa_name_l(1:3))
                    case ('ver')
                        md_algo_l = 1
                    case ('bee')
                        md_algo_l = 2
                    case ('lan')
                        md_algo_l = 3
                     case ('sla')
                        md_algo_l = 4
                   case default
                        print *, 'algorithm ', trim(mdpa_name_l), ' unknown'
                        stop
                end select
                allocate(pars_l(npars_l))
                call open_for_read(23,key_l)
                read(23,'(/)')
                do i = 1, npars_l
                    read(23,*) biffer, pars_l(i)
                end do
                close(23)
            case ('celldim')
                read(buffer, *, iostat=ios) celldim
            case ('pip')
                read(buffer, *, iostat=ios) pip_sign
                pos1 = scan(buffer, ' ')
                label = buffer(1:pos1)
                buffer = buffer(pos1+1:)
               select case(pip_sign)
                    case(0)
                        read(buffer, *, iostat=ios) height
                        if (ios .ne. 0) then
                            print *, 'Warning: You have not specified a &
                                               projectile height.'
                            print *, '         Projectile height set to 6.0 A.'
                        end if

                    case(1)
                        read(buffer, *, iostat=ios) key_p_pos, height
                        if (ios .ne. 0) then
                            read(buffer, *, iostat=ios) key_p_pos
                            if (ios == 0) then
                                print *, 'Warning: You have not specified a &
                                          projectile height.'
                                print *, '         Height set to 6.0 A.'
                            else
                                print *, 'Warning: You have not specified a &
                                                   projectile position.'
                                print *, '         Projectile position assigned to top &
                                                at 6.0 A'
                            end if
                        end if

                    case(2)
                        if (ios == 0) then
                            read(buffer, *, iostat=ios) key_p_pos
                        else
                            print *, 'Warning: You have not specified a &
                                               projectile position.'
                            print *, '         Projectile position assigned to top &
                                                at 6.0 A.'
                            pip_sign = 1
                        end if
                    case(-1)
                    case default
                        print '(A,I4,A)', 'Warning: Pip Flag value ', pip_sign, ' does not match any possible options.'
                        print *, '            Setting Pip to top at 6.0 A'
                        pip_sign = 1
                end select

            case ('rep')
                read(buffer, *, iostat=ios) rep
            case ('conf')
                read(buffer, *, iostat=ios) confname, confname_file
                if (confname == 'mxt') read(buffer, *, iostat=ios) confname, confname_file, n_confs

           case default
                print *, 'Skipping invalid label at line', line, label
            end select
        end if
    end do ! ios

    if (pip_sign == 1) then
        print *, 'Warning: You have selected option ', trim(key_p_pos),&
                 '         Your number of projectiles will be reduced to one.'
        n_p0 = 1
    end if

    close(38)

!------------------------------------------------------------------------------
!                       READ IN CONFIGURATION
!                       =====================
!------------------------------------------------------------------------------



    if (confname == 'POSCAR') then
        call open_for_read(38, confname_file)

        read(38,*) buffer
        read(38,*) cellscale
        read(38,*) c_matrix
!        read(38,*) n_l0, n_p
        read(38,'(A)') buffer
        read(buffer, *, iostat=ios) n_l0, n_p
        read(38,*) coord_sys


        a_lat = c_matrix(1,1)/celldim(1)*sqrt2

        ! Construct simulation cell matrix and its inverse
        d_matrix = 0.0d0
        d_matrix(1,1) = 1.0d0/c_matrix(1,1)
        d_matrix(2,2) = 1.0d0/c_matrix(2,2)
        d_matrix(3,3) = 1.0d0/c_matrix(3,3)
        d_matrix(1,2) = -d_matrix(2,2)*c_matrix(1,2)*d_matrix(1,1)

        cell_mat = 0.0d0
        cell_mat(1:2,1:2) = c_matrix(1:2,1:2)*(2*rep + 1)
        cell_mat(  3,  3) = c_matrix(3,3)

        cell_imat = 0.0d0
        cell_imat(1:2,1:2) = d_matrix(1:2,1:2)/(2*rep + 1)
        cell_imat(  3,  3) = d_matrix(3,3)

        ! Read in coordinates
        allocate(start_l(3,n_l0))
        read(38,*) start_l
        if (n_p > 0) then
            allocate(start_p(3,n_p))
            read(38,*) start_p
        end if

        ! Transform the read in coordinates into direct if they are cartesians:
        if (coord_sys == 'C' .or. coord_sys == 'c') then
            start_l=matmul(d_matrix,start_l)
            if (n_p > 0) start_p=matmul(d_matrix,start_p)
        end if

    !------------------------------------------------------------------------------
    !                    Replicate input cell
    !                    ====================
    !------------------------------------------------------------------------------
    ! 1. Allocate d_l-array:
    ! The array size is determined by the amount of repetitions of the cell image:
    ! the number of lattice atoms in the cell has to be multiplied by the number of
    ! permutations and then, the number of gold atoms has to be added again, since
    ! one also wants to keep the original image in the new array.
    ! temp is the number of gold atoms
    ! if the cell is translated by rep =1, then 8 new images are formed
    !                              rep =2, then 8 + 16
    !                              rep =3, then 8 + 16 + 24
    ! We have decided to replicate quadratically around the original lattice,
    ! otherwise, rep = 1 would not make much sense.
    !
        itemp=celldim(1)*celldim(2)
        n_l=itemp*celldim(3)*(2*rep+1)**2

        allocate(d_l(3,n_l))
        d_l = 0.0d0

        ! Replication
        i = 1
        s = 1
        do l = 1, celldim(3)
            do j =-rep, rep
                do k=-rep, rep
                    d_l(1,i:i+itemp-1) = start_l(1,(l-1)*itemp+1:l*itemp)+j
                    d_l(2,i:i+itemp-1) = start_l(2,(l-1)*itemp+1:l*itemp)+k
                    d_l(3,i:i+itemp-1) = start_l(3,(l-1)*itemp+1:l*itemp)
                    i = i+itemp
                end do
            end do
        end do


        allocate(pos_l(3,n_l), vel_l(3,n_l))
        pos_l = matmul(c_matrix,d_l)

        ! Exclude fixed atoms
        if (nlno == -1) then
            nlnofix = n_l/celldim(3)*(celldim(3)-1)    ! lowest layer fixed
        else
            nlnofix = n_l - nlno                       ! whatever number nlno of atoms fixed
        end if

        ! Sample velocities of lattice atoms from thermal distribution
        ! assuming the minimum energy configuration
        v_pdof = sqrt(2.0d0*kB*Tsurf/mass_l)


        vel_l = 0.0d0
        do i=1,nlnofix
            vel_l(1,i) = normal(0.0d0,v_pdof)
            vel_l(2,i) = normal(0.0d0,v_pdof)
            vel_l(3,i) = normal(0.0d0,v_pdof)
        enddo
        ! Set c.-of-m. velocity to zero
        vel_l(1,1:nlnofix) = vel_l(1,1:nlnofix) - sum(vel_l(1,1:nlnofix))/nlnofix
        vel_l(2,1:nlnofix) = vel_l(2,1:nlnofix) - sum(vel_l(2,1:nlnofix))/nlnofix
        vel_l(3,1:nlnofix) = vel_l(3,1:nlnofix) - sum(vel_l(3,1:nlnofix))/nlnofix

        if (md_algo_p > 0 .and. n_p == 0 .and. n_p0 == 0) then
            md_algo_p = 0
            print *, "Warning: Number of projectiles both in POSCAR and input file is 0."
            print *, "         Calculations will continue without projectile."
        end if

        itemp = n_p*(2*rep+1)**2

        ! Projectile initialization
        if (md_algo_p > 0) then    ! projectile existence justified

            if (n_p0 == 0) n_p0 = itemp

            allocate(d_p(3,n_p0))
            d_p = 0.0d0
            if (itemp < n_p0 ) then
                print *, "Warning: Number of projectiles larger than can be produced"
                print *, "         from repititions of POSCAR file."
                print *, "         All projectile positions set to zero."
            else
                i=1
                j=1
                l=1

                do r =-rep, rep
                do s =-rep, rep
                do l =   1, n_p

                    if (j > n_p0) exit
                    d_p(1,j) = start_p(1,l)+r
                    d_p(2,j) = start_p(2,l)+s
                    d_p(3,j) = start_p(3,l)

                    j=j+1

                end do
                end do
                end do
            end if

            allocate(pos_p(3,n_p0))
            pos_p = matmul(c_matrix,d_p)

        end if

        n_p = n_p0

        close(38)

        ! Create slab objects
        slab = atoms(n_l)
        ! Assign slab positions and velocities
        slab%r = pos_l
        slab%v = vel_l
        ! Assign the number of non-fixed atom
        slab%nofix = nlnofix

        ! Create projectile objects
        if (n_p > 0) then
            teil = atoms(n_p)
        else
            teil = atoms(1)
            teil%n_atoms = 0
            teil%nofix = 0
        end if

        if (n_p > 0) then
            teil%r = pos_p(:,1:n_p)
            deallocate(pos_p, d_p, start_p)
        end if

        deallocate(vel_l, pos_l, d_l, start_l)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                                  NOT POSCAR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    else if (confname == 'mxt') then

        open(unit=38, file=trim(confname_file)//'/mxt_conf00000001.bin', form='unformatted' )

        read(38) itemp
        if (step < 0.0d0) then
            read(38) step
        else
            read(38) rtemp
        end if
        read(38) rtemp
        read(38) rtemp
        if (.not. Tsurf_key) then
            Tsurf = rtemp
        else if ( rtemp .ne. Tsurf) then
            print '(A,f5.1,A)', 'Warning: The surface temperature from input (',Tsurf,' K) '
            print '(A,f5.1,A)', '         and configuration (',rtemp, ' K) files do not match.'
            print *,            '         Temperature set to that from input file.'
        end if
        read(38) nspec

        read(38) name_l, n_l, nlnofix, mass_l, pot_l, npars_l, key_l

        if (.not. allocated(pars_l)) allocate(pars_l(npars_l))

        read(38) pars_l, itemp
        if (.not. md_algo_l_key) md_algo_l = itemp

        read(38) a_lat        ! lattice constant makes us happy
        read(38) cell_mat     ! Cell matrix
        read(38) cell_imat    ! inverse cell matrix

        slab = atoms(n_l)
        slab%nofix = nlnofix

        read(38) slab%r, slab%v, slab%a, slab%dens

         if (.not. md_algo_p_key) then
            if (nspec > 1) then
                ! Projectile hasn't been defined in input file.
                ! Projectile defined from configuration file
                read(38) name_p, n_p, npnofix, mass_p, pot_p, npars_p, key_p
                if (.not. allocated(pars_p)) allocate(pars_p(npars_p))
                read(38) pars_p, itemp
                if (.not. md_algo_p_key) md_algo_p = itemp
                teil = atoms(n_p)
                read(38) teil%r, teil%v, teil%a, teil%dens
                teil%nofix = npnofix

            else
                ! Projectile hasn't been defined in input file.
                ! It does not exist in configuration file.
                n_p = 1
                teil = atoms(n_p)
                teil%n_atoms = 0
                teil%nofix = 0

            end if

        else
            ! Projectile has been defined in input file.
            ! Definition in input file always outrules other options.
            teil = atoms(n_p0)

        end if

        close(38)

    endif

    ! Create a directory for configuration data
    inquire(directory='conf',exist=exists)
    if (.not. exists) then
        call system('mkdir conf')
    end if
    ! Create a directory for trajectory data
    inquire(directory='traj',exist=exists)
    if (.not. exists) then
        call system('mkdir traj')
    end if

end subroutine simbox_init

subroutine traj_init(slab, teil)
!
! Purpose:
!           Initialise the entire system:
!               1. Geometry
!               2. Interaction Potential
!               3. Velocities
!
    implicit none

    type(atoms), intent(inout) :: slab, teil   ! hold r, v and f for atoms in the box

    integer :: traj_no, nspec, i
    real(8) :: dum
    integer :: ymm
    character(len=80) :: ddd
    character(len=8) :: str



!------------------------------------------------------------------------------
!                       READ IN CONFIGURATION
!                       =====================
!------------------------------------------------------------------------------


    ! Get random integer to open random trajectory for configurations
    write(str,'(I8.8)') int(ran1()*n_confs)+1

    open(unit=38, file=trim(confname_file)//'/mxt_conf'//str//'.bin', form='unformatted' )

    read(38) traj_no
    read(38) dum
    read(38) dum
    read(38) dum
    read(38) nspec

    read(38) ddd, ymm, ymm, dum, ddd, ymm, ddd
    read(38) (dum, i=1,ymm), ymm

    read(38) dum        ! lattice constant makes us happy
    read(38) (dum, i=1,9)
    read(38) (dum, i=1,9)

    read(38) slab%r, slab%v, slab%a, slab%dens

    if (.not.md_algo_p_key .and. nspec > 1) then
        read(38) ddd, ymm, ymm, dum, ddd, ymm, ddd
        read(38) (dum, i=1,ymm), ymm
        read(38) teil%r, teil%v, teil%a, teil%dens
    end if

    close(38)

end subroutine traj_init

subroutine particle_init(s)
    !
    ! Purpose:
    !           Assign random positions to the particle atom
    !           Furthermore set velocities
    !

    type(atoms) :: s
    integer :: i, j, n_p
    real(8) :: vinc
    real(8), dimension(2) :: cc1, cc2

    cc1 = (/a_lat*isqrt2,0.0d0/)
    cc2 = 0.5d0*cc1(1)*(/-1.0d0,sqrt3/)

    if (confname == 'mxt') then

        select case(pip_sign)

        case(0) ! Random positions

            do i=1,s%n_atoms
                s%r(1:2,i) = matmul((/ran1(), ran1()/),cell_mat(1:2,1:2))
                s%r(3,i)   = height
            end do

        case(1) ! Put atom in designated positions

            call lower_case(key_p_pos)
            select case (key_p_pos)
            case('top')
                s%r(1:2,1) = (/0.,0./)
            case('fcc')
                s%r(1:2,1) = (cc1 + 2.0d0*cc2)/3.0d0
            case('hcp')
                s%r(1:2,1) = (2.0d0*cc1 + cc2)/3.0d0
            case('bri')
                s%r(1:2,1) = (cc1 + cc2)*0.5d0
            end select

            s%r(1,1) = s%r(1,1) - height*tan(inclination)*cos(azimuth)
            s%r(2,1) = s%r(2,1) - height*tan(inclination)*sin(azimuth)
            s%r(3,1) = height

        case(2) ! Read in projectile positions from other file.

            call open_for_read(44,key_p_pos)
            read(44,*) n_p

            if (n_p < s%n_atoms) then

                print *, 'Error: Number of atoms in projectile configuration file'
                print *, '       is smaller than that in input file.'
                stop ' subroutine: particle_init()'

            else if (n_p > s%n_atoms) then

                print *, 'Warning: Number of atoms in projectile configuration file'
                print *, '         is larger than that in input file.'
                print '(A,I4,A)', '         Reading only first ', s%n_atoms, ' positions.'
            end if
                read(44,*) s%r

            close(44)
        case default

        end select
    end if

    if (pip_sign .ne. -1) then
    !     Assign projectile velocities
        vinc = sqrt(2.0d0*einc/mass_p)
        s%v(1,:) =  vinc*sin(inclination)*cos(azimuth)
        s%v(2,:) =  vinc*sin(inclination)*sin(azimuth)
        s%v(3,:) = -vinc*cos(inclination)
    end if

end subroutine particle_init

end module md_init

