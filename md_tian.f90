program md_tian
    ! Purpose:
    !       Do molecular dynamics calculations with the EMT potential.
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 18.02.2014    	Svenja M. Janke		Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach
    !
use atom_class
use md_init
use force
use mdalgo
use output
use fit4_tian

implicit none

type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

integer :: i, j, itraj, q, nwrites, ndata

real(8) :: imass_l, imass_p, rtemp, proj_downgone
real(8), dimension(:,:), allocatable :: output_info
real(8), dimension(:,:), allocatable :: rmin_p              ! lowest particle position
real(8), dimension(:,:,:), allocatable :: rbounce           ! lowest particle position
real(8), dimension(:,:), allocatable :: tock                ! collection array for momentum
integer, dimension(:), allocatable   :: col_start, col_end  ! collision time
integer, dimension(:), allocatable   :: imp, q_imp          ! collision number
logical :: exit_key, imp_switch
real(8), dimension(:), allocatable :: eed, eed_prec         ! embedded electron density for projectile
real(8) :: Eref, leng,mm, no                                            ! reference energy

!timing
!real(8) :: start, fin

! Construct simulation block and initialize everything
call simbox_init(slab, teil)

! Inverse Mass
imass_l = 1.0d0/mass_l
imass_p = 1.0d0/mass_p

! Lowest permitted projectile position
proj_downgone = slab%r(3,slab%n_atoms) - 3.0d0


! Allocate useful arrays
allocate(rmin_p(3,teil%n_atoms),rbounce(5,3,teil%n_atoms))
allocate(col_start(teil%n_atoms),col_end(teil%n_atoms))
nwrites=nsteps/wstep(2)
ndata= 4 + 7*teil%n_atoms ! Epot, Ekinl, Ekinp, density, r, v
allocate(output_info(ndata,nwrites))
allocate(tock(3,teil%n_atoms))
allocate(imp(teil%n_atoms), q_imp(teil%n_atoms))
allocate(eed(teil%n_atoms),eed_prec(teil%n_atoms))

!------------------------------------------------------------------------------
!
!                    CALCULATE THE REFERENCE ENERGY
!
!------------------------------------------------------------------------------

if (confname == 'poscar') then
    if (teil%n_atoms > 0) then
        select case (pes_key)
            case (0)
                call emt_e(slab,teil)
            case (1)
!                call  lj_e(slab,teil)
        end select
    else
        select case (pes_key)
            case (0)
                call emt1_e(slab)
            case (1)
!                call  lj1_e(slab)
        end select
    end if
    Eref = Epot
    print *, 'Eref = ', Eref
!    stop
    print *, slab%n_atoms
    call open_for_write(14,'trial.dat')
    write(14,'(3f15.10)') cell_mat
    write(14,*) slab%n_atoms
    write(14,'(3f15.10)') slab%r
!    write(14,'(3f15.10)') teil%r
    close(14)
!    print *, slab%n_atoms
    stop
end if

!!Check how much energy won by forming H-Au to vacuum
!! For that, the equilibrium bond-length of the H-Au molecule also has to be calculated
!mm=20
!do i=0,1! 200
!    leng = 1.0d0+ i/100.0d0
!    teil%r(:,1) = (/0.0d0, 0.0d0, 6.0d0/)
!    slab%r(:,17) = (/0.0d0, 0.0d0, 6.0d0+ 1.43 /)
!    ! 1.524 is the literature value for the Au-H bondlength (CRC)
!    ! 1.59 is the value we've been using so far. If you want to stay with the old way
!    ! of analysing the H-Au formation energy, you should be using the following line
!    ! instead of the above one:
!    ! slab%r(:,17) = (/0.0d0, 0.0d0, 6.0d0+ 1.59d0 /)
!    call emt_e(slab,teil)
!    !print *, Epot1, Epot-Eref
!    if (mm>Epot-Eref) then
!        mm=Epot-Eref
!        no = leng
!    end if
!end do
!call open_for_append(14,'Au_H_dist_minE.dat')
!!write(14,'(A5,15f15.5)') key_l(12:15),no, mm
!print *, key_l(12:15),no, mm
!close(14)
!stop


if (confname == 'fit') then
    call fit(slab,teil)
    stop
end if

!------------------------------------------------------------------------------
!
!                      LOOP OVER TRAJECTORIES
!
!------------------------------------------------------------------------------
do itraj = start_tr, ntrajs+start_tr-1

    if (sartre(itraj)) then
        print *, "Skipping already calculated trajectory ... ", itraj
        cycle
    end if

!------------------- TRAJECTORY INITIALISATION ROUTINE ------------------------

    print *, 'Running trajectory number ................ ', itraj

    call random_seed(put=randseed)  ! Seed random number generator
    do i = 1,100*(itraj - 1)
        call random_number(rtemp)   ! rotate it according to trajectory number
    end do

    exit_key    = .false.
    imp_switch  = .false.
    overwrite   = .true.
    rmin_p      = 6.10d0
    col_start   = 0
    col_end     = 0
    nwrites     = 0
    ndata       = 0
    imp         = 0
    q_imp       = 0
    eed_prec     = 0.0d0

! Skip initialisation routine for Annealing after 1st trajectory
    if ((sasteps > 0) .and. (itraj > 1)) then

    else
! Initialising all variables for start of new trajectory
        if (confname == 'mxt' .or. confname == 'geo') &
            call traj_init(slab, teil, Eref)

        if (teil%n_atoms > 0) then
            call particle_init(teil)
            select case (pes_key)
                case (0)
                    call emt(slab,teil)
                case (1)
    !                call  lj(slab,teil)
            end select
            if (md_algo_p == 3 .or. md_algo_p == 4 &
                .or. md_algo_p == 5) call ldfa(teil,imass_p)
            teil%a  = teil%f*imass_p
            teil%ao = teil%a
            teil%au = teil%ao

            if (md_algo_p == 5) then
                do i = 1, teil%n_atoms !Do we also need to multiply with the timestep here?
                    pEfric = teil%dens(i)*step*&
                            (teil%v(1,i)**2+teil%v(2,i)**2+teil%v(3,i)**2)/imass_p!*step
                end do
            end if

        else
            select case (pes_key)
                case (0)
                    call emt1(slab)
                case (1)
     !               call  lj1(slab)
            end select
        end if

        if (md_algo_l == 3 .or. md_algo_l == 4 ) call ldfa(slab, imass_l)
        slab%a  = slab%f*imass_l
        slab%ao = slab%a
        slab%au = slab%ao

        ! save initial state
        if (wstep(1)==-1) call out_short(slab, teil, Epot, Eref, &
                                         itraj, 0, rmin_p, col_end, imp, rbounce)
        tock = teil%v
    endif

    !print *, 'ntraj', itraj
    !timing
    !call cpu_time(start)


!------------------------------------------------------------------------------
!
!                         LOOP OVER TIME STEPS
!
!------------------------------------------------------------------------------
    do q = 1, nsteps
!        print *, q

!--------------------- SIMULATED ANNEALING ROUTINE ----------------------------

        if (sasteps > 0 .and. mod(q,sasteps) == 1) then
            Tsurf = Tmax - (Tmax-Tmin) * abs(2.0d0*(q + sasteps - 1)/nsteps - 1.0d0)
            if (Tsurf < 0) Tsurf = 0
        endif

!----------------------- PROPAGATION ROUTINE ----------------------------------

        call propagator_1(slab, md_algo_l, imass_l)         ! slab kick-drift

        if (teil%n_atoms > 0) then
           call propagator_1(teil, md_algo_p, imass_p)     ! projectile kick-drift
           select case (pes_key)                           ! slab-projectile forces
               case (0)
                   call emt(slab,teil)
               case (1)
!                   call  lj(slab,teil)
           end select
           eed = teil%dens                                 ! keep eed values
           call propagator_2(teil, md_algo_p, imass_p)     ! projectile kick
        else
            select case (pes_key)                          ! slab forces
                case (0)
                    call emt1(slab)
                case (1)
!                    call  lj1(slab)
            end select
        end if

        call propagator_2(slab, md_algo_l, imass_l)         ! slab kick


!-------------------- PROPERTIES CALCULATION ROUTINE --------------------------

        if (teil%n_atoms > 0) then
            do i = 1,teil%n_atoms

                ! Lowest projectile position
                if (teil%r(3,i) < rmin_p(3,i)) rmin_p(:,i) = teil%r(:,i)

                ! Collision time
                if (teil%r(3,i) < 2.0d0) then
                    if (col_start(i) == 0) col_start(i) = q
                    col_end(i) = q
                end if
                ! Collision number. Criterium:
                ! embedded projectile electron density has maximum
                if (eed(i) > eed_prec(i))  then
                    imp_switch = .true.
                else if (imp_switch .and. q > q_imp(i) + step_imp .and. eed(i)>0.25d0) then
                    imp(i) = imp(i) + 1
                    if (imp(i) < 6) rbounce(imp(i),:,i) = teil%r(:,i)
                    imp_switch = .false.
                    q_imp(i) = q
                end if
                eed_prec = eed
                ! Calculate post Electronic friction
                ! pef = Int(_t0^tend) eta_fric*v^2 dt
                if (md_algo_p == 5) then
                    call ldfa(teil, imass_p)
                    pEfric = pEfric+teil%dens(i)*step*&
                             (teil%v(1,i)**2+teil%v(2,i)**2+teil%v(3,i)**2)/imass_p
                end if

            end do
        end if

!---------------------------- OUTPUT ROUTINE ----------------------------------

        if (mod(q,wstep(2))==0) then

            ! write out
            select case(wstep(1))

                case(-1)

                case(0) ! save trajectory info
                        ! Epot, Ekinl, Ekinp, Etotal, density, r, v
                    nwrites = nwrites+1
                    output_info(1,nwrites) = Epot
                    output_info(2,nwrites) = E_kin(slab,mass_l)
                    output_info(3,nwrites) = E_kin(teil,mass_p)
                    output_info(4,nwrites) = Epot + output_info(2,nwrites)&
                                                  + output_info(3,nwrites)
                    j = teil%n_atoms
                    output_info(5    :4+  j,nwrites) = eed
                    output_info(5+  j:4+4*j,nwrites) = reshape(teil%r,(/3*j/))
                    output_info(5+4*j:4+7*j,nwrites) = reshape(teil%v,(/3*j/))
                    ndata = ndata + 1

                case(-2)
                    call out_all(slab, teil,itraj,Eref)

                case(-4)
                    call out_posvel(slab, teil, itraj, Eref)
                case(-5)
                    call out_pdb(slab, teil, q, Eref)

                case default ! full configuration of system
                    if (q > wstep(1)) call full_conf(slab, teil,itraj,Eref)

            end select

        end if

        if (mod(q,100000)==0 .and. wstep(1) == 3) then
            save_counter = save_counter - 1
            call full_conf(slab, teil,itraj,Eref)
        end if

        do i=1,teil%n_atoms
            if (teil%r(3,i) > proj_upgone .or. teil%r(3,i) < proj_downgone) then
                exit_key = .true.
            end if
        end do
        if (exit_key) exit

    end do ! steps

    col_end = col_end - col_start
    ! final state
    if (sasteps < 1) call full_conf(slab, teil,itraj,Eref)

    if (wstep(1)==-1) call out_short (slab, teil, Epot, Eref, itraj, q, rmin_p, &
                                      col_end, imp, rbounce)
    if (wstep(1)== 0) call out_detail(output_info, ndata, itraj, Eref)

    if (wstep(1)==-3 .or. wstep(1) == 0) call out_poscar(slab,teil,Epot, Eref, itraj)
!    if (wstep(1)== 0) then
!        call open_for_write(797,'/home/sjanke/git/md_tian/config.dat')
!        write(797,*) slab%r
!        write(797,*) teil%r
!        close(797)
!    end if


    !timing
    !call cpu_time(fin)
    !print *, fin - start, " seconds"
!    call open_for_write(111,'config.dat')
!    write(111,'(3f15.5)') slab%r
end do ! trajectories

if (allocated(pars_p)) deallocate(pars_p)
deallocate(rbounce,eed_prec, eed, imp, tock)
deallocate(output_info)
deallocate(col_end, col_start)
deallocate(rmin_p, pars_l)

end program md_tian
