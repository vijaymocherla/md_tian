module fit4_tian
!
! Purpose:
!           Do Fit for EMT procedure.
!
use force
use open_file
use atom_class
use useful_things
use dtrnlspbc

implicit none
save

contains

subroutine fit(slab, teil)

    !------------------------------------------------------------------------------
    !
    ! NLLSQ VARIABLES
    !
    !------------------------------------------------------------------------------
    ! todo make the arrays allocatable and automatically determine npts by reading to end of data
    !------------------------------------------------------------------------------

    type(atoms) :: slab, teil
    integer :: npts ! number of points to be fitted
    integer, dimension(14) :: IB    ! fixed parameters
    integer :: ip

    real(8), dimension(:,:,:), allocatable :: X
    real(8), dimension(:), allocatable :: Y
    integer, dimension(8) :: NARRAY ! Integer array of control parameters
    character(20) :: TITLE ! string of 20 characters used for title on printout

    ! Variables, non nllsq-related
    integer :: itemp3(3), natoms,q
    real(8) :: sumsq, ncells, Eref, se,pdens, sumsq_pure

    character(len=100) teil_nml_out, slab_nml_out ! directory in which parameter-files are.

    real(8) :: mm, no, leng, Epot1               ! calculate the H-Au bondlength
    real(8), dimension(:,:), allocatable :: blen ! array to calculate the H-Au bondlength

    ncells = 1.0d0/((2*rep(1)+1)*(2*rep(2)+1)) ! 1/number of atoms in layer
    IB = ibt(1:14)
    ip = ipc

    itemp3 = shape(x_all)   ! shape of array with all coordinates
    npts = itemp3(1)        ! number of configurations
    natoms = itemp3(3)      ! number of atoms in system

    allocate(X(npts,3,natoms), Y(npts))

    X(1:npts,:,1:natoms) = x_all    ! configurations
                                    ! the first configuration is H above the perfect surface
                                    ! to calculate the reference energy.
    Y(1:npts) = y_all               ! corresponding energies


    ! Name of output parameter files.
    teil_nml_out  = trim(fit_dir)//'emt'//trim(fitnum)//'_'//trim(name_p)//'.nml'
    slab_nml_out  = trim(fit_dir)//'emt'//trim(fitnum)//'_'//trim(name_l)//'.nml'


    !------------------------------------------------------------------------------------------------------------------
    ! CHECK EMT POTENTIAL SUBROUTINE AND WRITE RESULTS
    !------------------------------------------------------------------------------------------------------------------

    if (npts > 0 ) then
        if (confname == 'fit') then ! If fit is done to energies energies

            write(*,'(/(a))')'CHECK EMT ENERGY CALCULATION IS WORKING'
            write(*,*) 'site X Y Z EMT DFT'
            sumsq = 0

            ! Calculate the EMT-reference energy.
            slab%r = x_all(1,:,teil%n_atoms+1:natoms)
            if (teil%n_atoms > 0) then
                teil%r = x_all(1,:,1:teil%n_atoms)
                call emt_e(slab,teil)
            else
                call emt1_e(slab)
            end if
            Eref = Epot

            ! Calculate the rms between input and emt at the beginning of the fit.
            do q=2,npts
                slab%r = x_all(q,:,teil%n_atoms+1:natoms)
                if (teil%n_atoms > 0) then
                    teil%r = x_all(q,:,1:teil%n_atoms)
                    call emt_e(slab,teil)
                else
                    call emt1_e(slab)
                end if

                if(q<11) write( *,'(1X, 5F15.8)') X(q,1,1), X(q,2,1), X(q,3,1), (Epot-Eref)*ncells, Y(q)

                sumsq=sumsq+((Epot-Eref)*ncells-Y(q))**2
            end do
            write(*,*)
            write(*,*) 'rms error using starting parameters =',sqrt(sumsq/npts), ' Eref=', Eref


        ! Density Fit. Not recommended, because VASP and EMT density different
        elseif (confname == 'dens') then    ! If fit is done to densities

            write(*,'(/(a))')'CHECK EMT ENERGY CALCULATION IS WORKING'
            write(*,*) 'site X Y Z EMT DFT'
            sumsq = 0

            do q=1,npts
                call emt_dens_fit(x_all(q,:,:),Epot,pdens)

                if(q<11) write( *,'(1X, 5F15.8)') X(q,1,1), X(q,2,1), X(q,3,1), pdens, Y(q)

                sumsq=sumsq+(pdens-Y(q))**2

            end do
            write(*,*) 'rms error using starting parameters =',sqrt(sumsq/npts)


        end if
    end if

    ! Here, the fitting procedure starts. So, for debugging, you might want to comment in the 'stop' .
    !stop
    !------------------------------------------------------------------------------------------------------------------
    ! SETUP FOR NLLSQ
    !------------------------------------------------------------------------------------------------------------------
    ! Usage: CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)

    !------------------ SET UP B-----------------
    !
    ! The paramters are defined as follows:
    !
    !           teil    slab
    ! eta2      1       8
    ! n0        2       9
    ! E0        3       10
    ! lambda    4       11
    ! V0        5       12
    ! kappa     6       13
    ! s0        7       14
    !
    !B(1:7)  = pars_p    ! Parameters of the particle
    !B(8:14) = pars_l    ! Parameters of the slab

    !--------------------------------------------------------------------------
    ! SET UP NARRAY
    !--------------------------------------------------------------------------
    NARRAY(1) = npts ! number of data points
    NARRAY(2) = 3 ! number of independent variables (cartesian coordinates)
    NARRAY(3) = 14 ! number of parameters
    NARRAY(4) = IP ! number of parameters to be held constant during fit.
    NARRAY(5) = 1 ! intermediate printout option (0..3)
    NARRAY(6) = 1 ! final printout option (0..3)
    NARRAY(7) = 10 ! printout unit number
    NARRAY(8) = max_iterations ! maximum number of iterations

    TITLE = 'EMT_fit'
    CALL nllsqbc (Y, x_all , NARRAY , IB, natoms, ncells)

    ! Write  new parameter files
    call open_for_write(11,teil_nml_out)
    write(11, *) 'Projectile EMT Parameters'
    write(11, *) 'Name= '//trim(name_p)
    write(11, '(A8,X,f22.8)') (trim(parsname(q)), pars_p(q) , q=1,npars_p)
    close(11)

    call open_for_write(11,slab_nml_out)
    write(11, *) 'Slab EMT Parameters'
    write(11, *) 'Name= '//trim(name_l)
    write(11, '(A8,X,f22.8)') (trim(parsname(q)), pars_l(q) , q=1,npars_l)
    close(11)

    sumsq=0.0d0
    sumsq_pure = 0.0d0
    se=0.0d0

    !Calculate energy that would be won by abstracting Au-atom from bulk
    ! and forming H-Au.
    ! We only need one H-atom per cell
    allocate(blen(3,slab%n_atoms+1))
    blen(:,1) =(/0.0d0, 0.0d0, 6.0d0/)
    blen(:,2:1+slab%n_atoms)= x_all(1,:,teil%n_atoms+1:natoms)
    np_atoms =1
    call emt_e_fit(blen,Epot)
    mm=20
    ! Calculate minimum bond length
    ! and energy won/lost by abstracting Au-atom and forming H-Au
    do q=0,590
        leng = 0.1d0+ q/100.0d0
        blen(3,18) =  6.0d0+ leng
        call emt_e_fit(blen,Epot1)
        if (mm>Epot1-Epot) then
            mm=Epot1-Epot
            no = leng
        end if
    end do
    mm=mm
    call open_for_append(14,'Au_H_dist_min.dat')
    write(14,'(A5,15f15.5)') fitnum,no, mm
    print *, fitnum, ' Bondlength H-Au = ',no, 'A; H-Au formation energy =', mm, 'A'
    close(14)
    np_atoms= teil%n_atoms

    call emt_e_fit(x_all(1,:,:nl_atoms+np_atoms), Eref)
    call dev2eqdft(Eref)        ! EMT-Energy for Equilibrium-DFT-points
    call dev2aimddft(Eref,mm)   ! How does new fit reproduce AIMD? C44?
    call denseqdft()            ! Density at 10 Equilibrium sites

    if (confname == 'fit') then
        do q=2,npts
            call emt_e_fit(x_all(q,:,:nl_atoms+np_atoms), Epot)
            sumsq=sumsq+((Epot-Eref)*ncells-Y(q))**2
            if (q <= fracaimd(1)+1) then
                sumsq_pure = sumsq_pure+((Epot-Eref)*ncells-Y(q))**2
            endif
        end do
    elseif (confname == 'dens') then
        do q=2,npts
            call emt_dens_fit(x_all(q,:,:nl_atoms+np_atoms), Epot,pdens)
            sumsq=sumsq+(pdens-Y(q))**2
        end do

    end if
    sumsq=sqrt(sumsq/npts)*1000
    sumsq_pure=sqrt(sumsq_pure/(fracaimd(1)+1))*1000
    print *, 'combined_rms =',  sumsq, 'meV'
    print *, 'energy_grid_rms =',  sumsq_pure, 'meV'

end subroutine fit

subroutine dev2eqdft(Eref)
!
! Purpose:
!           Calculate points at equilibrium positions for all 10 sites.
!           This is not a comparison to the input data that was used, but
!           a general one to a number symmetry sites on the surface that
!           is read in from the top of the file.
!

integer :: npts, j, vgls
real(8), dimension(:,:),   allocatable :: fix_p
real(8), dimension(:,:,:), allocatable :: array
real(8) :: energy, Eref

call open_for_read(69,trim(fit_dir)//'Eq_points_dft.dat') ! Contains input geometries for all 10 sites
npts=lines_in_file(581, trim(fit_dir)//'Eq_points_dft.dat')
allocate(fix_p(npts,3))
read(69,*) vgls
do j = 1,npts-1
    read(69,*) fix_p(j,:)
end do
close(69)

allocate(array(npts,3,np_atoms+nl_atoms))
do j = 1,npts
    array(j,1,:np_atoms) = x_all(1,1,:np_atoms) + fix_p(j,1)
    array(j,2,:np_atoms) = x_all(1,2,:np_atoms) + fix_p(j,2)
    array(j,3,:np_atoms) = fix_p(j,3)
    array(j,:,np_atoms+1:) = x_all(1,:,np_atoms+1:)
end do

! Calculate energy and write into output file.
call open_for_write(17,trim(fit_dir)//'dev2Eqdft'//trim(fitnum)//'.dat')
do j=1,npts-1
    call emt_e_fit(array(j,:,:nl_atoms+np_atoms), energy)
    write(17,'(I4, 4f15.5)') (j+(npts/vgls)-1)/(npts/vgls), array(j,1,5), array(j,2,5),&
        array(j,3,5), (energy-Eref)/((2*rep(1)+1)*(2*rep(2)+1))
end do
close(17)
deallocate(array, fix_p)

end subroutine dev2eqdft

subroutine denseqdft()
!
! Purpose:
!           Calculate density for all 10 sites.
!           For comparison of new fit with input DFT-equilibirum-density.
!

integer :: j, npts,vgls
real(8), dimension(:,:),   allocatable :: fix_p
real(8), dimension(:,:,:), allocatable :: array
real(8) :: energy
real(8) :: pdens

call open_for_read(69,trim(fit_dir)//'/Eq_points_dft.dat')

npts = lines_in_file(581, trim(fit_dir)//'/Eq_points_dft.dat')

allocate(fix_p(npts-1,3))
read(69,*) vgls
do j = 1,npts-1
    read(69,*) fix_p(j,:)
end do
close(69)

allocate(array(npts,3,np_atoms+nl_atoms))
do j = 1,npts
    array(j,1,:np_atoms) = x_all(1,1,:np_atoms) + fix_p(j,1)
    array(j,2,:np_atoms) = x_all(1,2,:np_atoms) + fix_p(j,2)
    array(j,3,:np_atoms) = fix_p(j,3)
    array(j,:,np_atoms+1:) = x_all(1,:,np_atoms+1:)
end do

call open_for_write(17,trim(fit_dir)//'densEqdft'//trim(fitnum)//'.dat')
do j=1,npts-1
    call emt_dens_fit(array(j,:,:nl_atoms+np_atoms), energy,pdens)
    write(17,'(I4, 3f15.5, f20.7)') (j+(npts/vgls)-1)/(npts/vgls), array(j,1,5), array(j,2,5),&
        array(j,3,5), pdens
!deallocate(pdens)
end do
close(17)

deallocate(array, fix_p)

end subroutine denseqdft

subroutine dev2aimddft(Eref, mm)
!
! Purpose:
!           Calculate the energy for all AIMD trajectories with new EMT parameters.
!           For comparison between input AIMD and new fit.
!
real(8), intent(in) :: Eref, mm
real(8) :: energy,sumsq = 0.0d0,sum1 = 0.0d0,c44= 0.0d0,c12=0.0d0, c11=0.0d0
character(len=80) :: pos_l_p, energy_l_p
integer :: i, q, npts, col = 0
character(len=1) :: str
character(len=80) :: nr
real(8), dimension(:,:,:), allocatable :: d_l3, d_p3, array
real(8), dimension(:), allocatable :: E_dft1
!character(len=3), dimension(1) :: names     ! edit #trajs here

write(str,'(2I1)') rep(1)
nr=trim(fit_dir)//'traj'//trim(fitnum)//'_'//str
call open_for_append(1,trim(nr)//'.dat')
call open_for_append(141,'c44rms.dat')
!names = (/'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019'/)   ! edit #trajs here

do i = 1, trajn     ! edit #trajs here
    print*, 'Calculating traj', trajname(i)
    pos_l_p=trim(fit_dir)//'traj'//trajname(i)//'/XDATCAR.dat'
    energy_l_p=trim(fit_dir)//'traj'//trajname(i)//'/analyse.dat'
    call readinaimd(pos_l_p,energy_l_p, d_l3, d_p3, npts, E_dft1)
    allocate(array(npts,3,nl_atoms + np_atoms))
    do q = 1, npts
        array(q,:,:np_atoms) = d_p3(q,:,:)
        array(q,:,np_atoms+1:) = d_l3(q,:,:)
    end do
    deallocate(d_l3,d_p3)

    do q=1,npts
        call emt_e_fit(array(q,:,:nl_atoms+np_atoms), energy)
        write(1,'(I5, 2f20.10)') q, E_dft1(q)-evasp, (energy-Eref)/((2*rep(1)+1)*(2*rep(2)+1))
        !write(*,'(I5, 2f20.10)') q, E_dft1(q)-evasp, (energy-Eref)/((2*rep(1)+1)*(2*rep(2)+1))
        sumsq = sumsq + (((E_dft1(q)-evasp)-((energy-Eref)/((2*rep(1)+1)*(2*rep(2)+1))))**2)
    end do
    sum1 = sum1 + sumsq
    col = col + npts
    sumsq = 0.0d0
    deallocate(array, E_dft1)
end do
close(1)
sumsq = Sqrt(sum1/col)*1000.d0
!
print*,'aimd_rms = ',sumsq, ' meV'
c44 = 3*pars_l(5)*pars_l(6)*(beta*pars_l(1)-pars_l(6))/(8*pi*pars_l(7))*p2GPa
c11 = (3*pars_l(5)*(beta*pars_l(1)-pars_l(6))*pars_l(6)-pars_l(3)*pars_l(4)**2)&
    /(12*pi*pars_l(7))*p2GPa
c12 = (3*pars_l(5)*(-beta*pars_l(1)+pars_l(6))*pars_l(6)-2*pars_l(3)*pars_l(4)**2)&
    /(24*pi*pars_l(7))*p2GPa
print *, 'C44 = ', c44, ' GPa'
print *, 'C11 = ', c11, ' GPa'
print *, 'C12 = ', c12, ' GPa'
write(141,'(A6, 17f15.6)') fitnum, c44, mm, sumsq, pars_l(1), pars_l(2), pars_l(3),&
    pars_l(4), pars_l(5), pars_l(6), pars_l(7), pars_p(1), pars_p(2), pars_p(3),&
    pars_p(4), pars_p(5), pars_p(6), pars_p(7)

close(141)

end subroutine dev2aimddft

subroutine readinaimd(pos_l_p,energy_l_p, d_l3, d_p3, npts, E_dft1)
    !
    ! Purpose:
    !           Read in points and energies for routine dev2aimddft (or whenever else
    !           but for the input you might need them)

    character(len=80), intent(in) :: pos_l_p, energy_l_p
    character(len=80) :: buffer
    integer :: npts, itemp
    integer :: n_l0, n_p0
    integer :: i, j, ios,q, k, r,s, l
    real(8), dimension(:),     allocatable :: E_dft1
    real(8), dimension(:,:,:), allocatable :: aimd_l, aimd_p, d_l3, d_p3
    real(8) :: rtemp
    real(8) :: c_matrix(3,3) = 0.0d0

    call open_for_read(17, energy_l_p)
    call open_for_read(18, pos_l_p)

    ! The coefficients of the transformation matrix.
    read(18,'(A,/////)') buffer
    n_p0 = 0
    read(18,'(A)') buffer
    read(buffer, *) n_l0, n_p0

    ! According to the energy file, we define number of configurations
    i = 0
    do
        read(17,*,iostat=ios)
        if(ios <0) exit
        i = i + 1
    end do
    rewind(17)
    npts = i - 3

    allocate(E_dft1(npts))
    allocate(aimd_l(npts,3,n_l0))
    allocate(aimd_p(npts,3,n_p0))


    read(17,'(A)') buffer
    read(17,*) rtemp, rtemp, E_dft1(1)
    read(18,*) aimd_l(1,:,:)
    if (n_p0 > 0) read(18,*) aimd_p(1,:,:)

    j=2
    do i=2,npts
        read(17,*) rtemp, rtemp, E_dft1(j)
        read(18,*) aimd_l(j,:,:)
        if (n_p0 > 0) read(18,*) aimd_p(j,:,:)
        rtemp = E_dft1(j-1) - E_dft1(j)
        j = j + 1
    end do
    npts = j - 1
    close(17)
    close(18)

    itemp=celldim(1)*celldim(2)
    allocate(d_l3(npts,3,nl_atoms))

    c_matrix(1:2,1) = cell_mat(1:2,1)/(2*rep(1) + 1)
    c_matrix(1:2,2) = cell_mat(1:2,2)/(2*rep(2) + 1)
    c_matrix(3,3) = cell_mat(3,3)
    ! Replication
    do q = 1,npts
        i = 1
        do l = 1, celldim(3)
            do j =-rep(1), rep(1)
                do k =-rep(2), rep(2)
                    d_l3(q,1,i:i+itemp-1) = aimd_l(q,1,(l-1)*itemp+1:l*itemp)+j
                    d_l3(q,2,i:i+itemp-1) = aimd_l(q,2,(l-1)*itemp+1:l*itemp)+k
                    d_l3(q,3,i:i+itemp-1) = aimd_l(q,3,(l-1)*itemp+1:l*itemp)
                    i = i+itemp
                end do
            end do
        end do
        d_l3(q,:,:)= matmul(c_matrix,d_l3(q,:,:))
    end do

    ! Projectile initialization
    if (n_p0 > 0) then    ! projectile existence justified

        allocate(d_p3(npts,3,np_atoms))
        j=1
        do q = 1,npts
            j = 1
            do r =-rep(1), rep(1)
                do s =-rep(2), rep(2)
                    do l =   1, n_p0
                        d_p3(q,1,j) = aimd_p(q,1,l)+r
                        d_p3(q,2,j) = aimd_p(q,2,l)+s
                        d_p3(q,3,j) = aimd_p(q,3,l)
                        j=j+1
                    end do
                end do
            end do
            d_p3(q,:,:)= matmul(c_matrix,d_p3(q,:,:))
        end do

    end if
    deallocate(aimd_l, aimd_p)


end subroutine readinaimd

subroutine restricted_fit(ip, ib) ! This does not work.
    !
    ! Purpose:
    !           Restrict fitting parameters to 'save' shear modulus
    !           Parameters are modified via RANDOM NUMBERS
    !

    integer                 :: ip
    integer, dimension(20)  :: ib
    integer                 :: j
    real(8) :: thresh, randomnumber, other_randomnumber, third_randomnumber

    ! Create random numbers here.
    randomnumber = ran1()*3
    other_randomnumber = ran1()*10
    third_randomnumber = ran1()
    !print *, randomnumber, other_randomnumber, third_randomnumber


    ! if kappa is too large, lower kappa and increase eta.
    thresh = beta*pars_l(1)
    if (pars_l(6) > thresh ) then
        pars_l(6) = pars_l(6) - randomnumber
        pars_l(1) = pars_l(1) + other_randomnumber
    else
        ! if a parameter is restricted, add random number to this one. In case of kappa, lower it
        ! if all parameters run, add random number to all of them
        do j = 1,ip
            if (ib(j) == 8) then   ! eta
                pars_l(1) = pars_l(1) + randomnumber

            else if (ib(j) == 12) then  ! vo
                pars_l(5) = pars_l(5) + randomnumber

            elseif (ib(j) == 13) then  ! kappa
                pars_l(6) = pars_l(6) - randomnumber
                if (pars_l(6) > 0.0d0) pars_l(6) = pars_l(6) + other_randomnumber

            else
                pars_l(1) = pars_l(1) + randomnumber
                pars_l(5) = pars_l(5) + other_randomnumber
                pars_l(6) = pars_l(6) + third_randomnumber
            end if
        end do
    end if

    ! a) if certain parameter restricted, add random number to this one (vo, eta), lower kappa
    ! b) if all parameters are running, add random number to all of them. If difficulties here, make sure
    !    that kappa is lower than the other values
    ! c) If (kappa > beta*eta), lower kappa, increase eta, both with random numbers.

end subroutine restricted_fit


end module fit4_tian
