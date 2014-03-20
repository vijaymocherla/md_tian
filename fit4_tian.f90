module fit4_tian
    !
    ! Purpose:
    !           Do Fit for EMT procedure.
    !
    use force
    use open_file
    use atom_class
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
real(8), dimension(20) :: B ! parameters for the fit. 20 is the maximum hard coded in NLLSQ
integer, dimension(20) :: IB
integer :: ip
! At entry, B's are initial guesses for parameters
! On exit, B's are the results of the fit.
real(8) :: X(1000,3,1000),Y(1000),RRR(1000)
integer, dimension(8) :: NARRAY ! Integer array of control parameters
real(8), dimension(8) :: ARRAY = 0.0 ! Array of statistical parameters. Use 0.0 to get default values.

! number of parameters held constant = NARRAY(4)
character(20) :: TITLE ! string of 20 characters used for title on printout

! Variables, non nllsq-related
integer :: itemp3(3), natoms,q
real(8) :: sumsq, ncells, Eref, se

character(len=100) teil_nml_out, slab_nml_out

ncells = 1.0d0/(2*rep+1)**2
IB = ibt
ip = ipc

itemp3 = shape(x_all)
npts = itemp3(1)
natoms = itemp3(3)

X(1:npts,:,1:natoms) = x_all
Y(1:npts) = y_all

teil_nml_out  = trim(fit_dir)//'emt'//trim(fitnum)//'_'//trim(name_p)//'.nml'
slab_nml_out  = trim(fit_dir)//'emt'//trim(fitnum)//'_'//trim(name_l)//'.nml'

!------------------------------------------------------------------------------------------------------------------
! CHECK EMT POTENTIAL SUBROUTINE AND WRITE RESULTS
!------------------------------------------------------------------------------------------------------------------

if (npts > 0) then

    write(*,'(/(a))')'CHECK EMT ENERGY CALCULATION IS WORKING'
    write(*,*) 'site X Y Z EMT DFT'
!    Write(10,*) 'site X Y Z EMT DFT'
    sumsq = 0

    slab%r = x_all(1,:,teil%n_atoms+1:natoms)
    if (teil%n_atoms > 0) then
        teil%r = x_all(1,:,1:teil%n_atoms)
        call emt_e(slab,teil)
    else
        call emt1_e(slab)
    end if
    Eref = Epot

    do q=2,npts
        slab%r = x_all(q,:,teil%n_atoms+1:natoms)
        if (teil%n_atoms > 0) then
            teil%r = x_all(q,:,1:teil%n_atoms)
            call emt_e(slab,teil)
        else
            call emt1_e(slab)
        end if

        if(q<11) write( *,'(1X, 5F15.8)') X(q,1,1), X(q,2,1), X(q,3,1), (Epot-Eref)*ncells, Y(q)
!        write(10,'(1X, 5F15.8)')  X(q,1,1), X(q,2,1), X(q,3,1), (energy-e_ref)/(2*rep+1)**2, Y(q)
!        write(7, '(1X, 4F16.8)')  X(q,1,1), X(q,2,1), X(q,3,1), (energy-e_ref)/(2*rep+1)**2

        sumsq=sumsq+((Epot-Eref)*ncells-Y(q))**2
    end do
    write(*,*)
!    write(10,*)
    write(*,*) 'rms error using starting parameters =',sqrt(sumsq/npts), ' Eref=', Eref
!    write(10,*) 'rms error using starting parameters =',sqrt(sumsq/npts), ' Eref=', E_ref
!    write(*,*)
!    write(10,*)

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
B(1:7)  = pars_p
B(8:14) = pars_l

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
!*********************** ARRAY ***************************
! variable default
ARRAY(1) = 0.1d0 ! AL .1
ARRAY(2) = 0.00001d0 ! DELTA - FOR DERIVATIVES .00001
ARRAY(3) = 0.00005d0 ! E - CONVERGENCE CRITERION .00005
ARRAY(4) = 4.0d0 ! FF 4.0
ARRAY(5) = 45.0d0 ! GAMCR - CRITICAL ANGLE 45.
ARRAY(6) = 2.0d0 ! T 2.
ARRAY(7) = 0.001d0 ! TAU .001
ARRAY(8) = 1d-31 ! ZETA - CRITERION FOR 1E-31
                    ! SINGULAR MATRIX
!array = 0.0

TITLE = 'EMT_fit'

CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)

pars_l = B(8:14)
pars_p = B(1:7)

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
se=0.0d0
call emt_e_fit(x_all(1,:,:nl_atoms+np_atoms), Eref)
call dev2eqdft(Eref)
call dev2aimddft(Eref)
call denseqdft(Eref)

do q=2,npts
        call emt_e_fit(x_all(q,:,:nl_atoms+np_atoms), Epot)
        sumsq=sumsq+((Epot-Eref)*ncells-Y(q))**2
        se = se+Sqrt(((Epot-Eref)*ncells-Y(q))**2)
end do
sumsq=sqrt(sumsq/npts)*1000
se=(se/npts)*1000
print *, 'rms =',  sumsq, 'meV'
print*,  'SE =', se, 'meV'

end subroutine fit

subroutine dev2eqdft(Eref)
    !
    ! Purpose:
    !           Calculate points at equilibrium positions for all 10 sites.
    !           For comparison of new fit with input DFT-equilibirum points.
    !

    integer :: npts, j
    real(8), dimension(:,:),   allocatable :: fix_p
    real(8), dimension(:,:,:), allocatable :: array
    real(8) :: energy, Eref

    call open_for_read(69,trim(fit_dir)//'/Eq_points_dft.dat')
    npts=1500
    allocate(fix_p(npts,3))
    do j = 1,npts
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

    call open_for_write(17,trim(fit_dir)//'dev2Eqdft'//trim(fitnum)//'.dat')
    do j=1,1500
        call emt_e_fit(array(j,:,:nl_atoms+np_atoms), energy)
        write(17,'(I4, 4f15.5)') (j+150-1)/150, array(j,1,5), array(j,2,5),&
                                 array(j,3,5), (energy-Eref)/(2*rep+1)**2
    end do
    close(17)

    deallocate(array, fix_p)

end subroutine dev2eqdft

subroutine denseqdft(Eref)
    !
    ! Purpose:
    !           Calculate points at equilibrium positions for all 10 sites.
    !           For comparison of new fit with input DFT-equilibirum points.
    !

    integer :: npts, j
    real(8), dimension(:,:),   allocatable :: fix_p
    real(8), dimension(:,:,:), allocatable :: array
    real(8) :: energy, Eref
    real(8), dimension(:), allocatable :: pdens

    call open_for_read(69,trim(fit_dir)//'/Eq_points_dft.dat')
    npts=1500
    allocate(fix_p(npts,3))
    do j = 1,npts
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
    do j=1,1500
        call emt_dens_fit(array(j,:,:nl_atoms+np_atoms), energy,pdens)
        write(17,'(I4, 3f15.5, f20.7)') (j+150-1)/150, array(j,1,5), array(j,2,5),&
                                 array(j,3,5), pdens(5)
    deallocate(pdens)
    end do
    close(17)

    deallocate(array, fix_p)

end subroutine denseqdft


subroutine dev2aimddft(Eref)
    !
    ! Purpose:
    !           Calculate the energy for all AIMD trajectories with new EMT parameters.
    !           For comparison between input AIMD and new fit.
    !
    real(8), intent(in) :: Eref
    real(8) :: energy
    character(len=80) :: pos_l_p, energy_l_p
    integer :: i, q, npts
    character(len=1) :: str
    character(len=80) :: nr
    real(8), dimension(:,:,:), allocatable :: d_l3, d_p3, array
    real(8), dimension(:), allocatable :: E_dft1
    character(len=3), dimension(13) :: names

    write(str,'(I1)') rep
    nr=trim(fit_dir)//'traj'//trim(fitnum)//'_'//str
    call open_for_append(1,trim(nr)//'.dat')
    names = (/'005','010','801','814','817','818','820','821','825','831','832'&
             ,'833','858'/)

    do i = 1, 13
        print*, 'Calculating traj', names(i)
        pos_l_p=trim(fit_dir)//'traj'//names(i)//'/XDATCAR.dat'
        energy_l_p=trim(fit_dir)//'traj'//names(i)//'/analyse.dat'
        call readinaimd(pos_l_p,energy_l_p, d_l3, d_p3, npts, E_dft1)
        allocate(array(npts,3,nl_atoms + np_atoms))
        do q = 1, npts
            array(q,:,:np_atoms) = d_p3(q,:,:)
            array(q,:,np_atoms+1:) = d_l3(q,:,:)
        end do
        deallocate(d_l3,d_p3)

        do q=1,npts
            call emt_e_fit(array(q,:,:nl_atoms+np_atoms), energy)
            write(1,'(I5, 2f20.10)') q, E_dft1(q)-evasp, (energy-Eref)/(2*rep+1)**2
            !write(*,'(I5, 2f20.10)') q, E_dft1(q)-evasp, (energy-Eref)/(2*rep+1)**2
        end do
        deallocate(array, E_dft1)
    end do
    close(1)

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

    c_matrix(1:2,1:2) = cell_mat(1:2,1:2)/(2*rep + 1)
    c_matrix(3,3) = cell_mat(3,3)
    ! Replication
    do q = 1,npts
        i = 1
        do l = 1, celldim(3)
        do j =-rep, rep
        do k=-rep, rep
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
            do r =-rep, rep
            do s =-rep, rep
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


end module fit4_tian
