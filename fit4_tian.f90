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

end module fit4_tian
