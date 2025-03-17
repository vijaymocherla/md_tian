!!!!NONLINEAR LEAST SQUARE PROBLEM WITH BOUNDARY CONSTRAINTS

! Version 1.05

INCLUDE "mkl_rci.f90"
module dtrnlspbc_alloc
    use force

    implicit none

    real(8), dimension(:),     allocatable :: DFT_ref
    real(8), dimension(:,:,:), allocatable :: Xpoints
    integer, dimension(:),     allocatable :: held_constant
    real(8)                                :: ncells
    integer                                :: natoms

contains
end module

module dtrnlspbc
    implicit none
contains
    SUBROUTINE nllsqbc (Ydat , Xdat ,NARRAY , IB, other_natoms, other_ncells)
        use mkl_rci
        use mkl_rci_type
        use dtrnlspbc_alloc
        implicit none

        !!!!USER’S OBJECTIVE FUNCTION
        external :: obj_func

        real(8), dimension(:),     intent(in) :: Ydat
        real(8), dimension(:,:,:), intent(in) :: Xdat
        integer, dimension(8),     intent(in) :: NARRAY
        integer, dimension(14),    intent(in) :: IB
        integer,                   intent(in) :: other_natoms
        real(8),                   intent(in) :: other_ncells

        real(8), dimension(:), allocatable :: fvec, X
        integer :: M, N, i, j


        real(8), dimension(6) :: eps                     ! precisions for stop criteria
        real(8)               :: jac_eps                 ! jacobi calculation precision
        integer               :: rci_request, SUCCESSFUL ! rci parameters
        integer               :: iter=0, st_cr           ! iteration and stop-criterion counter
        integer               :: iter1, iter2            ! max number of total and trial-step iterations
        real(8)               :: rs = 100                ! inital step bound (keep fixed at 100)
        real(8)               :: r1, r2                  ! initial and final residuals
        type(handle_tr)       :: handle                  ! tr-solver handle
        integer               :: counter=0               ! work-around variable for a rare bug
        real(8), dimension(:), allocatable   :: LW, UP   ! lower and upper bounds
        real(8), dimension(:,:), allocatable :: fjac     ! jacobi matrix


        N      = 14 - NARRAY(4) ! N: number of fit variables to be fitted
        M      = NARRAY(1) - 1  ! M: number of energy data points - one reference energy
        natoms = other_natoms   ! Make them available to
        ncells = other_ncells   ! ++ fitting function

        if (N.le.0 .or. NARRAY(8).le.0) then             ! Check if all variables are held constant or
            print *, 'No fitting done.'                  !++ max number of iteration equals zero.
            return                                       ! If so, return to parent.
        end if

        allocate(fvec(M+1), Xpoints(M+1, 3, natoms), DFT_ref(M+1), X(N))
        allocate(LW(N), UP(N), fjac(M,N), held_constant(NARRAY(4)))

        ! put configurations and reference energies in objective function's scope
        Xpoints = Xdat(1:M+1, 1:3, 1:natoms)
        DFT_ref = Ydat(1:M+1)

        ! condense array holding constant parameters
        do i=1, 14
            if (IB(i) .ne. 0) then
                held_constant(i) = IB(i)
            end if
        end do
        !print *, 'held constant', held_constant


        eps (1:6) = 1.D-5     ! set stop criteria precision
        iter1     = NARRAY(8) ! set max number of iterations
        iter2     = 100       ! set max number of trial-step iterations
        jac_eps   = 1.D-8     ! set jacobi precision


        j = 1
        do i=1, 14
            if (all(held_constant .ne. i)) then
                ! copy fitting parameters into one array
                if ( i<8 ) then
                    !print *, 'copying pars_p', i, 'to X', j
                    X(j) = pars_p(i)
                elseif ( 8<=i<15 ) then
                    !print *, 'copying pars_l', i-7, 'to X', j
                    X(j) = pars_l(i-7)
                else
                    print *, 'Parameter problem in dtrnlspbc.f90'
                    stop
                end if

                ! set bounds
                if (i .eq. 3 .or. i .eq. 10) then
                    LW(j) = X(j) - abs(0.4*X(j))
	            UP(j) = 0.0
	        else
	            LW(j) = 0.0
                    UP(j) = X(j) + abs(0.4*X(j))
                end if

                ! restrict V0_H
                if (i .eq. 5) then
                    LW(j) = 0.5 * 0.427d0
                    UP(j) = 1.5 * 0.427d0
                end if

                j = j + 1
            end if
        end do

        if (j-1 .ne. N) then
            print *, 'Check number of fitconsts in *.inp file. I found', 14-j+1, &
                'fixed parameters but you specified', NARRAY(4)
            stop
        end if

        !        print *, 'j is', j
        !        print *, 'X', X
        !        print *, 'pars_l', pars_l
        !        print *, 'pars_p', pars_p
        !        print *, 'N', N
        !        print *, 'M', M
        !        print *, 'natoms', natoms
        !        print *, 'ncells', ncells
	!	write(*,*) 'X', X
        !        write(*,*) 'LW', LW
        !        write(*,*) 'UP', UP

        !!!!SET INITIAL VALUES
        fvec = 0.0
        fjac = 0.0
        !!!!INITIALIZE SOLVER (ALLOCATE MEMORY, SET INITIAL VALUES)
        !!!!handle IN/OUT: TR SOLVER handle
        !!!!N IN: NUMBER OF FUNCTION VARIABLES
        !!!!M IN: DIMENSION OF FUNCTION VALUE
        !!!!X IN: SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
        !!!!LW IN: LOWER BOUND
        !!!!UP IN: UPPER BOUND
        !!!!EPS IN: PRECISIONS FOR STOP-CRnum_iterIA
        !!!!num_iter1 IN: MAXIMUM NUMBER OF num_iterATIONS
        !!!!num_iter2 IN: MAXIMUM NUMBER OF num_iterATIONS OF CALCULATION OF TRIAL-STEP
        !!!!RS IN: INITIAL STEP BOUND
        !        print *, 'Initializing Intel Fitting'
        IF (DTRNLSPBC_INIT (handle, N, M, X, LW, UP, eps, iter1, iter2, rs) /= TR_SUCCESS) THEN
            !!!!IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN print ERROR MESSAGE
            print *, '| ERROR IN DTRNLSPBC_INIT'
            !!!!AND STOP
            call mkl_free_buffers
            STOP
        ENDIF
        !!!!SET INITIAL RCI CYCLE VARIABLES
        rci_request = 0
        SUCCESSFUL = 0
        !!!!RCI CYCLE
        !        print *, 'Starting Intel Fitting'
        DO WHILE (SUCCESSFUL == 0)
            !!!!CALL TR SOLVER
            !!!!handle IN/OUT: TR SOLVER handle
            !!!!fvec IN: VECTOR
            !!!!fjac IN: JACOBI MATRIX
            !!!!rci_request IN/OUT: RETURN NUMBER THAT DENOTES NEXT STEP FOR PERFORMING
            !            print *, 'call to DTRNLSPBC_SOLVE'
            !            print *, 'X', X
            !            print *, 'pars_l', pars_l
            !            print *, 'pars_p', pars_p

            IF (DTRNLSPBC_SOLVE (handle, fvec, fjac, rci_request) /= TR_SUCCESS) THEN
                !!!!IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN print ERROR MESSAGE
                print *, '| ERROR IN DTRNLSPBC_SOLVE'
                !!!!AND STOP
                call mkl_free_buffers
                STOP
            ENDIF

            ! introducing counter variable for infinite loop detection that rarely occurs
            if (rci_request.eq.0) then
                counter = counter + 1
            endif

            ! routine sometimes gets stuck if maximum number of iteration is reached
            ! iter1 is parameter passed in via NARRAY(8)
            ! iter  is counter for step displayed in stdout
            ! feel free to comment out and solve the bug yourself :)
            !if (iter.ge.iter1.and.rci_request_b4_solve.eq.0.and.rci_request.eq.1) then
            if (counter.ge.(iter1+20)) then
	    !  the algorithm has exceeded the maximal number of iterations
	        print *, "[ERROR]: dtrnlspbc.f90 failed. No fitting done. Change maxit parameter."
	        stop
	    endif

            !            print *, 'X', X
            !            print *, 'pars_l', pars_l
            !            print *, 'pars_p', pars_p
            !            print *, 'X after solve', X
            !        write(*,*) 'fjac', fjac
            !!!!rci_request IN/OUT: RETURN NUMBER THAT DENOTES NEXT STEP FOR PERFORMING
            !!!!ACCORDING TO rci_request VALUE WE DO NEXT STEPt
            SELECT CASE (rci_request)

                CASE (-1, -2, -3, -4, -5, -6)
                    !!!!STOP RCI CYCLE
                    SUCCESSFUL = 1
                CASE (1)
                    !!!!RECALCULATE FUNCTION VALUE
                    !!!!M IN: DIMENSION OF FUNCTION VALUE
                    !!!!N IN: NUMBER OF FUNCTION VARIABLES
                    !!!!X IN: SOLUTION VECTOR
                    !!!!RES OUT: FUNCTION VALUE F(X)
                    !                    print *, 'rci_request: call obj function'
                    !                    print *, 'BEFORE:'
                    !                    print *, 'X', X
                    !                    print *, 'pars_l', pars_l
                    !                    print *, 'pars_p', pars_p
                    !                    print *, 'Call to objective function'
                    !                    print *, 'AFTER:'
                    CALL  obj_func(M, N, X, fvec)
                !                    print *, 'fvec', fvec, 'has length', shape(fvec)
                !
                !                    write(*,*) 'X', X
                CASE (2)
                    !!!!COMPUTE JACOBI MATRIX
                    !!!!EXTENDET_POWELL IN: EXTERNAL OBJECTIVE FUNCTION
                    !!!!N IN: NUMBER OF FUNCTION VARIABLES
                    !!!!M IN: DIMENSION OF FUNCTION VALUE
                    !!!!fjac OUT: JACOBI MATRIX
                    !!!!X IN: SOLUTION VECTOR
                    !!!!jac_eps IN: JACOBI CALCULATION PRECISION
                    !                    print *, 'rci_request: call DJACOBI'
                    !                    print *, 'call to DJACOBI'
                    !                    print *, 'X before jacobi', X
                    !                    print *, 'shape jacobi', shape(fjac)
                    IF (DJACOBI (obj_func, N, M, fjac, X, jac_eps) /= TR_SUCCESS) THEN
                        !!!!IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN print ERROR MESSAGE
                        print *, '| ERROR IN DTRNLSPBC_SOLVE'
                        !!!!AND STOP
                        call mkl_free_buffers
                        STOP
                    ENDIF
                    !                    print *, 'after djacobi'
                    !                    print *, 'jacobi mat', fjac(1:M,1)
                    !                    print *, 'fvec', fvec
                    !                    print *, 'X after jacobi', X
                    ! first DFT reference is not printed b/c this is reference value
                    do i=1,M
                        if ((mod(i,10) .eq. 0) .or. (i .eq. M)) then
                            write(*,"(2i6, 2f12.4, 7f6.2 / 36x,7f6.2)") iter, i, DFT_ref(i+1), DFT_ref(i+1)-fvec(i), pars_p, pars_l
                        end if
                    end do
                    iter = iter+1
            ENDSELECT
        ENDDO

        !!!!GET SOLUTION STATUSES
        !!!!handle IN: TR SOLVER handle
        !!!!iter OUT: NUMBER OF ITERATIONS
        !!!!st_cr OUT: NUMBER OF STOP CRITERION
        !!!!r1 OUT: INITIAL RESIDUALS
        !!!!r2 OUT: FINAL RESIDUALS
        !        print *, 'Getting results'
        IF (DTRNLSPBC_GET (handle, iter, st_cr, r1, r2) /= TR_SUCCESS) THEN
            !!!!IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN print ERROR MESSAGE
            print *, '| ERROR IN DTRNLSPBC_GET'
            !!!!AND STOP
            STOP
            call mkl_free_buffers
        ENDIF
        !!!!FREE handle MEMORY
        !        print *, 'Freeing memory'
        IF (DTRNLSPBC_DELETE (handle) /= TR_SUCCESS) THEN
            !!!!IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN print ERROR MESSAGE
            print *, '| ERROR IN DTRNLSPBC_DELETE'
            !!!!AND STOP
            STOP
            call mkl_free_buffers
        ENDIF
        !!!!IF FINAL RESIDUAL IS LESS THAN REQUIRED PRECISION THEN print PASS
        !        IF (r2 < 1.D-1) THEN
        !            print *, '|DTRNLSPBC ............PASS'!, r1,r2
        !        !!!!ELSE print FAILED
        !        ELSE
        !            print *, '|DTRNLSPBC ............FAILED'!, r1,r2
        !        ENDIF
        !print *, 'st_cr', st_cr
        print *, 'initial residual', r1
        print *, 'final residual', r2

        deallocate (Xpoints, DFT_ref, X, LW, UP, fjac, fvec, held_constant)

    END SUBROUTINE nllsqbc
end module


SUBROUTINE obj_func (M, N, X, F)
    use dtrnlspbc_alloc
    implicit none

    integer, intent(in)                :: M ! M: number of energy data points - one reference energy
    integer, intent(in)                :: N ! N: number of fit variables to be fitted
    real(8), dimension(N), intent(in)  :: X ! contains N fitting parameter
    real(8), dimension(M), intent(out) :: F ! array that contains deviations from dft energies

    integer :: q      ! runs over configuarations
    integer :: i, j   ! loop variables
    real(8) :: Eref   ! holds reference energy
    real(8) :: energy ! holds energy of all other configurations

    ! copy fitted parameters in X to pars_* arrays
    !++ to make them available in emt routines
    j = 1
    do i=1, 14
        if (all(held_constant .ne. i)) then

            if ( i<8 ) then
!                print *, 'copying X', j, 'to pars_p', i
                pars_p(i) = X(j)
                j = j + 1
            elseif ( 8<=i<15 ) then
!                print *, 'copying X', j, 'to pars_l', i
                pars_l(i-7) = X(j)
                j = j + 1
            else
                print *, 'Parameter problem in dtrnlspbc.f90'
                stop
            end if
        end if
    end do
    !print *, 'in obj_func: pars_p and pars_l', pars_p, pars_l

    ! calculate reference energy
    call emt_e_fit(Xpoints(1,:,:natoms), Eref)

    ! put E(DFT)-E(EMT) in F which is the vector to be minimized
    do q=2, M+1
        call emt_e_fit(Xpoints(q,:,:natoms), energy)
        !        print *, 'energy', energy
        F(q-1) = DFT_ref(q)-(energy-Eref)*ncells
    end do

END SUBROUTINE obj_func
