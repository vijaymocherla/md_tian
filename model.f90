subroutine model( F, YDAT, XDAT, RRR, I, JP)
    use force
    implicit none

    real(8)             :: F
    real(8)             :: YDAT(1000)   !YDAT(N)
    real(8)             :: XDAT(1000,3,1000) !XDAT(N,M)
    real(8)             :: RRR(1000)
    integer             :: I, JP,ij

    real(8)             :: B, P, RES
    integer             :: IP,IB,J
    integer             :: N, M, KK

    integer             :: iteration
    real(8)             :: energy,  E_dref
    real(8), dimension(14):: denergy
    real(8), dimension(14) ::dEref
    real(8)                 :: ncells, Eref, delta

    integer :: q


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/BLK5/IB(20),IP
!    COMMON/debug/ debug(5)

    ncells = 1.0d0/((2*rep(1)+1)*(2*rep(2)+1))
    ! Select if derivatives shall be called or not.
    pars_l = B(8:14)
    pars_p = B(1:7)

    select case(jp)

    case(1)
        call emt_e_fit(xdat(1,:,:nl_atoms+np_atoms), Eref)
        call emt_e_fit(xdat(I,:,:nl_atoms+np_atoms), energy)
        energy=(energy-Eref)*ncells
        F   = energy
        RES = YDAT(I) - F

!        if ( ((mod(i,10)==0) .or. (i==N))) then
!            if (iteration > 0) write(*,1000) iteration, i, YDAT(i), F, B(1:14)
!            !write(*,1000) iteration, i, YDAT(i), YDAT(i)-F, B(1:14)
!        end if

    case(2)
        call emt_de_fit(xdat(1,:,:nl_atoms+np_atoms), Eref, dEref)
        call emt_de_fit(xdat(I,:,:nl_atoms+np_atoms), energy, denergy)


        energy=(energy-Eref)*ncells
        denergy=(denergy-dEref)*ncells
        F   = energy
        RES = YDAT(I) - F
        P(1:14)=denergy

        ! Set derivatives of those parameters zero that aren't supposed to contribute to fit.
        do ij=1,IP
            P(IB(ij)) = 0.0d0
        end do


!        !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------
!
!        if ( ((mod(i,10)==0) .or. (i==N))) then
!            !write(*,1000) iteration, i, YDAT(i), F, B(1:14)
!            write(*,1000) iteration, i, YDAT(i), YDAT(i)-F, B(1:14)
!        end if

    case(3)
        call emt_e_fit(xdat(1,:,:nl_atoms+np_atoms), Eref)
        call emt_e_fit(xdat(I,:,:nl_atoms+np_atoms), energy)
        energy=(energy-Eref)*ncells
        F   = energy
        RES = YDAT(I) - F

        if ( ((mod(i,10)==0) .or. (i==N))) then
            !if (iteration > 0) write(*,1000) iteration, i, YDAT(i), F, B(1:14)
            write(*,1000) iteration, i, YDAT(i), YDAT(i)-F, B(1:14)
        end if
    case(4)
        call emt_e_fit(xdat(1,:,:), Eref)
        call emt_e_fit(xdat(I,:,:), energy)
        energy=(energy-Eref)*ncells
        F   = energy
        RRR(I) = F

    end select


!    F   = energy
!    RES = YDAT(I) - F


   RETURN

    1000 format(2i4, 2f12.4, 7f6.2 / 32x,7f6.2)
    1010 format(2i4,3f6.2,3E12.3,14f6.2)

end subroutine MODEL
