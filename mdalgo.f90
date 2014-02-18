module mdalgo
!
!
! Purpose :
!           Contains all the algorithms needed for the md procedure
!
! Date          Author              History of Revision
! ====          ======              ===================
! 23.10.2013    Sascha & Svenja     Original

    use atom_class
    use useful_things
    use open_file
    use force
    implicit none

contains

subroutine propagator_1(s, md_algo, imass)
    !
    ! Purpose:
    !           select algorithm and calculate 1st and 2nd step
    !

    type(atoms) :: s
    integer :: md_algo
    real(8) :: imass

    select case (md_algo)

        case (1) ! velocity Verlet Algorithm
            call verlet_1(s)

        case (2) ! Refson-Beeman Algorithm
            call beeman_1(s)

        case (3) ! Langevin
            call langevin_1(s,imass)

        case (4) ! Langevin up to 2nd order
            call langevins_1(s,imass)

   end select

end subroutine propagator_1

subroutine propagator_2(s, md_algo, imass)
    !
    ! Purpose:
    !           select algorithm and calculate 3rd and 4th step
    !

    type(atoms) :: s
    integer     :: md_algo
    real(8)     :: imass, norm


        select case (md_algo)

            case (1) ! velocity Verlet Algorithm
                call newton(s, imass)
                call verlet_2(s)

            case (2) ! Refson-Beeman Algorithm
                do
                    call newton(s, imass)
                    call beeman_2(s)
                    call norm_dist(s%vp, s%vc, 3*s%n_atoms, norm)
                    if (norm < 1.0e-007) then
                        s%v  = s%vc
                        s%au = s%ao
                        s%ao = s%a
                        exit
                    else
                        s%vp = s%vc
                    end if
                end do

            case (3) ! Langevin

                call ldfa(s)
                call newton(s, imass)
                call langevin_2(s,imass)

            case (4) ! Langevin up to 2nd order
                call ldfa(s)
                call newton(s, imass)
                call langevins_2(s,imass)

       end select



end subroutine propagator_2

subroutine verlet_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of velocity Verlet algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    ! half-step velocities
    s%v(:,1:nf) = s%v(:,1:nf) + 0.5d0*step*s%a(:,1:nf)
    ! new positions
    s%r(:,1:nf) = s%r(:,1:nf) + step*s%v(:,1:nf)

end subroutine verlet_1

subroutine verlet_2(s)
    !
    ! Purpose:
    !           3rd step of velocity Verlet algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    ! new velocities
    s%v(:,1:nf) = s%v(:,1:nf) + 0.5d0*step*s%a(:,1:nf)

end subroutine verlet_2

subroutine beeman_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of Refson-Beeman algorithm,
    !           K. Refson, Physica 131B, (1985), 256.
    !           Moldy User's Manual.
    !
    type(atoms) :: s
    real(8)     :: step_sq
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    step_sq = step * step / 6.0d0

    ! new positions
    s%r(:,1:nf) = s%r(:,1:nf) + step*s%v(:,1:nf) &
                + step_sq*(4.0*s%ao(:,1:nf) - s%au(:,1:nf))
    ! predicted velocities
    s%vp(:,1:nf) = s%v(:,1:nf) + 0.5d0*step*(3.0d0*s%ao(:,1:nf)&
                                                 - s%au(:,1:nf))

end subroutine beeman_1


subroutine beeman_2(s)
    !
    ! Purpose:
    !           4th step of Refson-Beeman algorithm,
    !           K. Refson, Physica 131B, (1985), 256.
    !           Moldy User's Manual.
    !
    type(atoms) :: s
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    s%vc(:,1:nf) = s%v(:,1:nf) + step/6.0d0*(2.0d0*s%a(:,1:nf)&
                           + 5.0d0*s%ao(:,1:nf) - s%au(:,1:nf))

end subroutine beeman_2

subroutine langevin_1(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Dellago et al. JCP 108 (1998) 1964
    !

    type(atoms)                   :: s
    real(8)                       :: imass, temp
    real(8), dimension(  s%nofix) :: c0, c1, c2, xidt, ixidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%nofix) :: randy, cofm
    integer                       :: nf, i

    nf = s%nofix

    temp = kB*Tsurf*imass
    xidt = (s%dens(1:nf)*step)
    ixidt = step/xidt

    do i = 1, nf

        ! Preventing problems due to precision issues
        if (xidt(i) > 1.0d-2) then                      ! use precise expressions

            c0(i) = exp(-xidt(i))
            c1(i) = (1.0d0 - c0(i))*ixidt(i)
            c2(i) = (1.0d0 - c1(i)/step)*ixidt(i)

            sigma_r(i) = ixidt(i)*sqrt(temp*(2.d0*xidt(i) - 3.d0 + 4.d0*c0(i) - c0(i)**2))
            sigma_v(i) = sqrt(temp*(1.d0 - c0(i)**2))
            c_rv(i)    = ixidt(i)*temp*(1.d0 - c0(i))**2/(sigma_r(i)*sigma_v(i))

            randy(1,i) = normal(0.0d0,1.0d0)
            randy(2,i) = normal(0.0d0,1.0d0)
            randy(3,i) = normal(0.0d0,1.0d0)

        else                                            ! use series up to 2nd order in xi*dt

            c0(i) = 1.0d0 - xidt(i) + 0.50d0*xidt(i)**2
            c1(i) = (1.0d0 - 0.50d0*xidt(i) + 2.0d0*twelfth*xidt(i)**2)*step
            c2(i) = (0.5d0 - 2.0d0*twelfth*xidt(i) + 0.50d0*twelfth*xidt(i)**2)*step

            sigma_r(i) = step*sqrt(temp*(8.0d0*twelfth*xidt(i) - 0.50d0*xidt(i)**2))
            sigma_v(i) = sqrt(temp*2.0d0*(xidt(i) - xidt(i)**2))
            c_rv(i)    = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt(i))

            randy(1,i) = normal(0.0d0,1.0d0)
            randy(2,i) = normal(0.0d0,1.0d0)
            randy(3,i) = normal(0.0d0,1.0d0)

        end if

    end do

    cofm(1,:) = sigma_r*randy(1,:)
    cofm(2,:) = sigma_r*randy(2,:)
    cofm(3,:) = sigma_r*randy(3,:)

    ! propagate positions
    s%r(1,1:nf) = s%r(1,1:nf) + c1*s%v(1,1:nf) + c2*step*s%a(1,1:nf) + cofm(1,:)
    s%r(2,1:nf) = s%r(2,1:nf) + c1*s%v(2,1:nf) + c2*step*s%a(2,1:nf) + cofm(2,:)
    s%r(3,1:nf) = s%r(3,1:nf) + c1*s%v(3,1:nf) + c2*step*s%a(3,1:nf) + cofm(3,:)

    cofm(1,:) = sigma_v*c_rv*randy(1,:)
    cofm(2,:) = sigma_v*c_rv*randy(2,:)
    cofm(3,:) = sigma_v*c_rv*randy(3,:)

    ! partially propagate velocities
    s%v(1,1:nf) = c0*s%v(1,1:nf) + (c1 - c2)*s%a(1,1:nf) + cofm(1,:)
    s%v(2,1:nf) = c0*s%v(2,1:nf) + (c1 - c2)*s%a(2,1:nf) + cofm(2,:)
    s%v(3,1:nf) = c0*s%v(3,1:nf) + (c1 - c2)*s%a(3,1:nf) + cofm(3,:)

end subroutine langevin_1

subroutine langevin_2(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Dellago et al. JCP 108 (1998) 1964
    !

    type(atoms)                   :: s
    real(8)                       :: imass, temp
    real(8), dimension(  s%nofix) :: c0, c1, c2, xidt, ixidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%nofix) :: randy, cofm
    integer                       :: nf, i

    nf = s%nofix

    temp = kB*Tsurf*imass
    xidt = (s%dens(1:nf)*step)
    ixidt = step/xidt

    do i = 1, nf

        ! Preventing problems due to precision issues
        if (xidt(i) > 1.0d-2) then                      ! use precise expressions

            c0(i) = exp(-xidt(i))
            c1(i) = (1.0d0 - c0(i))*ixidt(i)
            c2(i) = (1.0d0 - c1(i)/step)*ixidt(i)

            sigma_r(i) = ixidt(i)*sqrt(temp*(2.d0*xidt(i) - 3.d0 + 4.d0*c0(i) - c0(i)**2))
            sigma_v(i) = sqrt(temp*(1.d0 - c0(i)**2))
            c_rv(i)    = ixidt(i)*temp*(1.d0 - c0(i))**2/(sigma_r(i)*sigma_v(i))

            randy(1,i) = normal(0.0d0,1.0d0)
            randy(2,i) = normal(0.0d0,1.0d0)
            randy(3,i) = normal(0.0d0,1.0d0)

        else                                            ! use series up to 2nd order in xi*dt

            c0(i) = 1.0d0 - xidt(i) + 0.50d0*xidt(i)**2
            c1(i) = (1.0d0 - 0.50d0*xidt(i) + 2.0d0*twelfth*xidt(i)**2)*step
            c2(i) = (0.5d0 - 2.0d0*twelfth*xidt(i) + 0.50d0*twelfth*xidt(i)**2)*step

            sigma_r(i) = step*sqrt(temp*(8.0d0*twelfth*xidt(i) - 0.50d0*xidt(i)**2))
            sigma_v(i) = sqrt(temp*2.0d0*(xidt(i) - xidt(i)**2))
            c_rv(i)    = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt(i))

            randy(1,i) = normal(0.0d0,1.0d0)
            randy(2,i) = normal(0.0d0,1.0d0)
            randy(3,i) = normal(0.0d0,1.0d0)

        end if

    end do

    cofm(1,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(1,:)
    cofm(2,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(2,:)
    cofm(3,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(3,:)


    ! partially propagate velocities
    s%v(1,1:nf) = s%v(1,1:nf) + c2*s%a(1,1:nf) + cofm(1,:)
    s%v(2,1:nf) = s%v(2,1:nf) + c2*s%a(2,1:nf) + cofm(2,:)
    s%v(3,1:nf) = s%v(3,1:nf) + c2*s%a(3,1:nf) + cofm(3,:)

end subroutine langevin_2


subroutine langevins_1(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !           Li & Wahnström, Phys. Rev. B (1992).
    !

    type(atoms)                   :: s
    real(8)                       :: imass, temp
    real(8), dimension(  s%nofix) :: c0, c1, c2, xidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%nofix) :: randy, cofm
    integer                       :: i
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    xidt = (s%dens(1:nf)*step)
    c0 = 1.0d0 - xidt + 0.50d0*xidt*xidt
    c1 = (1.0d0 - 0.50d0*xidt + 2.0d0*twelfth*xidt*xidt)*step
    c2 = (0.5d0 - 2.0d0*twelfth*xidt + 0.5d0*twelfth*xidt*xidt)*step

    temp = kB*Tsurf*imass
    sigma_r = step*sqrt(temp*(8.0d0*twelfth - 0.5d0*xidt)*xidt)
    sigma_v = sqrt(temp*2.0d0*(1.0d0 - xidt)*xidt)
    c_rv = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt)

    do i =1, nf
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm(1,:) = sigma_r*randy(1,:)
    cofm(2,:) = sigma_r*randy(2,:)
    cofm(3,:) = sigma_r*randy(3,:)

    ! propagate positions
    s%r(1,1:nf) = s%r(1,1:nf) + c1*s%v(1,1:nf) + c2*step*s%a(1,1:nf) + cofm(1,:)
    s%r(2,1:nf) = s%r(2,1:nf) + c1*s%v(2,1:nf) + c2*step*s%a(2,1:nf) + cofm(2,:)
    s%r(3,1:nf) = s%r(3,1:nf) + c1*s%v(3,1:nf) + c2*step*s%a(3,1:nf) + cofm(3,:)

    cofm(1,:) = sigma_v*c_rv*randy(1,:)
    cofm(2,:) = sigma_v*c_rv*randy(2,:)
    cofm(3,:) = sigma_v*c_rv*randy(3,:)

    ! partially propagate velocities
    s%v(1,1:nf) = c0*s%v(1,1:nf) + (c1 - c2)*s%a(1,1:nf) + cofm(1,:)
    s%v(2,1:nf) = c0*s%v(2,1:nf) + (c1 - c2)*s%a(2,1:nf) + cofm(2,:)
    s%v(3,1:nf) = c0*s%v(3,1:nf) + (c1 - c2)*s%a(3,1:nf) + cofm(3,:)


end subroutine langevins_1

subroutine langevins_2(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !           Li & Wahnström, Phys. Rev. B (1992).
    !Ŕ

    type(atoms)                   :: s
    real(8)                       :: imass, temp
    real(8), dimension(  s%nofix) :: c0, c1, c2, xidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%nofix) :: randy, cofm
    integer                       :: i
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    xidt = (s%dens(1:nf)*step)
    c0 = 1.0d0 - xidt + 0.50d0*xidt*xidt
    c1 = (1.0d0 - 0.50d0*xidt + 2.0d0*twelfth*xidt*xidt)*step
    c2 = (0.5d0 - 2.0d0*twelfth*xidt + 0.5d0*twelfth*xidt*xidt)*step

    temp = kB*Tsurf*imass
    sigma_r = step*sqrt(temp*(8.0d0*twelfth - 0.5d0*xidt)*xidt)
    sigma_v = sqrt(temp*2.0d0*(1.0d0 - xidt)*xidt)
    c_rv = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt)

    do i =1, nf
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm(1,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(1,:)
    cofm(2,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(2,:)
    cofm(3,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(3,:)

    ! partially propagate velocities
    s%v(1,1:nf) = s%v(1,1:nf) + c2*s%a(1,1:nf) + cofm(1,:)
    s%v(2,1:nf) = s%v(2,1:nf) + c2*s%a(2,1:nf) + cofm(2,:)
    s%v(3,1:nf) = s%v(3,1:nf) + c2*s%a(3,1:nf) + cofm(3,:)


end subroutine langevins_2

subroutine newton(s, minv)
    !
    ! Purpose:
    !           Newton equation
    !
    type(atoms) :: s
    real(8)     :: minv
    integer     :: nf     ! skip fixed atoms

    nf = s%nofix

    s%a(:,1:nf) = s%f(:,1:nf)*minv

end subroutine newton

end module mdalgo
