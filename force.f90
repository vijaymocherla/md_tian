module force
    !
    ! Purpose:
    !           Gather all potentials to calculate the forces for MD
    !           MAY THE FORCE BE WITH YOU!
    ! Definition of Zero of Energy:
    !                               The particle is at infinite distance to the
    !                               lattice.
    !
    ! Date          Author          History of Revison
    ! ====          ======          ==================
    ! 07.01.2014    Sascha&Svenja   subrutines emt and emt_e implemented
    ! 09.10.2013    Sascha&Svenja   Original
    !
    use atom_class
    use md_init

    implicit none
    save

    real(8), private, parameter                 :: beta    = 1.8093997906d0
    integer, private, dimension(3), parameter   :: b       = (/12, 6, 24/)

contains


subroutine emt(slab, teil)

!   Calculates energy and forces with EMT potential

    implicit none

    type(atoms), intent(inout)    :: teil, slab
    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l, chipl
    real(8) :: betas0_p, betaeta2_p, kappadbeta_p, chilp
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1p, igamma2p, igamma1l, igamma2l
    real(8) :: V_pl, V_lp, V_ll, V_pp, Ecoh_l, Ecoh_p, vref_l, vref_p

    real(8), dimension(3) :: rnnl, rnnp         ! nnn-distances
    real(8), dimension(3) :: xl, xp, r3temp
    real(8), dimension(3) :: dtheta

    real(8), dimension(:), allocatable :: sigma_ll, sigma_lp, s_l
    real(8), dimension(:), allocatable :: sigma_pp, sigma_pl, s_p

    real(8), dimension(:,:), allocatable :: dsigma_lp_l, dsigma_pl_p
    real(8), dimension(:,:), allocatable :: dEcoh_l_l, dEcoh_l_p, dEcoh_p_l, dEcoh_p_p
    real(8), dimension(:,:), allocatable :: dV_ll_l, dV_pp_p, dV_lp_l, dV_lp_p, dV_pl_p, dV_pl_l
    real(8), dimension(:,:), allocatable :: dvref_l_l, dvref_l_p, dvref_p_l, dvref_p_p

    real(8), dimension(:,:,:), allocatable :: dsigma_ll, dsigma_lp_p
    real(8), dimension(:,:,:), allocatable :: dsigma_pp, dsigma_pl_l
    real(8), dimension(:,:,:), allocatable :: ds_l_l, ds_p_p, ds_l_p, ds_p_l

!----------------------VALUES OF FREQUENT USE ---------------------------------

    ! beta * s0
    betas0_l = beta * pars_l(7)
    betas0_p = beta * pars_p(7)
    ! beta * eta2
    betaeta2_l = beta * pars_l(1)
    betaeta2_p = beta * pars_p(1)
    ! kappa / beta
    kappadbeta_l = pars_l(6) / beta
    kappadbeta_p = pars_p(6) / beta

    ! 'coupling' parameters between p and l
    chilp = pars_p(2) / pars_l(2) *exp(0.5d0/bohr2ang*(pars_l(7)-pars_p(7)))
    chipl = 1.0d0 / chilp

    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt2
    rnnp(3) = rnnp(1) * sqrt3

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next-nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4 * rcut / (sqrt3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)

    xl = b * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))
    xp = b * twelfth / (1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.

    r3temp = rnnl - betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(   -pars_l(1) * r3temp))
    igamma2l = 1.0d0 / sum(xl*exp(-kappadbeta_l * r3temp))

    r3temp = rnnp - betas0_p
    igamma1p = 1.0d0 / sum(xp*exp(   -pars_p(1) * r3temp))
    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))

!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(slab%n_atoms), sigma_pp(teil%n_atoms))
    allocate(sigma_lp(slab%n_atoms), sigma_pl(teil%n_atoms))
    allocate(     s_l(slab%n_atoms),      s_p(teil%n_atoms))
    allocate(  dsigma_ll(3, slab%n_atoms, slab%n_atoms))
    allocate(dsigma_lp_l(3, slab%n_atoms))
    allocate(dsigma_lp_p(3, teil%n_atoms, slab%n_atoms))
    allocate(dsigma_pl_l(3, slab%n_atoms, teil%n_atoms))
    allocate(dsigma_pl_p(3, teil%n_atoms))
    allocate(  dsigma_pp(3, teil%n_atoms, teil%n_atoms))
    allocate(     ds_l_l(3, slab%n_atoms, slab%n_atoms))
    allocate(     ds_p_p(3, teil%n_atoms, teil%n_atoms))
    allocate(     ds_l_p(3, teil%n_atoms, slab%n_atoms))
    allocate(     ds_p_l(3, slab%n_atoms, teil%n_atoms))
    allocate(dEcoh_l_l(3,slab%n_atoms), dEcoh_p_p(3, teil%n_atoms))
    allocate(dEcoh_p_l(3,slab%n_atoms), dEcoh_l_p(3, teil%n_atoms))
    allocate(  dV_ll_l(3,slab%n_atoms),   dV_pp_p(3, teil%n_atoms))
    allocate(  dV_lp_l(3,slab%n_atoms),   dV_lp_p(3, teil%n_atoms))
    allocate(  dV_pl_l(3,slab%n_atoms),   dV_pl_p(3, teil%n_atoms))
    allocate(dvref_l_l(3,slab%n_atoms), dvref_p_p(3, teil%n_atoms))
    allocate(dvref_p_l(3,slab%n_atoms), dvref_l_p(3, teil%n_atoms))

    ! initialize accumulators
    sigma_ll    = 0.0d0
    sigma_pp    = 0.0d0
    sigma_pl    = 0.0d0
    sigma_lp    = 0.0d0
    V_ll        = 0.0d0
    V_pp        = 0.0d0
    V_lp        = 0.0d0
    V_pl        = 0.0d0
    vref_l      = 0.0d0
    vref_p      = 0.0d0
    dsigma_ll   = 0.0d0
    dsigma_pp   = 0.0d0
    dsigma_lp_l = 0.0d0
    dsigma_pl_p = 0.0d0
    dV_ll_l     = 0.0d0
    dV_pp_p     = 0.0d0
    dV_lp_l     = 0.0d0
    dV_lp_p     = 0.0d0
    dV_pl_l     = 0.0d0
    dV_pl_p     = 0.0d0
    dvref_l_l   = 0.0d0
    dvref_l_p   = 0.0d0
    dvref_p_l   = 0.0d0
    dvref_p_p   = 0.0d0

    ! slab-slab
    do i = 1, slab%n_atoms
        do j = i+1, slab%n_atoms

            ! Applying PBCs
            r3temp = slab%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            dtheta = (pars_l(1) + rtemp1)*rtemp*r3temp
            dsigma_ll(:,i,i) = dsigma_ll(:,i,i) - dtheta    ! dsigma_i/dr_i
            dsigma_ll(:,j,j) = dsigma_ll(:,j,j) + dtheta
            dsigma_ll(:,j,i) = dtheta                       ! dsigma_i/dr_j
            dsigma_ll(:,i,j) =-dtheta                       ! dsigma_j/dr_i

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

            dtheta = (kappadbeta_l + rtemp1)*rtemp*r3temp
            dV_ll_l(:,i) = dV_ll_l(:,i) + dtheta
            dV_ll_l(:,j) = dV_ll_l(:,j) - dtheta

        end do
    end do
    ! projectile-projectile
    do i = 1, teil%n_atoms
        do j = i+1, teil%n_atoms

            ! Applying PBCs
            r3temp = teil%r(:,i) - teil%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)           ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_pp(i) = sigma_pp(i) + rtemp
            sigma_pp(j) = sigma_pp(j) + rtemp

            dtheta = (pars_p(1) + rtemp1)*rtemp*r3temp
            dsigma_pp(:,i,i) = dsigma_pp(:,i,i) - dtheta    ! dsigma_i/dr_i
            dsigma_pp(:,j,j) = dsigma_pp(:,j,j) + dtheta
            dsigma_pp(:,j,i) = dtheta                       ! dsigma_i/dr_j
            dsigma_pp(:,i,j) =-dtheta                       ! dsigma_j/dr_i

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))   ! V_ij*gamma2*V_0
            V_pp = V_pp + rtemp

            dtheta = (kappadbeta_p + rtemp1)*rtemp*r3temp
            dV_pp_p(:,i) = dV_pp_p(:,i) + dtheta
            dV_pp_p(:,j) = dV_pp_p(:,j) - dtheta

        end do
    end do

    ! projectile-slab
    do i = 1, teil%n_atoms
        do j = 1, slab%n_atoms

            ! Applying PBCs
            r3temp = teil%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)       ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            ! sigma_lp
            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_lp(j) = sigma_lp(j) + rtemp
            dtheta = (pars_p(1) + rtemp1)*rtemp*r3temp
            dsigma_lp_l(:,j) = dsigma_lp_l(:,j) + dtheta ! dsigma_i/dr_i
            dsigma_lp_p(:,i,j) = - dtheta

            ! sigma_pl
            rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )     ! sigma_ij*gamma1


            sigma_pl(i) = sigma_pl(i) + rtemp
            dtheta = (pars_l(1) + rtemp1)*rtemp*r3temp
            dsigma_pl_p(:,i) = dsigma_pl_p(:,i) - dtheta ! dsigma_i/dr_i
            dsigma_pl_l(:,j,i) = dtheta

            ! V_ij*gamma2*V_0
            rtemp = theta*exp(-kappadbeta_p*(r - betas0_p))
            V_lp = V_lp + rtemp
            dtheta = (kappadbeta_p + rtemp1)*rtemp*r3temp
            dV_lp_p(:,i) = dV_lp_p(:,i) + dtheta
            dV_lp_l(:,j) = dV_lp_l(:,j) - dtheta

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l))
            V_pl = V_pl + rtemp
            dtheta = (kappadbeta_l + rtemp1)*rtemp*r3temp
            dV_pl_p(:,i) = dV_pl_p(:,i) + dtheta
            dV_pl_l(:,j) = dV_pl_l(:,j) - dtheta

        end do
    end do

    ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
    V_ll     =     V_ll*igamma2l*pars_l(5)
    sigma_pp = sigma_pp*igamma1p
    V_pp     =     V_pp*igamma2p*pars_p(5)
    sigma_lp = sigma_lp*igamma1l
    sigma_pl = sigma_pl*igamma1p
    V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
    V_pl     =     V_pl*igamma2p*pars_p(5)*chipl

    dsigma_ll   = dsigma_ll  *igamma1l
    dsigma_pp   = dsigma_pp  *igamma1p
    dsigma_lp_l = dsigma_lp_l*igamma1l
    dsigma_lp_p = dsigma_lp_p*igamma1l
    dsigma_pl_p = dsigma_pl_p*igamma1p
    dsigma_pl_l = dsigma_pl_l*igamma1p
    dV_ll_l     =     dV_ll_l*igamma2l*pars_l(5)
    dV_pp_p     =     dV_pp_p*igamma2p*pars_p(5)
    dV_lp_l     =     dV_lp_l*igamma2l*pars_l(5)*chilp
    dV_lp_p     =     dV_lp_p*igamma2l*pars_l(5)*chilp
    dV_pl_l     =     dV_pl_l*igamma2p*pars_p(5)*chipl
    dV_pl_p     =     dV_pl_p*igamma2p*pars_p(5)*chipl

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = sigma_ll + chilp*sigma_lp
    s_p = sigma_pp + chipl*sigma_pl

    ds_l_l =-dsigma_ll
    do i = 1, slab%n_atoms

        ds_l_l(:,i,i) = ds_l_l(:,i,i) - chilp*dsigma_lp_l(:,i)
        ds_l_l(1,i,:) = ds_l_l(1,i,:)/(betaeta2_l*s_l)
        ds_l_l(2,i,:) = ds_l_l(2,i,:)/(betaeta2_l*s_l)
        ds_l_l(3,i,:) = ds_l_l(3,i,:)/(betaeta2_l*s_l)

        ds_p_l(1,i,:) =-chipl*dsigma_pl_l(1,i,:)/(betaeta2_p*s_p)
        ds_p_l(2,i,:) =-chipl*dsigma_pl_l(2,i,:)/(betaeta2_p*s_p)
        ds_p_l(3,i,:) =-chipl*dsigma_pl_l(3,i,:)/(betaeta2_p*s_p)

    end do

    ds_p_p =-dsigma_pp
    do i = 1, teil%n_atoms

        ds_p_p(:,i,i) = ds_p_p(:,i,i) - chipl*dsigma_pl_p(:,i)
        ds_p_p(1,i,:) = ds_p_p(1,i,:)/(betaeta2_p*s_p)
        ds_p_p(2,i,:) = ds_p_p(2,i,:)/(betaeta2_p*s_p)
        ds_p_p(3,i,:) = ds_p_p(3,i,:)/(betaeta2_p*s_p)

        ds_l_p(1,i,:) =-chilp*dsigma_lp_p(1,i,:)/(betaeta2_l*s_l)
        ds_l_p(2,i,:) =-chilp*dsigma_lp_p(2,i,:)/(betaeta2_l*s_l)
        ds_l_p(3,i,:) =-chilp*dsigma_lp_p(3,i,:)/(betaeta2_l*s_l)

    end do

    s_l = -log(s_l*twelfth)/betaeta2_l
    s_p = -log(s_p*twelfth)/betaeta2_p


!----------------------EMBEDDED ELECTRON DENSITY-------------------------------

    rtemp = 0.5d0/bohr2ang - betaeta2_l     ! -eta_l
    slab%dens = pars_l(2)*exp(rtemp*s_l)
    rtemp = 0.5d0/bohr2ang - betaeta2_p     ! -eta_l
    teil%dens = pars_p(2)*exp(rtemp*s_p)


!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

   ! dEcoh_l_l and dEcoh_p_l
    do i = 1, slab%n_atoms

        dEcoh_l_l(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(1,i,:))
        dEcoh_l_l(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(2,i,:))
        dEcoh_l_l(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(3,i,:))

        dEcoh_p_l(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(1,i,:))
        dEcoh_p_l(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(2,i,:))
        dEcoh_p_l(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(3,i,:))

    end do
    dEcoh_l_l = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_l
    dEcoh_p_l = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_l

    ! dEcoh_l_p and dEcoh_p_p
    do i = 1, teil%n_atoms

        dEcoh_l_p(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(1,i,:))
        dEcoh_l_p(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(2,i,:))
        dEcoh_l_p(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(3,i,:))

        dEcoh_p_p(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(1,i,:))
        dEcoh_p_p(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(2,i,:))
        dEcoh_p_p(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(3,i,:))

    end do
    dEcoh_l_p = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_p
    dEcoh_p_p = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_p

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    do i=1,slab%n_atoms

        rtemp = exp(-pars_l(6)*s_l(i))
        vref_l = vref_l + rtemp

        dvref_l_l = dvref_l_l + rtemp*ds_l_l(:,:,i)
        dvref_l_p = dvref_l_p + rtemp*ds_l_p(:,:,i)

    end do

    do i=1,teil%n_atoms

        rtemp = exp(-pars_p(6)*s_p(i))
        vref_p = vref_p + rtemp

        dvref_p_p = dvref_p_p + rtemp*ds_p_p(:,:,i)
        dvref_p_l = dvref_p_l + rtemp*ds_p_l(:,:,i)

    end do

    rtemp = 12.0d0 * pars_l(5)
    vref_l    =    vref_l*rtemp
    dvref_l_l = dvref_l_l*rtemp*pars_l(6)
    dvref_l_p = dvref_l_p*rtemp*pars_l(6)

    rtemp = 12.0d0 * pars_p(5)
    vref_p    =    vref_p*rtemp
    dvref_p_l = dvref_p_l*rtemp*pars_p(6)
    dvref_p_p = dvref_p_p*rtemp*pars_p(6)


!-------------------------------TOTAL ENERGY---------------------------------

    Epot = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

    ! minus sign was taken into account in calculation of separate contributions
    slab%f = dEcoh_l_l + dEcoh_p_l - dV_ll_l &
                 - 0.50d0*(dV_lp_l + dV_pl_l - dvref_l_l - dvref_p_l)
    teil%f = dEcoh_l_p + dEcoh_p_p - dV_pp_p &
                 - 0.50d0*(dV_lp_p + dV_pl_p - dvref_l_p - dvref_p_p)


    deallocate(dvref_l_p, dvref_p_l, dvref_p_p, dvref_l_l)
    deallocate(dV_pl_p, dV_pl_l, dV_lp_p, dV_lp_l, dV_pp_p, dV_ll_l)
    deallocate(dEcoh_p_l, dEcoh_l_p, dEcoh_p_p, dEcoh_l_l)
    deallocate(ds_l_p, ds_p_l, ds_p_p, ds_l_l)
    deallocate(dsigma_pl_p, dsigma_pl_l, dsigma_lp_p, dsigma_lp_l)
    deallocate(dsigma_pp, dsigma_ll)
    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)

end subroutine emt

subroutine emt_e(slab, teil)

! Calculates EMT-energy only

    implicit none

    type(atoms), intent(inout)    :: teil, slab

    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l, chipl
    real(8) :: betas0_p, betaeta2_p, kappadbeta_p, chilp
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1p, igamma2p, igamma1l, igamma2l
    real(8) :: V_pl, V_lp, V_ll, V_pp, Ecoh_l, Ecoh_p, vref_l, vref_p

    real(8), dimension(3) :: rnnl, rnnp         ! nnn-distances
    real(8), dimension(3) :: xl, xp, r3temp

    real(8), dimension(:), allocatable :: sigma_ll, sigma_lp, s_l
    real(8), dimension(:), allocatable :: sigma_pp, sigma_pl, s_p


!----------------------VALUES OF FREQUENT USE ---------------------------------

    ! beta * s0
    betas0_l = beta * pars_l(7)
    betas0_p = beta * pars_p(7)
    ! beta * eta2
    betaeta2_l = beta * pars_l(1)
    betaeta2_p = beta * pars_p(1)
    ! kappa / beta
    kappadbeta_l = pars_l(6) / beta
    kappadbeta_p = pars_p(6) / beta

    ! 'coupling' parameters between p and l
    chilp = pars_p(2) / pars_l(2) *exp(0.5d0/bohr2ang*(pars_l(7)-pars_p(7)))
    chipl = 1.0d0 / chilp

    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt2
    rnnp(3) = rnnp(1) * sqrt3

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next-nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4 * rcut / (sqrt3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)

    xl = b * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))
    xp = b * twelfth / (1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.

    r3temp = rnnl - betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(   -pars_l(1) * r3temp))
    igamma2l = 1.0d0 / sum(xl*exp(-kappadbeta_l * r3temp))

    r3temp = rnnp - betas0_p
    igamma1p = 1.0d0 / sum(xp*exp(   -pars_p(1) * r3temp))
    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))

!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(slab%n_atoms), sigma_pp(teil%n_atoms))
    allocate(sigma_lp(slab%n_atoms), sigma_pl(teil%n_atoms))
    allocate(     s_l(slab%n_atoms),      s_p(teil%n_atoms))

    ! initialize accumulators
    sigma_ll    = 0.0d0
    sigma_pp    = 0.0d0
    sigma_pl    = 0.0d0
    sigma_lp    = 0.0d0
    V_ll        = 0.0d0
    V_pp        = 0.0d0
    V_lp        = 0.0d0
    V_pl        = 0.0d0
    vref_l      = 0.0d0
    vref_p      = 0.0d0

    ! slab-slab
    do i = 1, slab%n_atoms
        do j = i+1, slab%n_atoms

            ! Applying PBCs
            r3temp = slab%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

        end do
    end do

    ! projectile-projectile
    do i = 1, teil%n_atoms
        do j = i+1, teil%n_atoms

            ! Applying PBCs
            r3temp = teil%r(:,i) - teil%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)           ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_pp(i) = sigma_pp(i) + rtemp
            sigma_pp(j) = sigma_pp(j) + rtemp

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))   ! V_ij*gamma2*V_0
            V_pp = V_pp + rtemp

        end do
    end do

    ! projectile-slab
    do i = 1, teil%n_atoms
        do j = 1, slab%n_atoms

            ! Applying PBCs
            r3temp = teil%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)       ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            ! sigma_lp
            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_lp(j) = sigma_lp(j) + rtemp

            ! sigma_pl
            rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )     ! sigma_ij*gamma1
            sigma_pl(i) = sigma_pl(i) + rtemp

            ! V_ij*gamma2*V_0
            rtemp = theta*exp(-kappadbeta_p*(r - betas0_p))
            V_lp = V_lp + rtemp

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l))
            V_pl = V_pl + rtemp

        end do
    end do

    ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
    V_ll     =     V_ll*igamma2l*pars_l(5)
    sigma_pp = sigma_pp*igamma1p
    V_pp     =     V_pp*igamma2p*pars_p(5)
    sigma_lp = sigma_lp*igamma1l
    sigma_pl = sigma_pl*igamma1p
    V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
    V_pl     =     V_pl*igamma2p*pars_p(5)*chipl

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = sigma_ll + chilp*sigma_lp
    s_p = sigma_pp + chipl*sigma_pl
    s_l = -log(s_l*twelfth)/betaeta2_l
    s_p = -log(s_p*twelfth)/betaeta2_p

!----------------------EMBEDDED ELECTRON DENSITY-------------------------------

    rtemp = 0.5d0/bohr2ang - betaeta2_l     ! -eta_l
    slab%dens = pars_l(2)*exp(rtemp*s_l)

    rtemp = 0.5d0/bohr2ang - betaeta2_p     ! -eta_l
    teil%dens = pars_p(2)*exp(rtemp*s_p)

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    do i=1,slab%n_atoms
        rtemp = exp(-pars_l(6)*s_l(i))
        vref_l = vref_l + rtemp
    end do

    do i=1,teil%n_atoms
        rtemp = exp(-pars_p(6)*s_p(i))
        vref_p = vref_p + rtemp
    end do

    rtemp = 12.0d0 * pars_l(5)
    vref_l    =    vref_l*rtemp
    rtemp = 12.0d0 * pars_p(5)
    vref_p    =    vref_p*rtemp


!-------------------------------TOTAL ENERGY---------------------------------

    Epot = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)

end subroutine emt_e

subroutine emt1(s)

!   Calculates energy and forces with EMT potential
!   in case of single species

    implicit none

    type(atoms), intent(inout)    :: s

    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1l, igamma2l
    real(8) :: V_ll, Ecoh_l, vref_l

    real(8), dimension(3) :: rnnl         ! nnn-distances
    real(8), dimension(3) :: xl, r3temp
    real(8), dimension(3) :: dtheta

    real(8), dimension(:), allocatable :: sigma_ll, s_l
    real(8), dimension(:,:), allocatable :: dEcoh_l_l, dV_ll_l, dvref_l_l
    real(8), dimension(:,:,:), allocatable :: dsigma_ll, ds_l_l

!----------------------VALUES OF FREQUENT USE ---------------------------------

    ! beta * s0
    betas0_l = beta * pars_l(7)
    ! beta * eta2
    betaeta2_l = beta * pars_l(1)
    ! kappa / beta
    kappadbeta_l = pars_l(6) / beta

    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next-nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4 * rcut / (sqrt3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)

    xl = b * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.

    r3temp = rnnl - betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(   -pars_l(1) * r3temp))
    igamma2l = 1.0d0 / sum(xl*exp(-kappadbeta_l * r3temp))

!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(s%n_atoms), s_l(s%n_atoms))
    allocate(dsigma_ll(3, s%n_atoms, s%n_atoms), ds_l_l(3, s%n_atoms, s%n_atoms))
    allocate(dEcoh_l_l(3,s%n_atoms), dV_ll_l(3,s%n_atoms), dvref_l_l(3,s%n_atoms))

    ! initialize accumulators
    sigma_ll    = 0.0d0
    V_ll        = 0.0d0
    vref_l      = 0.0d0
    dsigma_ll   = 0.0d0
    dV_ll_l     = 0.0d0
    dvref_l_l   = 0.0d0

    do i = 1, s%n_atoms
        do j = i+1, s%n_atoms

            ! Applying PBCs
            r3temp = s%r(:,i) - s%r(:,j)         ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            dtheta = (pars_l(1) + rtemp1)*rtemp*r3temp
            dsigma_ll(:,i,i) = dsigma_ll(:,i,i) - dtheta    ! dsigma_i/dr_i
            dsigma_ll(:,j,j) = dsigma_ll(:,j,j) + dtheta
            dsigma_ll(:,j,i) = dtheta                       ! dsigma_i/dr_j
            dsigma_ll(:,i,j) =-dtheta                       ! dsigma_j/dr_i

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

            dtheta = (kappadbeta_l + rtemp1)*rtemp*r3temp
            dV_ll_l(:,i) = dV_ll_l(:,i) + dtheta
            dV_ll_l(:,j) = dV_ll_l(:,j) - dtheta

        end do
    end do

    ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
    V_ll     =     V_ll*igamma2l*pars_l(5)

    dsigma_ll   = dsigma_ll  *igamma1l
    dV_ll_l     =     dV_ll_l*igamma2l*pars_l(5)

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = sigma_ll

    ds_l_l =-dsigma_ll
    do i = 1, s%n_atoms
        ds_l_l(1,i,:) = ds_l_l(1,i,:)/(betaeta2_l*s_l)
        ds_l_l(2,i,:) = ds_l_l(2,i,:)/(betaeta2_l*s_l)
        ds_l_l(3,i,:) = ds_l_l(3,i,:)/(betaeta2_l*s_l)
    end do

    s_l = -log(s_l*twelfth)/betaeta2_l

!----------------------EMBEDDED ELECTRON DENSITY-------------------------------

    rtemp = 0.5d0/bohr2ang - betaeta2_l     ! -eta_l
    s%dens = pars_l(2)*exp(rtemp*s_l)

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)

    do i = 1, s%n_atoms
        dEcoh_l_l(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(1,i,:))
        dEcoh_l_l(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(2,i,:))
        dEcoh_l_l(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(3,i,:))
    end do

    dEcoh_l_l = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_l

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    do i=1,s%n_atoms

        rtemp = exp(-pars_l(6)*s_l(i))
        vref_l = vref_l + rtemp

        dvref_l_l = dvref_l_l + rtemp*ds_l_l(:,:,i)

    end do

    rtemp = 12.0d0 * pars_l(5)
    vref_l    =    vref_l*rtemp
    dvref_l_l = dvref_l_l*rtemp*pars_l(6)

!-------------------------------TOTAL ENERGY---------------------------------

    Epot = Ecoh_l - V_ll + 0.50d0*vref_l

    ! minus sign was taken into account in calculation of separate contributions
    s%f = dEcoh_l_l - dV_ll_l + 0.50d0*dvref_l_l

    deallocate(dvref_l_l, dV_ll_l, dEcoh_l_l, ds_l_l, dsigma_ll)
    deallocate(s_l, sigma_ll)

end subroutine emt1

subroutine emt1_e(s)

!   Calculates EMT-energy only
!   in case of single species

    implicit none

    type(atoms), intent(inout)    :: s

    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1l, igamma2l
    real(8) :: V_ll, Ecoh_l, vref_l

    real(8), dimension(3) :: rnnl         ! nnn-distances
    real(8), dimension(3) :: xl, r3temp

    real(8), dimension(:), allocatable :: sigma_ll, s_l

!----------------------VALUES OF FREQUENT USE ---------------------------------

    ! beta * s0
    betas0_l = beta * pars_l(7)
    ! beta * eta2
    betaeta2_l = beta * pars_l(1)
    ! kappa / beta
    kappadbeta_l = pars_l(6) / beta

    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next-nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4 * rcut / (sqrt3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)

    xl = b * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.

    r3temp = rnnl - betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(   -pars_l(1) * r3temp))
    igamma2l = 1.0d0 / sum(xl*exp(-kappadbeta_l * r3temp))

!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(s%n_atoms), s_l(s%n_atoms))

    ! initialize accumulators
    sigma_ll    = 0.0d0
    V_ll        = 0.0d0
    vref_l      = 0.0d0

    do i = 1, s%n_atoms
        do j = i+1, s%n_atoms

            ! Applying PBCs
            r3temp = s%r(:,i) - s%r(:,j)         ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            rtemp1 = acut*rtemp*theta

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

        end do
    end do

    ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
    V_ll     =     V_ll*igamma2l*pars_l(5)

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = -log(sigma_ll*twelfth)/betaeta2_l

!----------------------EMBEDDED ELECTRON DENSITY-------------------------------

    rtemp = 0.5d0/bohr2ang - betaeta2_l     ! -eta_l
    s%dens = pars_l(2)*exp(rtemp*s_l)

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    vref_l = 12.0d0*pars_l(5)*sum(exp(-pars_l(6)*s_l))


!-------------------------------TOTAL ENERGY---------------------------------

    Epot = Ecoh_l - V_ll + 0.50d0*vref_l


    deallocate(s_l, sigma_ll)

end subroutine emt1_e

subroutine num_emt(s,t)
    !
    ! Purpose:
    !          Calculate numerical forces

    type(atoms) :: s, t
    integer :: i,j
    real(8) :: delta = 0.001d0, Epoti

    Epoti = Epot
    do j=1, 3
    do i=1, s%n_atoms
        s%r(j,i) = s%r(j,i) - delta
        call emt_e(s, t)
        s%f(j,i) = Epot/(2*delta)
        s%r(j,i) = s%r(j,i) + 2*delta
        call emt_e(s, t)
        s%r(j,i) = s%r(j,i) - delta
        s%f(j,i) =  s%f(j,i) - Epot/(2*delta)
    end do
    end do

    Epot = Epoti

end subroutine num_emt

subroutine num_emt1(s)
    !
    ! Purpose:
    !          Calculate numerical forces

    type(atoms) :: s
    integer :: i,j
    real(8) :: delta = 0.001d0, Epoti

    Epoti=Epot
    do j=1, 3
    do i=1, s%n_atoms
        s%r(j,i) = s%r(j,i) - delta
        call emt1_e(s)
        s%f(j,i) = Epot/(2*delta)
        s%r(j,i) = s%r(j,i) + 2*delta
        call emt1_e(s)
        s%r(j,i) = s%r(j,i) - delta
        s%f(j,i) =  s%f(j,i) - Epot/(2*delta)
    end do
    end do
    Epot=Epoti

end subroutine num_emt1



subroutine ldfa(s)
    !
    ! Purpose:
    !           Calculate the friction coefficient
    !
    type(atoms) :: s
    integer :: i, j
    real(8) :: fric, temp
    real(8), dimension(12), parameter :: coefs = (/ 0.0802484d0,-1.12851d0,&
                                                      9.28508d0, 2.10064d0,&
                                                     -843.419d0, 8.85354d3,&
                                                     -4.89023d4,  1.6741d5,&
                                                     -3.67098d5, 5.03476d5,&
                                                      -3.9426d5, 1.34763d5  /)

    !	12th order cubic spline fit interpolated from DFT data points of friction
    !   coefficient vs. electron density (calculated from DFT with VASP)

	do i = 1, s%n_atoms
        fric = s%dens(i)
        ! hbar*xi in eV
        if (fric >= -1.d-12 .or. fric <= 0.36) then  ! removed offset parameter
            fric = 0.0d0
            temp = s%dens(i)
            do j=1,12
                fric = fric + coefs(j)*temp
                temp = temp*s%dens(i)
            end do
        else if (fric > 0.36) then
            fric = 0.001*(4.7131-exp(-4.41305*fric))
        else
            print *, fric, s%dens(i), i, shape(s%dens)
            stop 'Error: dens takes negative values!'
        end if
        s%dens(i)=fric
    end do
    ! xi in 1/fs
    s%dens = s%dens / hbar

end subroutine ldfa


end module force

