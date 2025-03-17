module force
    !
    ! Purpose:
    !           Calculates energy and forces
    !           MAY THE FORCE BE WITH YOU!
    !
    ! Implemented Force Fields:
    !				EMT (effective medium theory)
    ! Definition of Zero of Energy:
    !                               A projectile is at infinite distance from the surface
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 18.02.2014    	Svenja M. Janke		Original
    !			Sascha Kandratsenka
    !			Dan J. Auerbach
    !
    use atom_class
    use md_init

    implicit none
    save
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
    dsigma_pl_l = 0.0d0
    dsigma_lp_p = 0.0d0
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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1, temp
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
        temp = V_ll
        do j = i+1, s%n_atoms

            ! Applying PBCs
            r3temp = s%r(:,i) - s%r(:,j)         ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance

            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle
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
            close(124)
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

subroutine emt_e_fit(xdata, energy)

! Calculates EMT-energy only

    implicit none

    real(8),dimension(:,:), intent(in)    :: xdata
    real(8)                               :: energy


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
    !a_lat = 4.2010d0
    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4.0d0 * rcut / (sqrt3 + 2.0d0)
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

    allocate(sigma_ll(nl_atoms), sigma_pp(np_atoms))
    allocate(sigma_lp(nl_atoms), sigma_pl(np_atoms))
    allocate(     s_l(nl_atoms),      s_p(np_atoms))

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
    do i = 1, nl_atoms
        do j = i+1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i+np_atoms) - xdata(:,j+np_atoms)   ! distance vector
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
    do i = 1, np_atoms
        do j = i+1, np_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j)   ! distance vector
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
    do i = 1, np_atoms
        do j = 1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j+np_atoms)   ! distance vector
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

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    do i=1,nl_atoms
        rtemp = exp(-pars_l(6)*s_l(i))
        vref_l = vref_l + rtemp
    end do

    do i=1,np_atoms
        rtemp = exp(-pars_p(6)*s_p(i))
        vref_p = vref_p + rtemp
    end do

    rtemp = 12.0d0 * pars_l(5)
    vref_l    =    vref_l*rtemp
    rtemp = 12.0d0 * pars_p(5)
    vref_p    =    vref_p*rtemp


!-------------------------------TOTAL ENERGY---------------------------------

    energy = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)


    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)

end subroutine emt_e_fit

subroutine emt_de_fit(xdata, energy, denergy)

! Calculates EMT-energy only

    implicit none

    real(8),dimension(:,:), intent(in)   :: xdata
    real(8),                intent(out)  :: energy
    real(8), dimension(14), intent(out)  :: denergy  ! derivatives
                                                        ! with respect to to
                                                        ! eta2, followed by no,
                                                        ! eo, lambda, vo, kappa
                                                        ! and so.

    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l, chipl
    real(8) :: betas0_p, betaeta2_p, kappadbeta_p, chilp
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1p, igamma2p, igamma1l, igamma2l
    real(8) :: V_pl, V_lp, V_ll, V_pp, Ecoh_l, Ecoh_p, vref_l, vref_p

    real(8), dimension(3) :: rnnl, rnnp         ! nnn-distances
    real(8), dimension(3) :: xl, xp, r3temp, r3temp1

    real(8), dimension(:), allocatable :: sigma_ll, sigma_lp, s_l, rn_ltemp
    real(8), dimension(:), allocatable :: sigma_pp, sigma_pl, s_p, rn_ptemp

!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
! Variables and Arrays for partial derivatives
! Apart from chi, all derivatives are 7 long. The first place denotes the
! derivative with respect to to eta2, followed by no, eo, lambda, vo, kappa and so.
! The general notation is:
! e.g. dV_lp_l(1) : the derivative of V_lp with respect to eta2_l.
    real(8), dimension(4) :: dchilp= 0.0d0, dchipl= 0.0d0     ! First element: p, then l
    real(8), dimension(7) :: dgamma1l= 0.0d0, dgamma2l= 0.0d0
    real(8), dimension(7) :: dgamma1p= 0.0d0, dgamma2p= 0.0d0
    real(8), dimension(7,nl_atoms) :: dsigma_ll
    real(8), dimension(7,np_atoms) :: dsigma_pp
    real(8), dimension(7,nl_atoms) :: dsigma_lp_l
    real(8), dimension(7,nl_atoms) :: dsigma_lp_p
    real(8), dimension(7, np_atoms) :: dsigma_pl_l, dsigma_pl_p
    real(8), dimension(7) :: dV_ll, dV_pp
    real(8), dimension(7) :: dV_lp_l, dV_lp_p
    real(8), dimension(7) :: dV_pl_l, dV_pl_p
    real(8), dimension(7,nl_atoms) :: ds_l_l, ds_l_p
    real(8), dimension(7,np_atoms) :: ds_p_l, ds_p_p
    real(8), dimension(7) :: dvref_l_l, dvref_l_p
    real(8), dimension(7) :: dvref_p_l, dvref_p_p
    real(8), dimension(7) :: dEcoh_l, dEcoh_p
    real(8), dimension(3) :: drnn, dxl, dxp


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
    rtemp = exp(0.5d0/bohr2ang*(pars_l(7)-pars_p(7)))
    chilp = pars_p(2) / pars_l(2) * rtemp
    chipl = 1.0d0 / chilp

! derivatives: (1) derivative over nop, (2) over nol
!
    dchilp(1) =  chilp/pars_p(2)        ! d chilp / d nop
    dchilp(2) = -chilp/pars_l(2)        ! d chilp / d nol
    dchilp(3) = -chilp * 0.5d0/bohr2ang ! d chilp / d sop
    dchilp(4) = -dchilp(3)              ! d chilp / d sol

    dchipl(1) = -chipl/pars_p(2)        ! d chipl / d nop
    dchipl(2) =  chipl/pars_l(2)        ! d chipl / d nol
    dchipl(3) =  chipl * 0.5d0/bohr2ang ! d chipl / d sop
    dchipl(4) = -dchipl(3)              ! d chipl / d sol
    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt2
    rnnp(3) = rnnp(1) * sqrt3

    drnn(1) = beta
    drnn(2) = drnn(1) * sqrt2
    drnn(3) = drnn(1) * sqrt3

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
    dxl = -xl*acut*drnn*exp(acut*(rnnl-rcut))/(1.0d0 + exp(acut*(rnnl-rcut)))
    dxp = -xp*acut*drnn*exp(acut*(rnnp-rcut))/(1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
! The derivative is not defined as the inverse and formed for each gamma
! individually.

    r3temp = rnnl - betas0_l

    r3temp1 = xl*exp(- pars_l(1)*r3temp)
    igamma1l = 1.0d0 / sum(r3temp1)
    dgamma1l = 0.d0
    dgamma1l(1) = - sum(r3temp*r3temp1)
    dgamma1l(7) = sum(r3temp1*betaeta2_l)+sum(dxl*exp(- pars_l(1)*r3temp)) -&
                  sum(r3temp1*pars_l(1)*drnn)

    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
    igamma2l = 1.0d0 / sum(r3temp1)
    dgamma2l = 0.d0
    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
    dgamma2l(7) = sum(r3temp1 * pars_l(6))
    dgamma2l(7) = sum(r3temp1*pars_l(6))+sum(dxl*exp(-kappadbeta_l*r3temp)) -&
                  sum(r3temp1*kappadbeta_l*drnn)

    r3temp = rnnp-betas0_p

    r3temp1 = xp*exp(- pars_p(1)*r3temp)
    igamma1p = 1.0d0 / sum(r3temp1)
    dgamma1p = 0.0d0
    dgamma1p(1) = - sum(r3temp*r3temp1)
    dgamma1p(7) = sum(r3temp1*betaeta2_p)+sum(dxp*exp(- pars_p(1)*r3temp)) -&
                  sum(r3temp1*pars_p(1)*drnn)

    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
    igamma2p = 1.0d0 / sum(r3temp1)
    dgamma2p = 0.d0
    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
    dgamma2p(7) = sum(r3temp1*pars_p(6))+sum(dxp*exp(-kappadbeta_p*r3temp)) -&
                  sum(r3temp1*kappadbeta_p*drnn)


!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(nl_atoms), sigma_pp(np_atoms))
    allocate(sigma_lp(nl_atoms), sigma_pl(np_atoms))
    allocate(     s_l(nl_atoms),      s_p(np_atoms))
    allocate(rn_ltemp(nl_atoms), rn_ptemp(np_atoms))

    ! initialize accumulators
    sigma_ll = 0.0d0
    sigma_pp = 0.0d0
    dsigma_ll = 0.0d0
    dsigma_pp = 0.0d0
    sigma_pl = 0.0d0
    sigma_lp=0.0d0
    dsigma_lp_l=0.0d0
    dsigma_lp_p=0.0d0
    dsigma_pl_l=0.0d0
    dsigma_pl_p=0.0d0
    V_ll = 0.0d0
    V_pp = 0.0d0
    dV_ll = 0.0d0
    dV_pp = 0.0d0
    V_lp = 0.0d0
    V_pl = 0.0d0
    dV_lp_l = 0.0d0
    dV_lp_p = 0.0d0
    dV_pl_l = 0.0d0
    dV_pl_p = 0.0d0
    ds_l_l = 0.0d0
    ds_l_p = 0.0d0
    ds_p_l = 0.0d0
    ds_p_p = 0.0d0
    dvref_l_l=0.0d0
    dvref_l_p=0.0d0
    dvref_p_l=0.0d0
    dvref_p_p=0.0d0
    dEcoh_l = 0.0d0
    dEcoh_p = 0.0d0
    denergy = 0.0d0
    energy = 0.0d0
    rn_ltemp = 0.0d0
    rn_ptemp=0.0d0
    s_l = 0.0d0
    s_p = 0.0d0

    ! slab-slab
    do i = 1, nl_atoms
        do j = i+1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i+np_atoms) - xdata(:,j+np_atoms)   ! distance vector
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

            !----------------------------SIGMA LATTICE-----------------------------
            ! Sigma is a contribution to the neutral sphere radius.
            ! It is a list in which for each lattice atom, the contributions of the
            ! others are summed up. To enforce the cut-off, it will be later
            ! corrected by gamma.

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            rtemp1 = rtemp*(r - betas0_l)
            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1


            !-----------------------PAIR POTENTIAL LATTICE-------------------------
            ! Will later be subjected to gamma to complete the cut-off.

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

            rtemp1 = rtemp*(r - betas0_l)
            dV_ll(6) = dV_ll(6) + rtemp1

        end do
    end do

    ! projectile-projectile
    do i = 1, np_atoms
        do j = i+1, np_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j)   ! distance vector
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

            !----------------------------SIGMA PARTICLE-----------------------------
            ! Sigma is a contribution to the neutral sphere radius.
            ! It is a list in which for each particle atom, the contributions of the
            ! others are summed up. To enforce the cut-off, it will be later
            ! corrected by gamma.

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_pp(i) = sigma_pp(i) + rtemp
            sigma_pp(j) = sigma_pp(j) + rtemp

            rtemp1 = rtemp*(r - betas0_p)
            dsigma_pp(1,i) = dsigma_pp(1,i) - rtemp1
            dsigma_pp(1,j) = dsigma_pp(1,j) - rtemp1

            rtemp1 = rtemp*betaeta2_p
            dsigma_pp(7,i) = dsigma_pp(7,i) + rtemp1
            dsigma_pp(7,j) = dsigma_pp(7,j) + rtemp1

            !-----------------------PAIR POTENTIAL LATTICE-------------------------
            ! Will later be subjected to gamma to complete the cut-off.

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))   ! V_ij*gamma2*V_0
            V_pp = V_pp + rtemp

            rtemp1 = rtemp*(r - betas0_p)
            dV_pp(6) = dV_pp(6) + rtemp1
            dV_pp(7) = dV_pp(7) + rtemp*pars_p(6)

       end do
    end do

    ! projectile-slab
    do i = 1, np_atoms
        do j = 1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j+np_atoms)   ! distance vector
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

            !-------------------------------MIXED SIGMA--------------------------------
            ! Contributions of both particle and lattice to neutral sphere radius
            ! To fully include the cut-off, we correct them later by gamma.
            ! Each of the mixed sigmas depends on both l and p parameters. Not all
            ! contributions need to be under the loop.
            ! sigma_lp

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_lp(j) = sigma_lp(j) + rtemp
            dsigma_lp_p(1,j) = dsigma_lp_p(1,j)-(r-betas0_p)*rtemp
            dsigma_lp_l(7,i) = dsigma_lp_l(7,i)+ betaeta2_p*rtemp

            ! sigma_pl
            rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )     ! sigma_ij*gamma1
            sigma_pl(i) = sigma_pl(i) + rtemp

            dsigma_pl_l(1,i) = dsigma_pl_l(1,i) - (r - betas0_l)*rtemp
            dsigma_pl_p(7,j) = dsigma_pl_p(7,j) + rtemp*betaeta2_l

            !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

            ! V_ij*gamma2*V_0
            rtemp = theta*exp(-kappadbeta_p*(r - betas0_p))
            V_lp = V_lp + rtemp

            dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
            dV_lp_p(7) = dV_lp_p(7) + pars_p(6)*rtemp

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l))
            V_pl = V_pl + rtemp
            dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
            dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l(6)

        end do
    end do


   ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
        dsigma_ll(7,:) = sigma_ll*pars_l(1)*beta-sigma_ll*igamma1l*dgamma1l(7)

    V_ll     =     V_ll*igamma2l*pars_l(5)
        dV_ll(5) = - V_ll/pars_l(5)
        dV_ll(6) = (dV_ll(6) * pars_l(5)/beta + V_ll*dgamma2l(6)) * igamma2l
        dV_ll(7) = V_ll*pars_l(6)-V_ll*igamma2l*dgamma2l(7)

    sigma_pp = sigma_pp*igamma1p
        dsigma_pp(1,:) = (dsigma_pp(1,:) - sigma_pp*dgamma1p(1))*igamma1p
        dsigma_pp(7,:) = (dsigma_pp(7,:) - sigma_pp*dgamma1p(7))*igamma1p

    V_pp     =     V_pp*igamma2p*pars_p(5)
        dV_pp(5) = - V_pp/pars_p(5)
        dV_pp(6) = (dV_pp(6) * pars_p(5)/beta + V_pp*dgamma2p(6)) * igamma2p
        dV_pp(7) = (pars_p(5)*dV_pp(7)-V_pp*dgamma2p(7))*igamma2p

    sigma_lp = sigma_lp*igamma1l
        ! Derivative with respect to l
        dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
        dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
        ! Derivative with respect to p
        dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
        dsigma_lp_p(7,:) = sigma_lp*pars_p(1)*beta

    sigma_pl = sigma_pl*igamma1p
        ! Derivative with respect to l
        dsigma_pl_l(1,:) = dsigma_pl_l(1,:)*igamma1p
        dsigma_pl_l(7,:) = sigma_pl*pars_l(1)*beta

        ! Derivative with respect to p
        dsigma_pl_p(1,:) = - sigma_pl*igamma1p*dgamma1p(1)
        dsigma_pl_p(7,:) = - sigma_pl*igamma1p*dgamma1p(7)

    V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
        ! Derivative with respect to l
        dV_lp_l(2) = V_lp/pars_l(2)
        dV_lp_l(5) = - V_lp/pars_l(5)
        dV_lp_l(7) = V_lp*igamma2l
        dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
        dV_lp_l(7) = V_lp*dchilp(4)*chipl - dV_lp_l(7)*dgamma2l(7)
        ! Derivative with respect to p
        dV_lp_p(2) = -V_lp/pars_p(2)
        dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l(5) / beta
        dV_lp_p(7) = V_lp*chipl*dchilp(3) + dV_lp_p(7)*chilp*pars_l(5)*igamma2l

    V_pl     =     V_pl*igamma2p*pars_p(5)*chipl
        ! Derivative with respect to l
        dV_pl_l(2) = - V_pl/pars_l(2)
        dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p(5)*chipl/beta
        dV_pl_l(7) = V_pl*chilp*dchipl(4) + dV_pl_l(7)*pars_p(5)*chipl*igamma2p
        ! Derivative with respect to p
        dV_pl_p(2) = V_pl/pars_p(2)
        dV_pl_p(5) = - V_pl/pars_p(5)
        dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
        dV_pl_p(7) = V_pl*chilp*dchipl(3) - V_pl*igamma2p*dgamma2p(7)

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = sigma_ll + chilp*sigma_lp
    rn_ltemp = 1.0d0/(s_l*betaeta2_l)
    s_l = -log(s_l*twelfth)/betaeta2_l

    ! Derivative with respect to l
    ds_l_l(1,:) = -s_l/pars_l(1) &
                  - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
    ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
    ds_l_l(7,:) = rn_ltemp*(dsigma_ll(7,:)+dchilp(4)*sigma_lp+chilp*dsigma_lp_l(7,:))
    ! Derivative with respect to p
    ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
    ds_l_p(2,:) = - sigma_lp*rn_ltemp*dchilp(1)
    ds_l_p(7,:) =  rn_ltemp*(dchilp(3)*sigma_lp+chilp*dsigma_lp_p(7,:))

    s_p = sigma_pp + chipl*sigma_pl
    rn_ptemp = 1.0d0 / (s_p*betaeta2_p)
    s_p = -log(s_p*twelfth)/betaeta2_p

    ! Derivative with respect to p
    ds_p_p(1,:) = -s_p/pars_p(1) &
                  - (dsigma_pp(1,:)+chipl*dsigma_pl_p(1,:))*rn_ptemp
    ds_p_p(2,:) = - sigma_pl*rn_ptemp*dchipl(1)
    ds_p_p(7,:) = rn_ptemp*(dsigma_pp(7,:)+dchipl(3)*sigma_pl+chipl*dsigma_pl_p(7,:))
    ! Derivative with respect to p
    ds_p_l(1,:) = - rn_ptemp*chipl*dsigma_pl_l(1,:)
    ds_p_l(2,:) = - sigma_pl*rn_ptemp*dchipl(2)
    ds_p_l(7,:) =   rn_ptemp*(dchipl(4)*sigma_pl+chipl*dsigma_pl_l(7,:))

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

    rn_ltemp = -pars_l(3)*pars_l(4)*s_l*exp(-pars_l(4)*s_l)
    rn_ptemp = -pars_p(4)*s_p*pars_p(3)*exp(-pars_p(4)*s_p)
        ! Derivative with respect to l
        dEcoh_l(1) =  sum(pars_l(4)*rn_ltemp*ds_l_l(1,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_l(1,:))
        dEcoh_l(2) =  sum(pars_l(4)*rn_ltemp*ds_l_l(2,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_l(2,:))
        dEcoh_l(3) =  sum((1.0d0 + pars_l(4)*s_l) * exp(-pars_l(4) * s_l)-1.0d0)
        dEcoh_l(4) =  sum(s_l*rn_ltemp)
        dEcoh_l(7) =  sum(pars_l(4)*rn_ltemp*ds_l_l(7,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_l(7,:))

        ! Derivative with respect to p
        dEcoh_p(1) =  sum(pars_l(4)*rn_ltemp*ds_l_p(1,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_p(1,:))
        dEcoh_p(2) =  sum(pars_l(4)*rn_ltemp*ds_l_p(2,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_p(2,:))
        dEcoh_p(3) =  sum((1.0d0 + pars_p(4)*s_p) * exp(-pars_p(4) * s_p)-1.0d0)
        dEcoh_p(4) =  sum(s_p*rn_ptemp)
        dEcoh_p(7) =  sum(pars_l(4)*rn_ltemp*ds_l_p(7,:))&
                    + sum(pars_p(4)*rn_ptemp*ds_p_p(7,:))

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    rn_ltemp = exp( -pars_l(6) * s_l)
    rtemp = -12.0d0 * pars_l(5) * pars_l(6)
    vref_l = 12.0d0 * pars_l(5) * sum(rn_ltemp)
    ! Derivative with respect to l
        dvref_l_l(1) = rtemp*sum(rn_ltemp*ds_l_l(1,:))
        dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l(2,:))
        dvref_l_l(5) = vref_l/pars_l(5)
        dvref_l_l(6) = - 12.0d0 * pars_l(5) * sum(rn_ltemp *s_l)
        dvref_l_l(7) = rtemp*sum(rn_ltemp*ds_l_l(7,:))
        ! Derivative with respect to p
        dvref_l_p(1) = rtemp*sum(rn_ltemp*ds_l_p(1,:))
        dvref_l_p(2) = rtemp*sum(rn_ltemp*ds_l_p(2,:))
        dvref_l_p(7) = rtemp*sum(rn_ltemp*ds_l_p(7,:))

    rn_ptemp = exp( -pars_p(6) * s_p)
    rtemp = -12.0d0 * pars_p(5) * pars_p(6)
    vref_p = 12.0d0 * pars_p(5) * sum(rn_ptemp)
    ! Derivative with respect to l
        dvref_p_p(1) = rtemp*sum(rn_ptemp*ds_p_p(1,:))
        dvref_p_p(2) = rtemp*sum(rn_ptemp*ds_p_p(2,:))
        dvref_p_p(5) = vref_p/pars_p(5)
        dvref_p_p(6) = - 12.0d0 * pars_p(5) * sum(rn_ptemp *s_p)
        dvref_p_p(7) = rtemp*sum(rn_ptemp*ds_p_p(7,:))
        ! Derivative with respect to p
        dvref_p_l(1) = rtemp*sum(rn_ptemp*ds_p_l(1,:))
        dvref_p_l(2) = rtemp*sum(rn_ptemp*ds_p_l(2,:))
        dvref_p_l(7) = rtemp*sum(rn_ptemp*ds_p_l(7,:))

!-------------------------------TOTAL ENERGY---------------------------------

    energy = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

    ! Derivative with respect to l
    denergy(8) = dEcoh_l(1) + 0.5d0*( dvref_l_l(1)+dvref_p_l(1))
    denergy(9) = dEcoh_l(2) &
                    +0.5d0*(dV_pl_l(2)+dvref_p_l(2)+dV_lp_l(2)+dvref_l_l(2))
    denergy(10) = dEcoh_l(3)
    denergy(11) = dEcoh_l(4)
    denergy(12) = dV_ll(5) + 0.5d0*(dvref_l_l(5)+dV_lp_l(5))
    denergy(13) = dV_ll(6) + 0.5d0*( dV_lp_l(6) + dV_pl_l(6) + dvref_l_l(6))
    denergy(14) = dEcoh_l(7) + dV_ll(7) + &
                  0.5d0*(dV_lp_l(7) + dV_pl_l(7) + dvref_l_l(7) + dvref_p_l(7))
    ! Derivative with respect to p (no correction by dEref since those do not
    ! contain any p-contribution)
    denergy(1) = dEcoh_p(1) + 0.5d0*(dvref_l_p(1)+dvref_p_p(1))
    denergy(2) = dEcoh_p(2) &
                   + 0.5d0*(dV_pl_p(2)+dV_lp_p(2)+dvref_l_p(2)+dvref_p_p(2))
    denergy(3) = dEcoh_p(3)
    denergy(4) = dEcoh_p(4)
    denergy(5) = dV_pp(5)+0.5d0*(dV_pl_p(5) + dvref_p_p(5))
    denergy(6) = dV_pp(6)+0.5d0*(dV_lp_p(6)+dV_pl_p(6)+dvref_p_p(6))
    denergy(7) = dEcoh_p(7) + dV_pp(7) + &
                 0.5d0*(dV_lp_p(7) + dV_pl_p(7) + dvref_l_p(7) + dvref_p_p(7))

    deallocate(rn_ltemp, rn_ptemp)
    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)


end subroutine emt_de_fit

subroutine emt_dens_fit(xdata, energy,pdens)

! Calculates EMT-energy only

    implicit none

    real(8),dimension(:,:), intent(in)    :: xdata
    real(8)                               :: energy, pdens


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
    !real(8), dimension(:), allocatable :: pdens


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
    !a_lat = 4.2010d0
    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4.0d0 * rcut / (sqrt3 + 2.0d0)
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

    allocate(sigma_ll(nl_atoms), sigma_pp(np_atoms))
    allocate(sigma_lp(nl_atoms), sigma_pl(np_atoms))
    allocate(     s_l(nl_atoms),      s_p(np_atoms))
!    allocate(                       pdens(np_atoms))

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
    do i = 1, nl_atoms
        do j = i+1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i+np_atoms) - xdata(:,j+np_atoms)   ! distance vector
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
    do i = 1, np_atoms
        do j = i+1, np_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j)   ! distance vector
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
    do i = 1, np_atoms
        do j = 1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j+np_atoms)   ! distance vector
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

!-----------------------------DENSITY CALCULATION------------------------------
! only for particle at present
    rtemp = 0.5d0/bohr2ang - betaeta2_p     ! -eta_l
    pdens = pars_p(2)*exp(rtemp*s_p(5))


!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    do i=1,nl_atoms
        rtemp = exp(-pars_l(6)*s_l(i))
        vref_l = vref_l + rtemp
    end do

    do i=1,np_atoms
        rtemp = exp(-pars_p(6)*s_p(i))
        vref_p = vref_p + rtemp
    end do

    rtemp = 12.0d0 * pars_l(5)
    vref_l    =    vref_l*rtemp
    rtemp = 12.0d0 * pars_p(5)
    vref_p    =    vref_p*rtemp


!-------------------------------TOTAL ENERGY---------------------------------

    energy = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)

end subroutine emt_dens_fit

subroutine emt_ddens_fit(xdata, energy, denergy)

! Calculates EMT-energy only

    implicit none

    real(8),dimension(:,:), intent(in)   :: xdata
    real(8),                intent(out)  :: energy
    real(8), dimension(14), intent(out)  :: denergy  ! derivatives
                                                        ! with respect to to
                                                        ! eta2, followed by no,
                                                        ! eo, lambda, vo, kappa
                                                        ! and so.

    integer :: i,j

    real(8) :: betas0_l, betaeta2_l, kappadbeta_l, chipl
    real(8) :: betas0_p, betaeta2_p, kappadbeta_p, chilp
    real(8) :: r, rcut, rr, acut, theta, rtemp, rtemp1
    real(8) :: igamma1p, igamma2p, igamma1l, igamma2l
    real(8) :: V_pl, V_lp, V_ll, V_pp

    real(8), dimension(3) :: rnnl, rnnp         ! nnn-distances
    real(8), dimension(3) :: xl, xp, r3temp, r3temp1

    real(8), dimension(:), allocatable :: sigma_ll, sigma_lp, s_l, rn_ltemp
    real(8), dimension(:), allocatable :: sigma_pp, sigma_pl, s_p, rn_ptemp

!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
! Variables and Arrays for partial derivatives
! Apart from chi, all derivatives are 7 long. The first place denotes the
! derivative with respect to to eta2, followed by no, eo, lambda, vo, kappa and so.
! The general notation is:
! e.g. dV_lp_l(1) : the derivative of V_lp with respect to eta2_l.
    real(8), dimension(4) :: dchilp= 0.0d0, dchipl= 0.0d0     ! First element: p, then l
    real(8), dimension(7) :: dgamma1l= 0.0d0, dgamma2l= 0.0d0
    real(8), dimension(7) :: dgamma1p= 0.0d0, dgamma2p= 0.0d0
    real(8), dimension(7,nl_atoms) :: dsigma_ll
    real(8), dimension(7,np_atoms) :: dsigma_pp
    real(8), dimension(7,nl_atoms) :: dsigma_lp_l
    real(8), dimension(7,nl_atoms) :: dsigma_lp_p
    real(8), dimension(7, np_atoms) :: dsigma_pl_l, dsigma_pl_p
    real(8), dimension(7) :: dV_ll, dV_pp
    real(8), dimension(7) :: dV_lp_l, dV_lp_p
    real(8), dimension(7) :: dV_pl_l, dV_pl_p
    real(8), dimension(7,nl_atoms) :: ds_l_l, ds_l_p
    real(8), dimension(7,np_atoms) :: ds_p_l, ds_p_p
    real(8), dimension(7) :: dvref_l_l, dvref_l_p
    real(8), dimension(7) :: dvref_p_l, dvref_p_p
    real(8), dimension(7) :: dEcoh_l, dEcoh_p
    real(8), dimension(3) :: drnn, dxl, dxp


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
    rtemp = exp(0.5d0/bohr2ang*(pars_l(7)-pars_p(7)))
    chilp = pars_p(2) / pars_l(2) * rtemp
    chipl = 1.0d0 / chilp

! derivatives: (1) derivative over nop, (2) over nol
!
    dchilp(1) =  chilp/pars_p(2)        ! d chilp / d nop
    dchilp(2) = -chilp/pars_l(2)        ! d chilp / d nol
    dchilp(3) = -chilp * 0.5d0/bohr2ang ! d chilp / d sop
    dchilp(4) = -dchilp(3)              ! d chilp / d sol

    dchipl(1) = -chipl/pars_p(2)        ! d chipl / d nop
    dchipl(2) =  chipl/pars_l(2)        ! d chipl / d nol
    dchipl(3) =  chipl * 0.5d0/bohr2ang ! d chipl / d sop
    dchipl(4) = -dchipl(3)              ! d chipl / d sol
    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt2
    rnnp(3) = rnnp(1) * sqrt3

    drnn(1) = beta
    drnn(2) = drnn(1) * sqrt2
    drnn(3) = drnn(1) * sqrt3

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
    dxl = -xl*acut*drnn*exp(acut*(rnnl-rcut))/(1.0d0 + exp(acut*(rnnl-rcut)))
    dxp = -xp*acut*drnn*exp(acut*(rnnp-rcut))/(1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
! The derivative is not defined as the inverse and formed for each gamma
! individually.

    r3temp = rnnl - betas0_l

    r3temp1 = xl*exp(- pars_l(1)*r3temp)
    igamma1l = 1.0d0 / sum(r3temp1)
    dgamma1l = 0.d0
    dgamma1l(1) = - sum(r3temp*r3temp1)
    dgamma1l(7) = sum(r3temp1*betaeta2_l)+sum(dxl*exp(- pars_l(1)*r3temp)) -&
                  sum(r3temp1*pars_l(1)*drnn)

    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
    igamma2l = 1.0d0 / sum(r3temp1)
    dgamma2l = 0.d0
    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
    dgamma2l(7) = sum(r3temp1 * pars_l(6))
    dgamma2l(7) = sum(r3temp1*pars_l(6))+sum(dxl*exp(-kappadbeta_l*r3temp)) -&
                  sum(r3temp1*kappadbeta_l*drnn)

    r3temp = rnnp-betas0_p

    r3temp1 = xp*exp(- pars_p(1)*r3temp)
    igamma1p = 1.0d0 / sum(r3temp1)
    dgamma1p = 0.0d0
    dgamma1p(1) = - sum(r3temp*r3temp1)
    dgamma1p(7) = sum(r3temp1*betaeta2_p)+sum(dxp*exp(- pars_p(1)*r3temp)) -&
                  sum(r3temp1*pars_p(1)*drnn)

    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
    igamma2p = 1.0d0 / sum(r3temp1)
    dgamma2p = 0.d0
    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
    dgamma2p(7) = sum(r3temp1*pars_p(6))+sum(dxp*exp(-kappadbeta_p*r3temp)) -&
                  sum(r3temp1*kappadbeta_p*drnn)


!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(nl_atoms), sigma_pp(np_atoms))
    allocate(sigma_lp(nl_atoms), sigma_pl(np_atoms))
    allocate(     s_l(nl_atoms),      s_p(np_atoms))
    allocate(rn_ltemp(nl_atoms), rn_ptemp(np_atoms))

    ! initialize accumulators
    sigma_ll = 0.0d0
    sigma_pp = 0.0d0
    dsigma_ll = 0.0d0
    dsigma_pp = 0.0d0
    sigma_pl = 0.0d0
    sigma_lp=0.0d0
    dsigma_lp_l=0.0d0
    dsigma_lp_p=0.0d0
    dsigma_pl_l=0.0d0
    dsigma_pl_p=0.0d0
    V_ll = 0.0d0
    V_pp = 0.0d0
    dV_ll = 0.0d0
    dV_pp = 0.0d0
    V_lp = 0.0d0
    V_pl = 0.0d0
    dV_lp_l = 0.0d0
    dV_lp_p = 0.0d0
    dV_pl_l = 0.0d0
    dV_pl_p = 0.0d0
    ds_l_l = 0.0d0
    ds_l_p = 0.0d0
    ds_p_l = 0.0d0
    ds_p_p = 0.0d0
    dvref_l_l=0.0d0
    dvref_l_p=0.0d0
    dvref_p_l=0.0d0
    dvref_p_p=0.0d0
    dEcoh_l = 0.0d0
    dEcoh_p = 0.0d0
    denergy = 0.0d0
    energy = 0.0d0
    rn_ltemp = 0.0d0
    rn_ptemp=0.0d0
    s_l = 0.0d0
    s_p = 0.0d0

    ! slab-slab
    do i = 1, nl_atoms
        do j = i+1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i+np_atoms) - xdata(:,j+np_atoms)   ! distance vector
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

            !----------------------------SIGMA LATTICE-----------------------------
            ! Sigma is a contribution to the neutral sphere radius.
            ! It is a list in which for each lattice atom, the contributions of the
            ! others are summed up. To enforce the cut-off, it will be later
            ! corrected by gamma.

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            rtemp1 = rtemp*(r - betas0_l)
            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1


            !-----------------------PAIR POTENTIAL LATTICE-------------------------
            ! Will later be subjected to gamma to complete the cut-off.

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

            rtemp1 = rtemp*(r - betas0_l)
            dV_ll(6) = dV_ll(6) + rtemp1

        end do
    end do

    ! projectile-projectile
    do i = 1, np_atoms
        do j = i+1, np_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j)   ! distance vector
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

            !----------------------------SIGMA PARTICLE-----------------------------
            ! Sigma is a contribution to the neutral sphere radius.
            ! It is a list in which for each particle atom, the contributions of the
            ! others are summed up. To enforce the cut-off, it will be later
            ! corrected by gamma.

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_pp(i) = sigma_pp(i) + rtemp
            sigma_pp(j) = sigma_pp(j) + rtemp

            rtemp1 = rtemp*(r - betas0_p)
            dsigma_pp(1,i) = dsigma_pp(1,i) - rtemp1
            dsigma_pp(1,j) = dsigma_pp(1,j) - rtemp1

            rtemp1 = rtemp*betaeta2_p
            dsigma_pp(7,i) = dsigma_pp(7,i) + rtemp1
            dsigma_pp(7,j) = dsigma_pp(7,j) + rtemp1

            !-----------------------PAIR POTENTIAL LATTICE-------------------------
            ! Will later be subjected to gamma to complete the cut-off.

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))   ! V_ij*gamma2*V_0
            V_pp = V_pp + rtemp

            rtemp1 = rtemp*(r - betas0_p)
            dV_pp(6) = dV_pp(6) + rtemp1
            dV_pp(7) = dV_pp(7) + rtemp*pars_p(6)

       end do
    end do

    ! projectile-slab
    do i = 1, np_atoms
        do j = 1, nl_atoms

            ! Applying PBCs
            r3temp = xdata(:,i) - xdata(:,j+np_atoms)   ! distance vector
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

            !-------------------------------MIXED SIGMA--------------------------------
            ! Contributions of both particle and lattice to neutral sphere radius
            ! To fully include the cut-off, we correct them later by gamma.
            ! Each of the mixed sigmas depends on both l and p parameters. Not all
            ! contributions need to be under the loop.
            ! sigma_lp

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_lp(j) = sigma_lp(j) + rtemp
            dsigma_lp_p(1,j) = dsigma_lp_p(1,j)-(r-betas0_p)*rtemp
            dsigma_lp_l(7,i) = dsigma_lp_l(7,i)+ betaeta2_p*rtemp

            ! sigma_pl
            rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )     ! sigma_ij*gamma1
            sigma_pl(i) = sigma_pl(i) + rtemp

            dsigma_pl_l(1,i) = dsigma_pl_l(1,i) - (r - betas0_l)*rtemp
            dsigma_pl_p(7,j) = dsigma_pl_p(7,j) + rtemp*betaeta2_l

            !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

            ! V_ij*gamma2*V_0
            rtemp = theta*exp(-kappadbeta_p*(r - betas0_p))
            V_lp = V_lp + rtemp

            dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
            dV_lp_p(7) = dV_lp_p(7) + pars_p(6)*rtemp

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l))
            V_pl = V_pl + rtemp
            dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
            dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l(6)

        end do
    end do


   ! divide by cut-off scaling factors
    sigma_ll = sigma_ll*igamma1l
        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
        dsigma_ll(7,:) = sigma_ll*pars_l(1)*beta-sigma_ll*igamma1l*dgamma1l(7)

    V_ll     =     V_ll*igamma2l*pars_l(5)
        dV_ll(5) = - V_ll/pars_l(5)
        dV_ll(6) = (dV_ll(6) * pars_l(5)/beta + V_ll*dgamma2l(6)) * igamma2l
        dV_ll(7) = V_ll*pars_l(6)-V_ll*igamma2l*dgamma2l(7)

    sigma_pp = sigma_pp*igamma1p
        dsigma_pp(1,:) = (dsigma_pp(1,:) - sigma_pp*dgamma1p(1))*igamma1p
        dsigma_pp(7,:) = (dsigma_pp(7,:) - sigma_pp*dgamma1p(7))*igamma1p

    V_pp     =     V_pp*igamma2p*pars_p(5)
        dV_pp(5) = - V_pp/pars_p(5)
        dV_pp(6) = (dV_pp(6) * pars_p(5)/beta + V_pp*dgamma2p(6)) * igamma2p
        dV_pp(7) = (pars_p(5)*dV_pp(7)-V_pp*dgamma2p(7))*igamma2p

    sigma_lp = sigma_lp*igamma1l
        ! Derivative with respect to l
        dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
        dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
        ! Derivative with respect to p
        dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
        dsigma_lp_p(7,:) = sigma_lp*pars_p(1)*beta

    sigma_pl = sigma_pl*igamma1p
        ! Derivative with respect to l
        dsigma_pl_l(1,:) = dsigma_pl_l(1,:)*igamma1p
        dsigma_pl_l(7,:) = sigma_pl*pars_l(1)*beta

        ! Derivative with respect to p
        dsigma_pl_p(1,:) = - sigma_pl*igamma1p*dgamma1p(1)
        dsigma_pl_p(7,:) = - sigma_pl*igamma1p*dgamma1p(7)

    V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
        ! Derivative with respect to l
        dV_lp_l(2) = V_lp/pars_l(2)
        dV_lp_l(5) = - V_lp/pars_l(5)
        dV_lp_l(7) = V_lp*igamma2l
        dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
        dV_lp_l(7) = V_lp*dchilp(4)*chipl - dV_lp_l(7)*dgamma2l(7)
        ! Derivative with respect to p
        dV_lp_p(2) = -V_lp/pars_p(2)
        dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l(5) / beta
        dV_lp_p(7) = V_lp*chipl*dchilp(3) + dV_lp_p(7)*chilp*pars_l(5)*igamma2l

    V_pl     =     V_pl*igamma2p*pars_p(5)*chipl
        ! Derivative with respect to l
        dV_pl_l(2) = - V_pl/pars_l(2)
        dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p(5)*chipl/beta
        dV_pl_l(7) = V_pl*chilp*dchipl(4) + dV_pl_l(7)*pars_p(5)*chipl*igamma2p
        ! Derivative with respect to p
        dV_pl_p(2) = V_pl/pars_p(2)
        dV_pl_p(5) = - V_pl/pars_p(5)
        dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
        dV_pl_p(7) = V_pl*chilp*dchipl(3) - V_pl*igamma2p*dgamma2p(7)

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------

    s_l = sigma_ll + chilp*sigma_lp
    rn_ltemp = 1.0d0/(s_l*betaeta2_l)
    s_l = -log(s_l*twelfth)/betaeta2_l

    ! Derivative with respect to l
    ds_l_l(1,:) = -s_l/pars_l(1) &
                  - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
    ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
    ds_l_l(7,:) = rn_ltemp*(dsigma_ll(7,:)+dchilp(4)*sigma_lp+chilp*dsigma_lp_l(7,:))
    ! Derivative with respect to p
    ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
    ds_l_p(2,:) = - sigma_lp*rn_ltemp*dchilp(1)
    ds_l_p(7,:) =  rn_ltemp*(dchilp(3)*sigma_lp+chilp*dsigma_lp_p(7,:))

    s_p = sigma_pp + chipl*sigma_pl
    rn_ptemp = 1.0d0 / (s_p*betaeta2_p)
    s_p = -log(s_p*twelfth)/betaeta2_p

    ! Derivative with respect to p
    ds_p_p(1,:) = -s_p/pars_p(1) &
                  - (dsigma_pp(1,:)+chipl*dsigma_pl_p(1,:))*rn_ptemp
    ds_p_p(2,:) = - sigma_pl*rn_ptemp*dchipl(1)
    ds_p_p(7,:) = rn_ptemp*(dsigma_pp(7,:)+dchipl(3)*sigma_pl+chipl*dsigma_pl_p(7,:))
    ! Derivative with respect to p
    ds_p_l(1,:) = - rn_ptemp*chipl*dsigma_pl_l(1,:)
    ds_p_l(2,:) = - sigma_pl*rn_ptemp*dchipl(2)
    ds_p_l(7,:) =   rn_ptemp*(dchipl(4)*sigma_pl+chipl*dsigma_pl_l(7,:))

!-----------------------------DENSITY CALCULATION------------------------------
! only for particle at present
    rtemp = 0.5d0/bohr2ang - betaeta2_p     ! -eta_l
    energy = pars_p(2)*exp(rtemp*s_p(5))

    denergy(1)  = -energy*(rtemp*ds_p_p(1,5)-beta*s_p(5))
    denergy(2)  = -energy*(1.d0/pars_p(2)+rtemp*ds_p_p(2,5))
    denergy(3)  = 0.0d0
    denergy(4)  = 0.0d0
    denergy(5)  = 0.0d0
    denergy(6)  = 0.0d0
    denergy(14) = energy*rtemp
    denergy(7)  = denergy(14)*ds_p_p(7,5)
    denergy(8)  = -denergy(14)*ds_p_l(1,5)
    denergy(9)  = -denergy(14)*ds_p_l(2,5)
    denergy(10) = 0.0d0
    denergy(11) = 0.0d0
    denergy(12) = 0.0d0
    denergy(13) = 0.0d0
    denergy(14) = denergy(14)*ds_p_l(7,5)

    deallocate(rn_ltemp, rn_ptemp)
    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)


end subroutine emt_ddens_fit

subroutine emt1nn(s)

!   Calculates energy and forces with nearest-neighbour EMT potential
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

    xl = nneighs * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))

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
            call pbc_dist( s%r(:,i), s%r(:,j), cell_mat, cell_imat, r)
            ! drops atoms outside (cutoff*rcut)-sphere
            if (r > cutoff*rcut) cycle

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

end subroutine emt1nn

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

subroutine ldfa(s,imass)
    !
    ! Purpose:
    !           Calculate the friction coefficient
    !
    type(atoms) :: s
    integer :: i, j
    real(8) :: fric, temp
    real(8) :: imass
    ! As implemented here, the friction coefficient is only applicable for the
    ! H-atom as calculated by Li and Wahnstrom(PRB46(1992)14528)
    ! according to Puska and Nieminen (PRB, 27, 1983, 6121), the mass still
    ! needs to be applied
    ! hbar*eta = hbar**2/mass Q(kf) (conversion between PN and LW)
    real(8), parameter :: convert= 1.00794d0*amu2mass
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
    s%dens = s%dens*convert*imass / hbar
    
    ! For simulated annealing, the Langevin dynamics are used as a heat bath
    ! But using the Au-atomic densities is too inefficient, so in this case
    ! a friction coefficient is set that is of the order of magnitude of
    ! the friction an H-atom experiences when running through Au.
    ! The friction coefficient is now 0.003 1/fs.
    if (sasteps > 0) s%dens = 0.003d0 !0.000015231d0/(imass*convert)

end subroutine ldfa

subroutine lj_e(slab, teil)

! Calculates Lennard-Jones energy

    implicit none

    type(atoms), intent(inout)    :: teil, slab

    integer :: i,j

    real(8) :: r, r6, r12, rcut
    real(8) :: eps4_l, eps4_p, eps4_pl, sigma_pl

    Epot = 0.0d0

!----------------------VALUES OF FREQUENT USE ---------------------------------

    eps4_l  = 4.0d0*pars_l(1)
    eps4_p  = 4.0d0*pars_p(1)

    eps4_pl  = sqrt(eps4_l*eps4_l)
    sigma_pl = 0.50d0*(pars_l(2) + pars_p(2))

!------------------------------ CUT-OFF ---------------------------------------

    rcut = 0.5*cell_mat(1,1)


    ! slab-slab
    do i = 1, slab%n_atoms
        do j = i+1, slab%n_atoms

            ! Applying PBCs
            call pbc_dist( slab%r(:,i), slab%r(:,j), cell_mat, cell_imat, r)
            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r6 = (pars_l(2)/r)**6
            r12 = r6*r6
            Epot = Epot + eps4_l*(r12 - r6)

        end do
    end do

    ! projectile-projectile
    do i = 1, teil%n_atoms
        do j = i+1, teil%n_atoms

            ! Applying PBCs
            call pbc_dist( teil%r(:,i), teil%r(:,j), cell_mat, cell_imat, r)
            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r6 = (pars_p(2)/r)**6
            r12 = r6*r6
            Epot = Epot + eps4_p*(r12 - r6)

        end do
    end do

    ! projectile-slab
    do i = 1, teil%n_atoms
        do j = 1, slab%n_atoms

            ! Applying PBCs
            call pbc_dist(teil%r(:,i), slab%r(:,j), cell_mat, cell_imat, r)
            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r6 = (sigma_pl/r)**6
            r12 = r6*r6
            Epot = Epot + eps4_pl*(r12 - r6)

        end do
    end do

end subroutine lj_e

subroutine lj(slab, teil)

! Calculates Lennard-Jones energy and forces

    implicit none

    type(atoms), intent(inout)    :: teil, slab

    integer :: i,j

    real(8) :: r, r6, r12, rcut, dvdr
    real(8) :: eps4_l, eps4_p, eps4_pl, sigma_pl
    real(8) :: eps4f_l, eps4f_p, eps4f_pl
    real(8), dimension(3) :: uvec

    Epot = 0.0d0
    slab%f = 0.0d0
    teil%f = 0.0d0

!----------------------VALUES OF FREQUENT USE ---------------------------------

    eps4_l  = 4.0d0*pars_l(1)
    eps4_p  = 4.0d0*pars_p(1)
    eps4f_l = 6.0d0*eps4_l/pars_l(2)
    eps4f_p = 6.0d0*eps4_p/pars_p(2)

    sigma_pl = 0.50d0*(pars_l(2) + pars_p(2))
    eps4_pl  = sqrt(eps4_l*eps4_l)
    eps4f_pl = 6.0d0*eps4_pl/sigma_pl

!------------------------------ CUT-OFF ---------------------------------------

    rcut = 0.5*cell_mat(1,1)


    ! slab-slab
    do i = 1, slab%n_atoms
        do j = i+1, slab%n_atoms

            ! Applying PBCs. Unit vector directs from i to j
            call pbc_distdir(slab%r(:,i),slab%r(:,j),cell_mat,cell_imat,r, uvec)

            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r   = pars_l(2)/r
            r6  = r**6
            r12 = r6*r6
            Epot = Epot + eps4_l*(r12 - r6)

            dvdr = eps4f_l*r*(r6 - 2.0d0*r12)
            slab%f(:,i) = slab%f(:,i) + dvdr*uvec
            slab%f(:,j) = slab%f(:,j) - dvdr*uvec

        end do
    end do

    ! projectile-projectile
    do i = 1, teil%n_atoms
        do j = i+1, teil%n_atoms

            ! Applying PBCs. Unit vector directs from i to j
            call pbc_distdir(teil%r(:,i),teil%r(:,j),cell_mat,cell_imat,r, uvec)

            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r   = pars_p(2)/r
            r6  = r**6
            r12 = r6*r6
            Epot = Epot + eps4_p*(r12 - r6)

            dvdr = eps4f_p*r*(r6 - 2.0d0*r12)
            teil%f(:,i) = teil%f(:,i) + dvdr*uvec
            teil%f(:,j) = teil%f(:,j) - dvdr*uvec

        end do
    end do

    ! projectile-slab
    do i = 1, teil%n_atoms
        do j = 1, slab%n_atoms

            ! Applying PBCs. Unit vector directs from i to j
            call pbc_distdir(teil%r(:,i),slab%r(:,j),cell_mat,cell_imat,r, uvec)

            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r   = sigma_pl/r
            r6  = r**6
            r12 = r6*r6
            Epot = Epot + eps4_pl*(r12 - r6)

            dvdr = eps4f_pl*r*(r6 - 2.0d0*r12)
            teil%f(:,i) = teil%f(:,i) + dvdr*uvec
            slab%f(:,j) = slab%f(:,j) - dvdr*uvec
        end do
    end do

end subroutine lj

subroutine lj1_e(s)

! Calculates Lennard-Jones energy

    implicit none

    type(atoms), intent(inout)    :: s

    integer :: i,j

    real(8) :: r, r6, r12, rcut
    real(8) :: eps4_l

    Epot = 0.0d0

!----------------------VALUES OF FREQUENT USE ---------------------------------

    eps4_l  = 4.0d0*pars_l(1)

!------------------------------ CUT-OFF ---------------------------------------

    rcut = 0.5*cell_mat(1,1)

    do i = 1, s%n_atoms
        do j = i+1, s%n_atoms

            ! Applying PBCs
            call pbc_dist( s%r(:,i), s%r(:,j), cell_mat, cell_imat, r)
            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r6 = (pars_l(2)/r)**6
            r12 = r6*r6
            Epot = Epot + eps4_l*(r12 - r6)

        end do
    end do

end subroutine lj1_e

subroutine lj1(s)

! Calculates Lennard-Jones energy and forces

    implicit none

    type(atoms), intent(inout)    :: s

    integer :: i,j

    real(8) :: r, r6, r12, rcut, dvdr
    real(8) :: eps4_l
    real(8) :: eps4f_l
    real(8), dimension(3) :: uvec

    Epot = 0.0d0
    s%f  = 0.0d0

!----------------------VALUES OF FREQUENT USE ---------------------------------

    eps4_l  = 4.0d0*pars_l(1)
    eps4f_l = 6.0d0*eps4_l/pars_l(2)

!------------------------------ CUT-OFF ---------------------------------------

    rcut = 0.5*cell_mat(1,1)

    do i = 1, s%n_atoms
        do j = i+1, s%n_atoms

            ! Applying PBCs. Unit vector directs from i to j
            call pbc_distdir(s%r(:,i),s%r(:,j),cell_mat,cell_imat,r, uvec)

            ! drops atoms outside cutoff
            if (r > rcut) cycle

            r   = pars_l(2)/r
            r6  = r**6
            r12 = r6*r6
            Epot = Epot + eps4_l*(r12 - r6)

            dvdr = eps4f_l*r*(r6 - 2.0d0*r12)
            s%f(:,i) = s%f(:,i) + dvdr*uvec
            s%f(:,j) = s%f(:,j) - dvdr*uvec

        end do
    end do

end subroutine lj1

end module force
