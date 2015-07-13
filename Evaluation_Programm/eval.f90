program trajeval
    use open_file
    use eloss_arrays
    implicit none

! Purpose:
!       This program reads in from trajectories and creates files for
!   1. AIELD (Angle Integrated Energy Loss Distribution), a dat-file (done)
!   2. An table that shows where atoms have ended up after the end of traj (done)
!   3. Angular distributions (done)
!   4. What scattering results look like under different angles (for comparison to experiment) (done)
!   5. All for both tof and Eloss. (done)
!
! Author:   SM Janke
! Authored: 2015/7
!

! The programm will furthermore ask before it overwrites files.
character(len=100)          :: pos_init_file, buffer, label, filename,f1
character(len=80)           :: empty
character(len=8)            :: str

! read in from the input file or related to it
integer                     :: ios =0, pos1, line
integer                     :: i, j, k
integer                     :: l, x ! integers to check if calculations completed
integer                     :: ntraj, ang
integer, dimension(3)       :: cell
real(8)                     :: Einc
real(8)                     :: alat, d0
real(8)                     :: ttraj = 1000
real(8), allocatable,dimension(:) :: theta, phi

! related to read-in from mxt_fin-files
real(8)                     :: ek_ps, ek_ss, epot_s, epef_s, etot_s
real(8)                     :: Ekin_p, Ekin_l, Epot, pEF, Etot
real(8)                     :: Eref, Eloss
real(8), dimension(3)       :: rp_s, rp, vp_s, vp, rmin
real(8), dimension(5,3)     :: r_bounce = 0.0d0
integer                     :: time, tatsurf, nbounce, b

! related to evaluation
integer                                 :: nbin, nlay, tbin, nth, nph
real(8), allocatable, dimension(:,:)     :: angar ! for selection of angles
integer, allocatable, dimension(:,:)    :: eaild ! array to bin Eloss for AIELD into
integer, allocatable, dimension(:,:)    :: taild,dtaild ! array to bin tof for
                                                        ! AITOFD into and density conversion
integer, allocatable, dimension(:,:,:)  :: angelos ! array to bin Eloss after angles
integer, allocatable, dimension(:,:,:)  :: angetof ! array to bin tof after angles
integer, allocatable, dimension(:,:,:)  :: dangetof ! array tof for density conversion
real(8), allocatable, dimension(:)      :: ival, tval,thar,phar  ! array with binning intervals
real(8), allocatable, dimension(:)      :: iscat ! array with numbers of possible
                                             ! whereabout layers
integer, allocatable, dimension(:)      :: elay  ! array with whereabout according
                                             ! to layers at end of trajectory
integer, allocatable, dimension(:)      :: zlay  ! array with lowest z
integer, allocatable, dimension(:,:)    :: wlay  ! where did they come from?
real(8)                                 :: ni,tot
real(8)                                 :: phio, thetao ! exit-angles
integer, allocatable, dimension(:,:,:)   :: curu ! First step to circular angluar dist.
real(8)                                 :: time1 ! distance divided by velocity = time

! output related
character(len=4)                        :: thout, phiout
character(len=1)                        :: over
character(len=200)                      :: sys, outname, nout
logical                                 :: exists,exists1

!dummies for debugging
integer :: bb=0
logical :: used=.true.

! Read in from input file
if (iargc() == 0) stop " I need an input file"
call getarg(1, pos_init_file)

call open_for_read(38, pos_init_file)
do while (ios == 0)
    read(38, '(A)', iostat=ios) buffer
    if (ios == 0) then
        line = line + 1
    ! Find the first instance of whitespace.  Split label and data.
        pos1 = scan(buffer, ' ')
        label = buffer(1:pos1)
        buffer = buffer(pos1+1:)

        select case (label)
        case('name')
            read(buffer,*,iostat=ios) filename
        case('outname')
            read(buffer,*,iostat=ios) outname
        case('celldim')
            read(buffer,*,iostat=ios) cell
        case('ntraj')
            read(buffer,*,iostat=ios) ntraj
        case('Einc')
            read(buffer,*,iostat=ios) Einc
        case('ttraj')
            read(buffer,*,iostat=ios) ttraj
        case('alat')
            read(buffer,*,iostat=ios) alat
        case('azimuth')
            read(buffer,*,iostat=ios) ang
            allocate(phi(ang+1)) ! needs one more because first
                                    ! element is going to be number of angles
            read(buffer,*,iostat=ios) phi
        case('inclination')
            read(buffer,*,iostat=ios) ang
            allocate(theta(ang+1)) ! needs one more because first
                                    ! element is going to be number of angles
            read(buffer,*,iostat=ios) theta
        case default
            print *, 'Skipping invalid label at line', line, label
        end select
    end if
end do ! ios
close(38)

!------------------------------------------------------------------------------
!
!                           PREPARE CALCULATION'S OUTPUT
!
!------------------------------------------------------------------------------
! See if output directory already exists
nout='output/'//trim(outname)
inquire(directory=trim(nout),exist=exists)
if (.not. exists) then
    call system('mkdir '//trim(nout))
    else
        inquire(file=trim(nout)//'/AIELD.dat',exist=exists1)
        if (exists1) then
            write(*,*) "File already exists. If you continue, the old files &
                        will be replace."
            write(*,*) "Do you want to proceed? (y/n)"
            read(*,*) over
            if (over == 'n' .or. 'N') stop "File won't be overwritten."
        end if
end if

sys = 'tar xfz '//trim(filename)//&
    '.tar.gz -C ../output_analysis'
inquire(file=trim(filename)//'.tar.gz',exist=exists)
if (exists) then
    print *, 'extracting file...'
    call system(sys)
    print *, 'extraction ended.'
end if


! Prepare array for selection of angles. For comparison to experiment
! This array is structured such:
! angar(a,b)
! a: number of angular conditions
! b1: thmin, b2: thmax, b3: phimin, b4: phimax
allocate(angar(ang,4))
do i=2,ang+1
    angar(i-1,1) = theta(i)-Dang
    angar(i-1,2) = theta(i)+Dang
    angar(i-1,3) = phi(i)-Dang
    angar(i-1,4) = phi(i)+Dang
    ! What if back-scattering is monitored?
    if (angar(i-1,1) < -Dang .and. angar(i-1,2) < -Dang) then
        angar(i-1,1) = Abs(angar(i-1,1))
        angar(i-1,2) = Abs(angar(i-1,2))
        angar(i-1,3) = -180+angar(i-1,3)
        angar(i-1,4) = -180+angar(i-1,4)
    end if
end do

! Arrays for angular distribution
nth=90.d0/Dangd
nph=360.d0/Dangd
allocate(thar(nth),phar(nph),curu(7,nth,nph))
thar = 0.0d0
phar=0.0d0
curu =0
ni=-Dangd
do i=1,nth
    ni = ni+Dangd
    thar(i)=ni
end do
ni=-180-Dangd
do i=1,nph
    ni = ni+Dangd
    phar(i)=ni
end do

! Loop over trajectories
do i=1, ntraj
    if (i == 50000) print *, i
    if (mod(i,50000)==0) print *, i
    write(str,'(I8.8)') i
    f1 = trim(filename)//'/traj/mxt_fin'//str//'.dat'
    call open_for_read(54,f1)
    read(54,'(A15,f15.5)', iostat=ios) empty, ek_ps
    ! Check if the file exists
    if (ios .ne. 0) then
        bb=bb+1
        cycle !cycle if file not there
    end if
    ! Check if the file contains all information
    l=0
    x=0
    do
        read(54,*,iostat=ios) empty
        l = l +ios
        if (l<-1) exit
        x=x+1
    end do
    rewind(54)
    read(54,'(A15,f15.5)',iostat=ios) empty, ek_ps
    if (x < 28) then
        print *,i
        bb = bb+1
        cycle   ! cycle if file too short
    end if
    read(54,'(A15,f15.5)') empty, ek_ss
    read(54,'(A15,f15.5)') empty, epot_s
    read(54,'(A15,f15.5)') empty, epef_s
    read(54,'(A15,f15.5)') empty, Eref
    read(54,'(A15,f15.5)') empty, etot_s
    read(54,*)
    read(54,*) rp_s(1), rp_s(2), rp_s(3)
    read(54,*)
    read(54,*) vp_s(1), vp_s(2), vp_s(3)
    read(54,*)
    read(54,'(A15,I8)'), empty, time
    read(54,'(A15,f15.5)') empty, Ekin_p
    read(54,'(A15,f15.5)') empty, Ekin_l
    read(54,'(A15,f15.5)') empty, Epot
    read(54,'(A15,f15.5)') empty, pEF
    read(54,*)
    read(54,'(A15,f15.5)') empty, Etot
    read(54,*)
    read(54,*) rp(1), rp(2), rp(3)
    read(54,*)
    read(54,*) vp(1), vp(2), vp(3)
    read(54,*)
    read(54,*) rmin(1), rmin(2), rmin(3)
    read(54,*)
    read(54,*) tatsurf
    read(54,*)
    read(54,*) nbounce
    read(54,*)
    b = nbounce
    if (nbounce > 5) b = 5
    read(54,*) r_bounce(1:b,1), r_bounce(1:b,2), r_bounce(1:b,3)
    close(54)

    ! Begin evalulation
    Eloss = ek_ps-Ekin_p

    !--------------------------------------------------------------------------
    !
    !                                ELOSS and TOF
    !
    !--------------------------------------------------------------------------
    ! calculate the number of binning intervals.
    nbin=(Einc+0.5d0)/DE
    tbin=up/Dt
    ! build array that collects energy bins & array with binning intervals
    if(.not. allocated(ival)) then
        allocate(eaild(7,nbin),ival(nbin),angelos(ang,7,nbin))
        allocate(taild(7,tbin),tval(tbin),angetof(ang,7,tbin))
        allocate(dtaild(7,tbin),dangetof(ang,7,tbin))
        ival =0.0d0
        tval =0.0d0
        eaild=0.0d0
        taild=0.0d0
        dtaild=0.0d0
        angelos = 0
        angetof = 0
        dangetof = 0
        ! Build arrays with binning intervals
        ni = 0.d0
        do j=1,nbin
            ival(j) = -0.5d0+ni
            ni = ni+DE
        end do
        ni = 0.0d0
        do j=1,tbin
            tval(j) = ni
            ni = ni+Dt
        end do
    end if
    ! build array for binning intervals for tof
    if (.not. allocated(tval)) then
        allocate(tval(tbin))
    end if

    ! Angle Integrated Energy Loss Distributions

    call beineloss(Eloss, rp, nbounce, rmin, nbin, ival, eaild)
    ! Angle Integrated tof distribution, Flux

    call bintof   (vp   , rp, nbounce, rmin, tbin, tval, taild)

    ! Energy loss with regards to angle
    ! First, calculate the exit-angles for all scattered
    if (rp(3) < 7.0d0 .and. rp(3) > 4.0d0) then
        !Needs to be done this way, otherwise no neg angles for phi
        phio = atan2(vp(2),vp(1))/deg2rad
        thetao = acos(vp(3)/sqrt(vp(1)**2+vp(2)**2+vp(3)**2))/deg2rad

        ! Select now according to selection-angle

        do k = 1, ang
            if (angar(k,1) > 0.0d0 .and. angar(k,2) > 0.0d0 ) then
            if (thetao >= angar(k,1) .and. thetao <= angar(k,2)&
                .and. phio >= angar(k,3) .and. phio <= angar(k,4)) then
                ! Eloss
                call beineloss(Eloss,rp, nbounce, rmin, nbin, ival,&
                               angelos(k,:,:))
                ! TOF

                call bintof   (vp   , rp, nbounce, rmin, tbin, tval,&
                               angetof(k,:,:))
            end if
            ! What if normal incidence conditions or something close to it?
            else
            if (thetao >= angar(k,1) .and. thetao <= angar(k,2)&
            .and. phio >= angar(k,3) .and. phio <= angar(k,4)) then
                ! Eloss
                call beineloss(Eloss,rp, nbounce, rmin, nbin, ival,&
                               angelos(k,:,:))
                ! TOF
                call bintof   (vp   , rp, nbounce, rmin, tbin, tval,&
                               angetof(k,:,:))
            elseif (thetao >= angar(k,1) .and. thetao <= angar(k,2)&
            .and. phio >= angar(k,3)-180 .and. phio <= angar(k,4)-180) then
                ! Eloss
                call beineloss(Eloss,rp, nbounce, rmin, nbin, ival,&
                               angelos(k,:,:))
                ! TOF
                call bintof   (vp   , rp, nbounce, rmin, tbin, tval,&
                               angetof(k,:,:))
            end if
            end if
        end do
    end if

    !--------------------------------------------------------------------------
    !
    !                           WHERE DO ATOMS END UP?
    !
    !--------------------------------------------------------------------------
    ! z-position of H-atom and verdict at end of trajectory
    ! lowest(A)     highest(A)         verdict
    !   3.0             6.1             scattered
    !   0.0             < 3.0           on surface
    !   -a0/sqrt(3)     < 0.0           1st subsurface
    !   -2 a0/sqrt(3)   < a0/sqrt(3)    2nd subsurface
    !   -3 a0/sqrt(3)   < 2 a0/sqrt(3)  3rd subsurface
    ! etc until thru (which is < z lowest level - a0/sqrt(3))
    nlay = cell(3) + 2
    d0 = -alat/sqrt(3.)
    if (allocated(iscat)) then
    else
        allocate(iscat(nlay),elay(nlay),zlay(nlay),wlay(nlay,nlay))
        elay=0.0d0
        zlay=0.0d0
        wlay=0.0d0
        iscat(1) = 4.0d0
        iscat(2) = 0.0d0
        do j= 3, nlay
            iscat(j) = (j-2)*d0
        end do
    end if

    ! Where is the H-atom at the end of the trajectory with respect to its
    ! z-coordinate
    ! And where did they come from?
    ! wlay(a,b)
    ! a: where are they now
    ! b: where did they have their lowest

    ! scattered
    !-----------
    used=.true.
    if(rp(3) .gt. iscat(1)) then   ! did they scatter?
        elay(1) = elay(1) +1
        used=.false.

        ! Where did they come from?
        ! did they scatter above surface? This should not give a hit.
        ! If it does, something wrong.
        if(rmin(3) .gt. iscat(1)) then
            wlay(1,1) = wlay(1,1) +1
        ! did they pass through? This should not give a hit.
        ! If it does, something wrong.
        elseif (rmin(3) .lt. iscat(nlay)) then
            wlay(1,nlay) = wlay(1,nlay) +1
        end if
        ! Where did they scatter in slab?
        do j=2,nlay-1
            if (rmin(3) .gt. iscat(j) .and. rmin(3) < iscat(j-1)) then
                wlay(1,j) = wlay(1,j) +1
            end if
        end do

    ! Did they pass through?
    ! Questions about lowest z are irrelavant here, since, if they go through,
    ! their lowest z is obviously outside of the slab.
    !-----------------------

    elseif (rp(3) < iscat(nlay)) then
        elay(nlay) = elay(nlay) +1
        wlay(nlay,nlay) = wlay(nlay,nlay) +1
        used=.false.
    end if

    ! Remaining in surface
    ! --------------------
    do j=2,nlay
        ! Are they between layer j-1 and j?
        if (rp(3) > iscat(j) .and. rp(3) < iscat(j-1)) then
            used=.false.
            elay(j) = elay(j) +1
            ! Where was their lowest z?
            ! On surface?
            if(rmin(3) > iscat(2)) then
                wlay(j,1) = wlay(j,1) +1
            ! Through?
            elseif (rmin(3) < iscat(nlay)) then
                wlay(j,nlay) = wlay(j,nlay) +1
            end if
            ! or somewhere in slab?
            do k=3,nlay
                if (rmin(3) > iscat(k) .and. rmin(3) < iscat(k-1)) then
                    wlay(j,k-1) = wlay(j,k-1) +1
                end if
            end do
        end if
    end do

    if (used==.true.) then
        print *, rp

    endif

    !--------------------------------------------------------------------------
    !
    !                          ANGULAR DISTRIBUTIONS
    !
    !--------------------------------------------------------------------------
    ! Circular Angular Distribution
    !
    ! Bin all angles into the intervals into which they belong
    ! This time, I need a matrix (c, a,b)
    ! a: how many counts in th from 0 to 90 in Dang
    ! b: how many counts in phi from -180 to 180 in Dang
    ! c: 1 all, 2 single, 3 double, 4 triple, 5 multi 6 non-penetrating
    !    7 penetrating

    ! First, choose only such trajectories for this that actually leave the slab
    ! then, calculate their exit angles.
    if (rp(3) < 7.0d0 .and. rp(3) > iscat(1)) then
        !Needs to be done this way, otherwise no neg angles for phi
        phio = atan2(vp(2),vp(1))/deg2rad
        thetao = acos(vp(3)/sqrt(vp(1)**2+vp(2)**2+vp(3)**2))/deg2rad

    do j=1,nth
    do k=1,nph
        if(thetao< thar(j)+Dangd .and. thetao>thar(j) .and. &
           phio< phar(k)+Dangd .and. phio> phar(k)) then
           !all
           curu(1,j,k)=curu(1,j,k)+1
           if (nbounce == 1) curu(2,j,k)=curu(2,j,k)+1
           if (nbounce == 2) curu(3,j,k)=curu(3,j,k)+1
           if (nbounce == 3) curu(4,j,k)=curu(4,j,k)+1
           if (nbounce  > 3) curu(5,j,k)=curu(5,j,k)+1
           if (rmin(3) >= 0) curu(6,j,k)=curu(6,j,k)+1
           if (rmin(3)  < 0) curu(7,j,k)=curu(7,j,k)+1
        end if
    end do
    end do
    end if

end do

!------------------------------------------------------------------------------
!
!                         TIME OF FLIGHT, FLUX TO DENSITY
!
!------------------------------------------------------------------------------
! For TOF we have everything in flux. Here, we convert flux into density
! We can do that by dividing by the velocity
! I(t) dt Dt / distance = I(t) dt / v

do j=1,7
    do i=1,tbin
        dtaild(j,i) = (tval(i)+Dt/2.d0)*taild(j,i)/distance
    end do
end do

do k=1,ang
    do j=1,7
        do i=1,tbin
            dangetof(k,j,i) = (tval(i)+Dt/2.d0)*angetof(k,j,i)/distance
        end do
    end do
end do


!------------------------------------------------------------------------------
!
!                                   OUTPUT
!
!------------------------------------------------------------------------------

! LASTLY:
! 2. untar input directories
! 6. Tar input directory again.
! 7. Check if we get similar output to mathematica.

! write AIELD
!-------------
call open_for_write(77,trim(nout)//'/AIELD.dat')
write(77,'(7A8)') 'Eloss/eV','all','single','double', 'triple', 'multiple', &
                   'non-pene', 'pene'
do i=1,nbin
    write(77,'(f8.5,7I6)') (ival(i)+DE/2.d0),eaild(1,i),eaild(2,i),eaild(3,i),&
                           eaild(4,i), eaild(5,i), eaild(6,i), eaild(7,i)
end do
 close(77)


! angelos(a,b,c) a: angle, b: condition, c: nbin
!write Angle-resolved Eloss spectra
!-----------------------------------
do i=1,ang
    write(thout,'(f4.1)') theta(i+1)
    write(phiout,'(f4.1)') phi(i+1)
    call open_for_write(90,trim(nout)//'/Eloss_spectra_tho'//trim(thout)//&
         '_phio'//trim(phiout)//'.dat')
    write(90,'(7A8)') 'Eloss/eV','all','single','double', 'triple', 'multiple', &
                      'non-pene', 'pene'
    do j=1, nbin
        write(90,'(f8.5,7I20)')(ival(j)+DE/2.d0),angelos(i,1,j),angelos(i,2,j),&
                               angelos(i,3,j), angelos(i,4,j), angelos(i,5,j),&
                               angelos(i,6,j), angelos(i,7,j)
    end do
end do
close(90)

! write TOF
!----------
call open_for_write(77,trim(nout)//'/AITOFD.dat')
write(77,'(7A8)') 'tof/µs','all','single','double', 'triple', 'multiple', &
                   'non-pene', 'pene'
do i=1,tbin
    write(77,'(f8.5,7I6)')(tval(i)+Dt/2.d0),dtaild(1,i),dtaild(2,i),dtaild(3,i),&
                           dtaild(4,i), dtaild(5,i), dtaild(6,i), dtaild(7,i)
end do
 close(77)
! angelos(a,b,c) a: angle, b: condition, c: nbin
!write Angle-resolved TOF spectra
!-----------------------------------
do i=1,ang
    write(thout,'(f4.1)') theta(i+1)
    write(phiout,'(f4.1)') phi(i+1)
    call open_for_write(90,trim(nout)//'/TOF_spectra_tho'//trim(thout)//&
         '_phio'//trim(phiout)//'.dat')
    write(90,'(7A8)') 'TOF/µs','all','single','double', 'triple', 'multiple',&
                      'non-pene', 'pene'
    do j=1, tbin
        write(90,'(f8.5,7I20)')(tval(j)+Dt/2.d0),dangetof(i,1,j),&
                                dangetof(i,2,j),dangetof(i,3,j), &
                                dangetof(i,4,j), dangetof(i,5,j),&
                                dangetof(i,6,j), dangetof(i,7,j)
    end do
end do
close(90)

! write statistics.
!------------------
! Also write statistics of where they scattered.
tot=sum(elay)
call open_for_write(89,trim(nout)//'/statistics.dat')
write(89,'(A20,1f10.1)') 'Total number of trajectories:', tot
write(89,*) 'scattered:', elay(1)/tot
write(89,*) 'remaining:', sum(elay(2:nlay))/tot
write(89,*) 'In which layer did the trajectories scatter?'
write(89,*) 'scattered, surface, in', nlay-3, 'subsurface layers, through'
write(89,'(20I8)') elay
write(89,*) 'origin, surface, in', nlay-3, 'subsurface layers, through'
write(89,'(A20, 20I8)') 'surface ', wlay(2,:)
do i = 3, nlay-1
    write(89,'(I2,A18, 20I8)') i-2,' subsurface layer ', wlay(i,:)
end do
write(89,'(A20, 20I8)') 'through ', wlay(nlay,:)
close(89)

!write Angular Distributions
!----------------------------
call open_for_write(57,trim(nout)//'/circ_ang_dist.dat')
write(57,*) 'Angualar distributions'
write(57,*) 'Theta_out are the lines and phi_out the columns in counts.'
write(57,*) "For plots, don't forget to convert angles to radiants."
write(57,*) 'All bounce events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(1,:,i)
end do
write(57,*) ' '
write(57,*) 'Single bounce events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(2,:,i)
end do
write(57,*) ' '
write(57,*) 'Double bounce events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(3,:,i)
end do
write(57,*) ' '
write(57,*) 'Triple bounce events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(4,:,i)
end do
write(57,*) ' '
write(57,*) 'Multiple bounce events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(5,:,i)
end do
write(57,*) ' '
write(57,*) 'penetrating events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(6,:,i)
end do
write(57,*) ' '
write(57,*) 'non-penetrating events'
write(57,'(20f8.2)') 0.0, thar+Dangd/2.d0
do i=1,nph
    write(57,'(f8.2,50I8)') phar(i)+Dangd/2.d0, curu(7,:,i)
end do
write(57,*) ' '
close(57)


write(*,*) 'All files produced. Cleaning up now.'
sys = 'rm -r '//trim(filename)
!'tar cfz '//trim(filename)//'.tar.gz '//trim(filename)
print *, sys
call system(sys)
print *, 'The number of files that were damaged or not present is:', bb
deallocate(ival,eaild,angelos,elay,wlay,angar,theta, phi)
deallocate(tval,taild,dtaild,dangetof,thar,phar)
!deallocate(tval,taild,dtaild,angetof,dangetof,thar,phar) !somehow, this one causes problems with the program
end program trajeval

