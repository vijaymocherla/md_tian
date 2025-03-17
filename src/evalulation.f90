program output_analysis
    use open_file

    implicit none
!
!   Purpose:
!           Read in and sort data from trajectory calculations
!
    real(8), parameter          :: pi       = 3.14159265359d0
    real(8), parameter          :: kB       = 8.61733238496d-5
    real(8), parameter          :: amu2mass = 103.6382d0
    real(8), parameter          :: deg2rad  = pi/180.0d0
    real(8), parameter          :: bohr2ang = 0.529177211d0

    character(len=80) :: filename, empty, folder, energy, sys, name1
    character(len=80) :: file_read_in, directory
    character(len=80) :: generala, energya, bouncesa, coorda
    character(len=80) :: generalb, energyb, bouncesb, coordb
    character(len=80) :: buffer, label
    character(len=80) :: name_p, name_l
    character(len=3)  :: pot_p, pot_l
    character(len=8)  :: str

    integer :: i,j = 0, q, k, l, b
    integer :: ios=0, itemp, line=0, pos1
    integer :: start_tr, ntrajs
    integer, dimension(2) :: wstep
    integer :: nsteps
    integer :: n_p
    integer :: time
    integer, allocatable, dimension(:) :: bounce, b_site, layer
    integer, allocatable, dimension(:) :: rl0, rl1, rl2, rl3, rl4, rl5, rl6
    integer, allocatable, dimension(:) :: rl7, rl8

    real(8) :: rtemp, a_lat, d
    real(8) :: einc, inclination, azimuth, Tsurf, step, time_at_surf
    real(8) :: mass_p, mass_l

    real(8), allocatable, dimension(:)    :: Ekin_p_s, Ekin_l_s, Epot_s, Etot_s
    real(8), allocatable, dimension(:)    :: Ekin_p, Ekin_l, Epot, Etot
    real(8), allocatable, dimension(:)    :: Eloss, theta, phi
    real(8), allocatable,dimension(:,:)   :: r_p_start, v_p_start, r_p, v_p, rmin
    real(8), allocatable,dimension(:,:,:) :: r_bounce


    ! Prompt the user to enter necessary information
    write(*,*) 'Please enter the lattice constant.'
    !read(*,*) a_lat
    a_lat=4.201
    write(*,*) 'Please enter the directory name.'
    !file_read_in = 'non_ad_pstroem_ver_6x6x4_T300_2.8_20_60'
    !name1 = 'pstroem_ver_6x6x4_T300_2.8_45_90'
    read(*,*) file_read_in
    write(*,*) 'Please enter name of unzip directory'
    !read(*,*) name1
    write(*,*) 'Please enter the number of finished trajectories'
    !read(*,*) j
    j = 1000000

!    sys = 'tar xf '//trim(file_read_in)//'.tar'
!    call system(sys)
!    sys = 'mv '//trim(name1)//' '//trim(file_read_in)
!    call system(sys)
!    sys = 'mkdir '//trim(file_read_in)//'/results'
!    call system(sys)



    ! Calculate lattice properties. Surf and a_lat need to be updated as lattice changes
    d = a_lat/sqrt(3.0d0)


    filename = trim(file_read_in)//'/md_tian.inp.1'


    directory = trim(file_read_in)//'/results/'
    generala  = trim(directory)//'scatter_all.dat'
    energya   = trim(directory)//'scatter_E.dat'
    bouncesa  = trim(directory)//'scatter_bounces.dat'
    coorda    = trim(directory)//'scatter_coord.dat'

    generalb  = trim(directory)//'noscatter_all.dat'
    energyb   = trim(directory)//'noscatter_E.dat'
    bouncesb  = trim(directory)//'noscatter_bounces.dat'
    coordb    = trim(directory)//'noscatter_coord.dat'


    call open_for_read(38,filename)


        do while (ios == 0)
        read(38, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1


        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)


            select case (label)
            case('start')
                read(buffer,*,iostat=ios) start_tr
            case('ntrajs')
                read(buffer,*,iostat=ios) ntrajs
                ntrajs = ntrajs*400
            case('Einc')
                read(buffer,*,iostat=ios) einc
            case('inclination')
                read(buffer,*,iostat=ios) inclination
                inclination = inclination*deg2rad
            case('azimuth')
                read(buffer,*,iostat=ios) azimuth
                azimuth = azimuth*deg2rad
            case('Tsurf')
                read(buffer,*,iostat=ios) Tsurf
            case('step')
                read(buffer,*,iostat=ios) step
            case('nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('wstep')
                read(buffer,*,iostat=ios) wstep
                if (wstep(1)==-1) wstep(2) = nsteps + 1
            case ('projectile')
                read(buffer, *, iostat=ios) name_p, mass_p, pot_p, itemp, &
                                            empty, empty, n_p
                mass_p=mass_p*amu2mass
            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, itemp, &
                                            empty, empty, rtemp
           case default
                print *, 'Skipping unrelated information at line', line, label
            end select
        end if
    end do ! ios
    close(38)

    allocate(r_p(3,n_p),v_p(3,n_p),rmin(3,n_p))
    allocate(r_p_start(3,n_p),v_p_start(3,n_p))
    allocate(r_bounce(5,3,n_p))
    allocate(Ekin_p_s(n_p), Ekin_l_s(n_p), Epot_s(n_p), Etot_s(n_p))
    allocate(Ekin_p(n_p), Ekin_l(n_p), Epot(n_p), Etot(n_p), Eloss(n_p), theta(n_p))
    allocate(bounce(n_p), phi(n_p), b_site(n_p), layer(n_p))
    allocate(rl0(n_p),rl1(n_p),rl2(n_p),rl3(n_p),rl4(n_p),rl5(n_p),rl6(n_p),rl7(n_p),rl8(n_p))


! How many trajectories have been calculated indeed

    if (j == 0) then
        do i= 1,ntrajs
            write(str,'(I8.8)') i
            folder   = trim(file_read_in)//'/traj/mxt_fin'//str//'.dat'
            call open_for_read(1,folder)
            read(1,'(A)',iostat=ios) buffer
            if (ios == 0) then
                j = j + 1
            end if
            close(1)
        end do
    end if

    print *, j

    ! Prepare Output files
    call open_for_append(2,generala)
    write(2,*) 'total traj = ',    j
    write(2,*) 'itraj end time (fs) No. bounces site 1st bounce &
               theta_end (°) azimuth_end (°) E_loss (eV)'
    ! site of first bounce: 0 top
    !                       1 hollow

    call open_for_append(3, energya)
    write(3,*) 'itraj Ekin_p (eV) Ekin_l (eV) Epot (eV) Etot (eV) &
               Eloss (eV)'

    call open_for_append(4, bouncesa)
    write(4,*) 'itraj time at surface (s) bounces bounce sites: x, y, z (A)'

    call open_for_append(5, coorda)
    write(5,*) 'itraj x_start (A) y_start (A) z_start (A) x_end (A)',&
            'y_end (A)', 'z_end (A)','lowest_x (A)', 'lowest_y (A)', 'lowest_z (A)'

    call open_for_append(9, generalb)
    write(9,*) 'itraj E_loss (eV) site 1st bounce layer'
    ! layer 0 : surface
    !       1 : 1st subsurface
    !       2 : 2nd subsurface etc

    call open_for_append(7,energyb)
    write(7,*) 'itraj Ekin_p (eV) Ekin_l (eV) Epot (eV) Etot (eV) &
               Eloss (eV)'

    call open_for_append(8, coordb)
    write(8,*) 'itraj x_start (A) y_start (A) z_start (A) x_end (A) &
            y_end (A) z_end (A) lowest_x (A) lowest_y (A) lowest_z (A)'

    call open_for_write(10,bouncesb)
    !write(10,*) 'itraj, bounces bounce sites: x, y, z (A)'
    write(10,*) 'itraj Eloss (eV) site 1st bounce remain on layer 1 2 3 4 5 6 7 8'


    ! Evaluate each trajectory.
    b_site = 0
    rl0 = 0
    rl1 = 0
    rl2 = 0
    rl3 = 0
    rl4 = 0
    rl5 = 0
    rl6 = 0
    rl7 = 0
    rl8 = 0
    ios = 0
    q = 0

    do i= 1,j!ntrajs
        write(str,'(I8.8)') i
        print *, i
        folder   = trim(file_read_in)//'/traj/mxt_fin'//str//'.dat'
        call open_for_read(1,folder)
        read(1,'(A15,f15.5)',iostat=ios) empty, Ekin_p_s
        if(ios .ne. 0 ) print *, i
        if (ios == 0) then
            q = q + 1
            !print *, 'Attending to traj', q
            read(1,'(A15,f15.5)',iostat=ios) empty, Ekin_l_s
            read(1,'(A15,f15.5)',iostat=ios) empty, Epot_s
            ! These lines are necessary because we included the reference energy in later calculations
            if (mod(i+1,2500)==0) read(1,*) empty, empty
            if (mod(i,2500)==0) read(1,*) empty, empty
            read(1,'(A15,f15.5)',iostat=ios) empty, Etot_s
            read(1,*)
            read(1,*) r_p_start
            read(1,*)
            read(1,*) v_p_start
        end if
            read(1,*,iostat=ios)
        if (ios == 0) then
            read(1,'(A15,I8)',iostat=ios) empty, time
            read(1,'(A15,f15.5)',iostat=ios) empty, Ekin_p
            read(1,'(A15,f15.5)',iostat=ios) empty, Ekin_l
            read(1,'(A15,f15.5)',iostat=ios) empty, Epot
            ! These lines are necessary because we included the reference energy in later calculations
            if (mod(i+1,2500)==0) read(1,*) empty, empty
            if (mod(i,2500)==0) read(1,*) empty, empty
            read(1,'(A15,f15.5)',iostat=ios) empty, Etot
            read(1,*)
            read(1,*) r_p
            read(1,*)
            read(1,*) v_p
            read(1,*)
            read(1,*) rmin
            read(1,*)
            read(1,*) time_at_surf
            read(1,*)
            read(1,*) bounce
            read(1,*)
            do k = 1, n_p
                b = bounce(k)
                if (b > 5) b = 5
                read(1,*) r_bounce(1:b,1,k), r_bounce(1:b,2,k), r_bounce(1:b,3,k)
            end do


            Eloss = Ekin_p_s - Ekin_p
            ! site at which first bounce happened. Either hollow (2) or top (1). But still needs deciding how done.
            ! Perhaps via height? If bounced below certain height, hollow?
            b_site = 1


            do k = 1, n_p
                !phi(k)   = acos(v_p(1,k)/sqrt(v_p(1,k)**2+v_p(2,k)**2))/deg2rad
                phi(k) = atan2(v_p(2,k),v_p(1,k))/deg2rad !Needs to be done this way, otherwise no neg angles.
                theta(k) = acos(v_p(3,k)/sqrt(v_p(1,k)**2+v_p(2,k)**2+v_p(3,k)**2))/deg2rad



                if (r_p(3,k) < 7.0d0 .and. r_p(3,k) > 6.0d0) then
                    ! Write into different files
                    write(2,'(4I8, 3f20.5)') q, time, bounce(1), b_site(1), theta(1), phi(1), Eloss(1)
                    write(3,'(I8, 5f20.5)') q, Ekin_p(1), Ekin_l(1), Epot(1), Etot(1), Eloss(1)
                    write(4,'(I8, f20.5, I8, 100f20.5)') q, time_at_surf, bounce, r_bounce(:,:,1)
                    write(5,'(I8, 9f20.5)') q, r_p_start(:,1), r_p(:,1), rmin(:,1)


                    if (n_p > 1) then
                        write(2,'(2A8,2I8,3f20.5)') '', '', bounce(k), b_site(k), theta(k), phi(k), Eloss(k)
                        write(4,'(3A8,100f20.5)') '', r_bounce(:,:,k)
                        write(3,'(A8,5f20.5)') '', Ekin_p(k), Ekin_l(k), Epot(k), Etot(k), Eloss(k)
                        write(5,'(A8,9f20.5)') '', r_p_start(:,k), r_p(:,k), rmin(:,k)
                    end if

                else


                    ! Which layer does the particle end up at?
                    ! How many atoms remain per layer?
                    if (r_p(3,k) > 0.0d0 .and. r_p(3,k) < 3.0d0 ) then
                        layer(k) = 0
                        rl0(k) = rl0(k) + 1
                    end if
                    if (r_p(3,k) < 0.0d0 .and. r_p(3,k) > -d) then
                        layer(k) = 1
                        rl1(k) = rl1(k) + 1
                    end if
                    if (r_p(3,k) < -d .and. r_p(3,k) > -2*d) then
                        layer(k) = 2
                        rl2(k) = rl2(k) + 1
                    end if
                    if (r_p(3,k) < -2*d .and. r_p(3,k) > -3*d) then
                        layer(k) = 3
                        rl3(k) = rl3(k) + 1
                    end if
                    if (r_p(3,k) < -3*d .and. r_p(3,k) > -4*d) then
                        layer(k) = 4
                        rl4(k) = rl4(k) + 1
                    end if
                    if (r_p(3,k) < -4*d .and. r_p(3,k) > -5*d) then
                        layer(k) = 5
                        rl5(k) = rl5(k) + 1
                    end if
                    if (r_p(3,k) < -5*d .and. r_p(3,k) > -6*d) then
                        layer(k) = 6
                        rl6(k) = rl6(k) + 1
                    end if
                    if (r_p(3,k) < -6*d .and. r_p(3,k) > -7*d) then
                        layer(k) = 7
                        rl7(k) = rl7(k) + 1
                    end if
                    if (r_p(3,k) < -7*d .and. r_p(3,k) > -8*d) then
                        layer(k) = 8
                        rl8(k) = rl8(k) + 1
                    end if

                    write(9,'(2I8,f20.5, I8)') q, b_site(1), Eloss(1), layer(1)
                    write(7,'(I8, 5f20.5)') q, Ekin_p(1), Ekin_l(1), Epot(1), Etot(1), Eloss(1)
                    write(8,'(I8, 9f20.5)') q, r_p_start(:,1), r_p(:,1), rmin(:,1)
                    !write(10,'(I8, f15.5, I8, 100f15.5)') q, bounce, r_bounce(:,:,1)

                    if (n_p > 1) then
                            write(9,'(A8,I8,f20.5,I8)') '', b_site(k), Eloss(k), layer(k)
                            write(7,'(A8,5f20.5)') '', Ekin_p(k), Ekin_l(k), Epot(k), Etot(k), Eloss(k)
                            write(8,'(A8,9f20.5)') '', r_p_start(:,k), r_p(:,k), rmin(:,k)
                    !       write(10,'(3A8,100f15.5)') '', r_bounce(:,:,k)
                    end if

                end if

            end do


        end if
        close(1)
    end do
    ! write how many atoms remain in each layer
    write(10,*) rl0, rl1, rl2, rl3, rl4, rl5, rl6, rl7, rl8


    close(2)
    close(3)
    close(4)
    close(5)
    close(7)
    close(8)
    close(9)
    close(10)
    deallocate(r_p,v_p,rmin)
    deallocate(r_p_start,v_p_start)
    deallocate(Ekin_p_s, Ekin_l_s, Epot_s, Etot_s, bounce, b_site)
    deallocate(Ekin_p, Ekin_l, Epot, Etot, Eloss, theta, phi)
    deallocate(rl0,rl1,rl2,rl3,rl4,rl5,rl6,rl7,rl8)
end program output_analysis

