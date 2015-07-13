module eloss_arrays
    !
    ! Purpose:
    !           Contains eloss binning.
use open_file
implicit none
save

real(8), parameter          :: pi       = 3.14159265359d0
real(8), parameter          :: kB       = 8.61733238496d-5
real(8), parameter          :: amu2mass = 103.6382d0
real(8), parameter          :: deg2rad  = pi/180.0d0
real(8), parameter          :: bohr2ang = 0.529177211d0
real(8), parameter          :: DE       = 0.01d0 ! size of energy binning interval (eV)
real(8), parameter          :: Dang     = 5.d0 ! size of angle-selection-interval
                                            ! CAREFUL: if you change these, you also need to change the output
real(8), parameter          :: Dangd    = 2.d0 ! size of angle-selection-interval for angle-distribution
real(8), parameter          :: up       = 79.5d0 ! upper timelimit (µs)
real(8), parameter          :: distance = 2.5d0 ! distance to detector (dm)
real(8), parameter          :: Dt       = 0.1d0 ! time binning step (µs)

contains

subroutine beineloss(Eloss,rp, nbounce, rmin, nbin, ival, array)
!
! Purpose: Energy loss binning
!
real(8), intent(in)                     :: Eloss
real(8), dimension(3), intent(in)       :: rp, rmin
integer, intent(in)                     :: nbounce, nbin
real(8), dimension(:), intent(in)       :: ival  ! array with binning intervals
integer, dimension(:,:), intent(inout)  :: array ! array to bin Eloss

integer :: j, k

        do j=1,nbin
        if (Eloss > ival(j) .and. Eloss < ival(j+1) .and. &
            rp(3) < 7.0d0 .and. rp(3) > 5.0d0) then
            ! all
            array(1,j)=array(1,j)+1
            ! Single
            if (nbounce == 1) array(2,j)=array(2,j)+1
            ! Double
            if (nbounce == 2) array(3,j)=array(3,j)+1
            ! Triple
            if (nbounce == 3) array(4,j)=array(4,j)+1
            ! Multiple
            if (nbounce  > 3) array(5,j)=array(5,j)+1
            ! Non-Penetrating
            if (rmin(3) >= 0) array(6,j)=array(6,j)+1
            ! Penetrating
            if (rmin(3)  < 0) array(7,j)=array(7,j)+1
        end if
        end do

end subroutine beineloss

subroutine bintof(vp,rp, nbounce, rmin, tbin, tval, tarray)
!
! Purpose: tof binning
!          Here, I first calculate (distance/vp) = t and then bin with regards to the
!          timestep Dt. So, we get I(t)dt here
!
!
real(8), dimension(3), intent(in)       :: vp, rp, rmin
integer, intent(in)                     :: nbounce, tbin
real(8), dimension(:), intent(in)       :: tval  ! array with binning intervals
integer, dimension(:,:), intent(inout)  :: tarray ! array to bin Eloss
real(8)                                 :: time ! distance divided by velocity = time

integer :: j, k

time= distance/Sqrt(vp(1)**2+vp(2)**2+vp(3)**2)
! Bin according to intervals
! Bin all
do j=1,tbin
    if (time > tval(j) .and. time < tval(j+1) .and. &
        rp(3) < 7.0d0 .and. rp(3) > 5.0d0) then
        ! all
        tarray(1,j)=tarray(1,j)+1
        ! single bounce
        if (nbounce == 1) tarray(2,j)=tarray(2,j)+1
        ! double bounce
        if (nbounce == 2) tarray(3,j)=tarray(3,j)+1
        ! triple bounce
        if (nbounce == 3) tarray(4,j)=tarray(4,j)+1
        ! multiple bounce
        if (nbounce  > 3) tarray(5,j)=tarray(5,j)+1
        ! non-penetrating
        if (rmin(3) >= 0) tarray(6,j)=tarray(6,j)+1
        ! Penetrating
        if (rmin(3)  < 0) tarray(7,j)=tarray(7,j)+1
    end if
end do


end subroutine bintof
end module eloss_arrays
