module useful_things

use atom_class

    contains



function ran1()  !returns random number between 0 - 1
 implicit none
        real(8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1 = x
end function ran1

function normal(mean,sigma) !returns a normal distribution
 implicit none
        real(8) normal,tmp
        real(8) mean,sigma   ! Sigma is the velocity we want to achieve
        integer flag
        real(8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0d0
            do while(rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do
                r1=2.0d0*ran1()-1.0d0
                r2=2.0d0*ran1()-1.0d0
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0d0*log(rsq)/rsq)
            gsave=r1*fac        ! shouldn't those two values be below zero?
            tmp=r2*fac          !
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
end function normal

subroutine lower_case(str)
    character(*), intent(in out) :: str
    integer :: i

    do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
    end do
end subroutine lower_case


subroutine norm_dist(vec1, vec2, length, norm)
    !
    ! Purpose: normalised distance between 2 vectors
    !

    integer :: length
    real(8), dimension(length) :: vec1, vec2
    real(8) :: norm, n1, n2

    norm = dot_product(vec1 - vec2, vec1 - vec2)
    n1   = dot_product(vec1,vec1)
    n2   = dot_product(vec2,vec2)

    n1 = max(n1,n2)
    if (n1 .eq. 0.0d0) then
        norm = 0.0d0
    else
        norm=Sqrt(norm/n1)
    end if

end subroutine norm_dist

function E_kin(s,mass)

    type(atoms) :: s
    real(8) :: mass, E_kin

     E_kin = 0.5d0*sum(s%v*s%v)*mass


end function E_kin

end module useful_things
