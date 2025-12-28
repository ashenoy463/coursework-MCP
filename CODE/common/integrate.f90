module INTEGRATE

    use iso_fortran_env, only: IK => int32, RK => real64

    implicit none
    public::TRAPEZOID
    public::SIMPSON

contains

    subroutine TRAPEZOID(F,A,B,N,INTEGRAL)
        ! Composite Trapezoidal Rule
        !
        ! Integrates function F on [A,B] with N subintervals

        implicit none

        real, external      :: F            ! Integrand
        real, intent(in)    :: A            ! Interval start
        real, intent(in)    :: B            ! Interval end
        integer, intent(in) :: N            ! No. of subintervals
        real, intent(out)   :: INTEGRAL     ! Result

        integer :: I 
        real    :: H

        H = (B-A)/N

        INTEGRAL = 0.5*(F(A) + F(B))

        do I = 1,N-1
            INTEGRAL = INTEGRAL + F(A+H*I)
        end do 
        
        INTEGRAL = H*INTEGRAL
        
    end subroutine TRAPEZOID


    subroutine SIMPSON(F,A,B,N,INTEGRAL)
        ! Composite Simpson's 1/3rd Rule
        !
        ! Integrates function F on [A,B] with N subintervals

        implicit none

        real, external      :: F            ! Integrand
        real, intent(in)    :: A            ! Interval start
        real, intent(in)    :: B            ! Interval end
        integer, intent(in) :: N            ! No. of subintervals
        real, intent(out)   :: INTEGRAL     ! Result

        integer         :: I 
        real            :: H
        real, parameter :: R = 4.0/3.0

        H = (B-A)/N

        INTEGRAL = (F(A) + F(B))/3

        do I = 1,N-1,2
            INTEGRAL = INTEGRAL + R*F(A+H*I)
        end do 

        do I = 2,N-1,2
            INTEGRAL = INTEGRAL + 0.5*R*F(A+H*I)
        end do 
        
        INTEGRAL = H*INTEGRAL
        
    end subroutine SIMPSON


end module INTEGRATE

