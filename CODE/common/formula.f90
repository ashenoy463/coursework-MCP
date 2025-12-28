module FORMULA

    use iso_fortran_env, only: IK => int32, RK => real64

    implicit none
    public::NBALL_VOLUME
    double precision,parameter::PI = 4*atan(1.d0)
    
    contains

    real function NBALL_VOLUME(N,R)
        ! Compute the volume of a N-dimensional hypersphere of radius R

        integer,intent(in)::N
        real,intent(in)::R

        NBALL_VOLUME = (R**N)*(PI**(0.5*N))/gamma((0.5*N)+1)
        
    end function NBALL_VOLUME

    pure function GAUSSIAN(X,K) result(FUNC)

        real,dimension(:),intent(in)::X
        real,intent(in)::K
        real,dimension(SIZE(X))::FUNC

        FUNC(:) = exp(-K*(X(:)**2))

    end function GAUSSIAN

end module FORMULA

