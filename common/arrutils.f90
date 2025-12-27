
module ARRUTILS

    use iso_fortran_env, only: IK => int32, RK => real64

    implicit none
    public::LINSPACE

contains

    pure function LINSPACE(A,B,D,H)
        
        real,intent(in)::A
        real,intent(in)::B
        integer,optional,intent(in)::D
        real,optional,intent(in)::H

        integer::I
        integer::N
        real::EPS
        real,allocatable::LINSPACE(:)

        if(present(H)) then
            EPS = H
            N = nint((B-A)/H)+1
        else if (present(D)) then
            EPS = (B-A)/D
            N = D+1
        end if
        
        allocate(LINSPACE(N))

        do I=0,N-1
            LINSPACE(I+1) = A + I*EPS
        end do
        
    end function LINSPACE

end module ARRUTILS

