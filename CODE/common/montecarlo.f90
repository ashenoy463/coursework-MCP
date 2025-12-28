module MONTECARLO

    use iso_fortran_env, only: IK => int32, RK => real64

    implicit none
    public::ACCREJ_SPHERE
    
    contains

    subroutine ACCREJ_SPHERE(RAD,DIMM,NRAND,ESTVOL)
        ! Compute volume of DIMM-dimensional hypersphere of radius RAD
        ! by acceptance-rejection method
       
        implicit none 

        real,intent(in)::RAD        ! Scaling factor
        integer,intent(in)::DIMM    ! Dimension of RN vector
        integer,intent(in)::NRAND   ! Number of throws
        real,intent(out)::ESTVOL    ! Estimated Volume

        integer::I
        real::NIN
        double precision,allocatable::R(:,:)

        ! Generate and shift RN vector
        allocate(R(NRAND,DIMM))
        call RANDOM_NUMBER(R)
        R = (R - 0.5)*(2*RAD)

        ! Count internal points
        NIN = 0
        do I=1,NRAND
            if (sum(R(I,:)**2) <= RAD**DIMM) then
                NIN = NIN+1 
            endif
        end do
        
        ! Compute volume
        ESTVOL = (NIN/NRAND)*(2**DIMM)

        deallocate(R)

    end subroutine ACCREJ_SPHERE

end module MONTECARLO

