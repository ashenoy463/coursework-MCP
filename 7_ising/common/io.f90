module IO

    implicit none

    public  ::  PRINT_LATTICE
    public  ::  PRINT_OBSERVABLES

contains

    subroutine PRINT_LATTICE(LATT,U)

        integer, intent(in) :: U            ! Output unit (6 for stdout)
        integer, intent(in) :: LATT(:,:)    ! Lattice to print
        integer             :: I
        integer             :: J
        integer             :: L_I,L_J
        
        character(len=32) :: FRMT
        character(len=32) :: BUF
       
        L_I = size(LATT(:,1))
        L_J = size(LATT(1,:))

        ! Generate format
        write(BUF,"(I32)") L_J
        FRMT = "("//trim(BUF)//"(I2,2X))"

        do I = 1,L_I
            write(U,FRMT) (LATT(I,J), J=1,L_J)
            write(U,*)
        end do
        write(U,*)

    end subroutine

    subroutine PRINT_OBSERVABLES(TIME,E,E_SQ,M,M_SQ,U)
        ! Write observables to output file

        integer, intent(in) :: TIME
        real,    intent(in) :: E
        real,    intent(in) :: E_SQ
        real,    intent(in) :: M
        real,    intent(in) :: M_SQ
        integer, intent(in) :: U

        write(U,202) TIME, E ,E_SQ, M, M_SQ

        202 format (2X,I8,4X,4(F20.10,4X))      ! Time Series (COMMON/IO.f90)

    end subroutine

end module IO
