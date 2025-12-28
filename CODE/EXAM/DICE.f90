program DICE

    implicit none

    real::ROLL(4)
    integer::I,J,N
    real::M,SQM
    
    read(*,*) N

    M = 0
    SQM = 0 

    do I=1,N
        call RANDOM_NUMBER(ROLL)
        ROLL = ceiling(ROLL*6)
        write(*,*) (int(ROLL(J)), J=1,4)
        M = M + ROLL(1)
        SQM = SQM + ROLL(1)**2
    end do
    
    M = M/N
    SQM = SQM/N

    write(*,*)
    write(*,*) "Mean: ", M
    write(*,*) "Uniform mean: ", 3.5
    write(*,*) "Var :" , SQM - M**2
    write(*,*) "Uniform Var", 2.917

end program DICE
