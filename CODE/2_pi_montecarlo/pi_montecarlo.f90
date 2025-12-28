! AYUSH PRAVIN SHENOY
! 24021014
!
! To estimate the value of π to within a given precision using the acceptance-rejection
! method along with ensemble-averaging

program PIMONTEC

    use iso_fortran_env, only: IK => int32, RK => real64
    
    use MONTECARLO

    implicit none
    
    integer ::  N        ! ensemble size
    real    ::  THR      ! ERR < 1.0E-(THR)
    integer ::  M        ! Number of throws

    integer ::  J
    real    ::  GUESS
    real    ::  ENSAV

    double precision,parameter::PI=4*atan(1.d0) 

    read(*,*) THR
    read(*,*) N
    read(*,*) M

    write(*,*) "Error Threshold (1E): ", -THR
    write(*,*)
    write(*,1) "N_THROWS", "GUESS", "ABS ERROR", "PERCENT ERROR"
    write(*,1) ("----------------", J=1,4)

    THR = 10**(-THR)

    do

        ! Populate ensemble
        ENSAV=0
        do J=1,N
            call ACCREJ_SPHERE(1.0,2,M,GUESS)
            ENSAV = ENSAV + GUESS
        end do
        ENSAV = ENSAV/N
        
        ! Output results
        write(*,2) M, ENSAV, ENSAV-PI, 100*(ENSAV-PI)/PI

        ! Check for convergence
        if (abs(ENSAV - PI) .lt. THR) exit
    
        M = 10*M

    end do

    write(*,*)
    write(*,*) "Converged value of PI is: ", ENSAV

1 format (A15,4X,3(A15,4X))
2 format (I15,4X,3(F15.10,4X))

end program PIMONTEC
