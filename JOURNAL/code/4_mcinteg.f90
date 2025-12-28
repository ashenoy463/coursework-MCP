! AYUSH PRAVIN SHENOY
! 24021014
!
! To estimate the value of a gaussian integral to within a given precision using
! uniform-sample Monte Carlo integration along with ensemble-averaging


program MCINTEG

    use FORMULA
    
    implicit none
    
    real    ::  A           ! Gaussian parameter
    real    ::  T           ! Interval Endpoint
    integer ::  N           ! Number of samples
    real    ::  THR         ! Convergence threshold
    integer ::  M           ! Ensemble size

    real::INTEG
    integer::I
    real::TRUVAL
    real,allocatable::R(:)


    read(*,*) A
    read(*,*) T
    read(*,*) N
    read(*,*) THR
    read(*,*) M

    write(*,*) "Error Threshold (1E): ", -THR
    write(*,*)
    write(*,1) "ENSEMBLE SIZE", "INTEGRAL", "ABS ERROR", "PERCENT ERROR"
    write(*,1) ("----------------", I=1,4)

    allocate(R(N))

    THR = 10**(-THR)
    TRUVAL = 0.5*sqrt(PI/A)*erf(T)

    do
        INTEG = 0
        do I=1,M
            call RANDOM_NUMBER(R)
            R = R*T                              ! Scale dist to T
            R = GAUSSIAN(R,A)                    ! Evaluate function
            INTEG = INTEG + sum(R)*T/N           ! Calculate integral
        end do
        INTEG = INTEG/M

        write(*,2) M, INTEG, INTEG-TRUVAL, 100.0*(INTEG - TRUVAL)/TRUVAL

        ! Check convergence
        if (abs(INTEG-TRUVAL) .lt. THR) exit

        ! Adjust ensemble size
        M = 2*M

    end do

    write(*,*) 
    write(*,*) "The converged value is: ", INTEG
    write(*,*) "The tabulated value is: ", TRUVAL

    1 format (A15,4X,3(A15,4X))
    2 format (I15,4X,3(F15.10,4X))
    
end program MCINTEG
