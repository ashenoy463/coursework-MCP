! AYUSH PRAVIN SHENOY
! 24021014
!
! To estimate the gaussian integral by the Trapezoidal method to a given precision
! (by calculating the truncation error from the analytical value) as well as through
! adaptive subintervals

program GINTEG

    use INTEGRATE
    use FORMULA
   
    implicit none

    real :: ALPHA           ! Gaussian parameters
    real :: THR             ! Error tolerance
    
    real :: B               ! Integration Upper Limit
    real    ::  TRUVAL      ! Analytical value
    real    ::  INTEG       ! Numerical result
    real    ::  H           ! Subinterval size
    integer ::  N           ! No of subintervals
    real    ::  DELTA
    
    integer::I

    read(*,*) ALPHA
    read(*,*) THR

    B = 4.0/sqrt(2*ALPHA)           ! 4 Sigma
    TRUVAL = 0.5*sqrt(PI/ALPHA)     ! Half-Integral

    write(*,*) "Alpha               : ", ALPHA
    write(*,*) "Error Threshold (1E): ", -THR
    write(*,*)

    ! Trapezoidal with analytical error bound
    DELTA = 10**(-THR)
    H = sqrt(3.0/sqrt(2*ALPHA))*sqrt(DELTA)
    N = B/H

    write(*,*) "DELTA :", DELTA
    write(*,*) "    H :", H
    write(*,*) "    N :", N
    write(*,*)

    call TRAPEZOID(F,0.0,B,N,INTEG)

    write(*,*) "# ANALYTICAL ERROR BOUND"
    call PRINT_HEADER
    write(*,2) N, INTEG, abs(INTEG-TRUVAL), 100.0*(INTEG-TRUVAL)/TRUVAL
    call PRINT_RESULT(INTEG,TRUVAL)

    ! Trapezoidal with adaptive subintervals
    write(*,*) "# TRAPEZOIDAL METHOD"
    call PRINT_HEADER

    N=1
    do
        call TRAPEZOID(F,0.0,B,N,INTEG)
        write(*,2) N, INTEG, abs(INTEG-TRUVAL), 100.0*(INTEG-TRUVAL)/TRUVAL
        if (abs(2*(INTEG-TRUVAL)) <= 1E-3) exit
        N = N*2
    end do
    call PRINT_RESULT(INTEG,TRUVAL)

    ! Simpson with adaptive subintervals
    write(*,*) "# SIMPSONS 1/3 RULE"
    call PRINT_HEADER

    N=1
    do
        call SIMPSON(F,0.0,B,N,INTEG)
        write(*,2) N, INTEG, abs(INTEG-TRUVAL), 100.0*(INTEG-TRUVAL)/TRUVAL
        if (abs(2*(INTEG-TRUVAL)) <= 1E-3) exit
        N = N*2
    end do
    call PRINT_RESULT(INTEG,TRUVAL)
    !write(*,*) H, 2*INTEG, TRUVAL, abs(2*INTEG - TRUVAL)
    
    2 format (2X,I15,4X,3(F15.10,4X))    ! Values

    contains

        real function F(X)
        
            real, intent(in) :: X

            F = exp(-ALPHA*(X**2))
        end function F

        subroutine PRINT_HEADER()
            write(*,1) "#", "N_INTERVALS", "INTEGRAL", "ABS ERROR", "% ERROR"
            write(*,1) "#", ("---------------", I=1,4)

            1 format (A1,1X,4(A15,4X))    ! Table header
        end subroutine

        subroutine PRINT_RESULT(INTEG,TRUVAL)

            real, intent(in) :: INTEG
            real, intent(in) :: TRUVAL

            write(*,*)
            write(*,*) "Converged value  :", 2*INTEG
            write(*,*) "Analytical value :", 2*TRUVAL
            write(*,*)

        end subroutine

end program GINTEG
