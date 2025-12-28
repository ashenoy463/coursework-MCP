! AYUSH SHENOY
! 24021014
! INTERNAL EXAM (QUESTION 3)

subroutine RAND_EXP(X,LAMBDA)
    ! Sample exponential distribution with mean LAMBDA

    real, intent(inout) :: X
    real, intent(in)    :: LAMBDA

    call RANDOM_NUMBER(X)
    X = 1.0 - X
    X = - log(X)/LAMBDA

end subroutine

subroutine MC_TRIAL(N,M,ENSAV,ENSERR)

    ! Conduct MC trial of M experiments, each with N throws

    integer, intent(in) :: N
    integer, intent(in) :: M
    real, intent(out)   :: ENSAV
    real, intent(out)   :: ENSERR

    integer         :: I,J
    real            ::  R1,R2,T1,T2,P1,P2,CBETA,DENOM
    real,parameter  :: PI = acos(-1.0)
    real            :: ENS(M)

    ENSAV = 0
    ERR = 0

    do J = 1,M

        ESTIMATE = 0
        I = 1

        do while (I<N)
            ! Sample random variables
            call RAND_EXP(R1,4.0)
            call RAND_EXP(R2,4.0)
            call RANDOM_NUMBER(T1)
            call RANDOM_NUMBER(T2)
            call RANDOM_NUMBER(P1)
            call RANDOM_NUMBER(P2)
            T1 = PI*T1
            T2 = PI*T2
            P1 = 2*PI*P1
            P2 = 2*PI*P2
       
            ! Calculate integrand
            CBETA = cos(T1) + cos(T2) + sin(T1)*sin(T2)*cos(P1-P2)
            DENOM = R1**2 + R2**2 - 2*R1*R2*CBETA
            DENOM = sqrt(DENOM)

            ! If not divergent, compute expectation value
            if (DENOM > 1E-2) then
                ESTIMATE = ESTIMATE + (R1**2)*(R2**2)*sin(T1)*sin(T2)/DENOM
                I = I + 1
            end if
            end do
            ENS(J) = (0.25*PI**4)*ESTIMATE/N
    end do 
    ENSAV = sum(ENS)/M

    ENS = ENS - ENSAV
    ENSERR = sqrt(sum(ENS**2)/(real(M)*(real(M)-1)))

end subroutine MC_TRIAL

program CORREL
    
    implicit none
    
    integer :: N ! Number of throws in a experiment
    integer :: M ! Number of experiments in an ensemble
    real    :: INTEG
    real    :: ERR

    N = 100000    ! Reasonably large

    call MC_TRIAL(N,1,INTEG,ERR) 
    write(*,*) "Result for 1 trial : ", INTEG
    write(*,*)

    ! Run MC simulation
    write(*,"(A4,4X,3(A15,4X))") "M", "1/SQRT(M)", "INTEGRAL", "ERROR"
    write(*,"(A4,4X,3(A15,4X))") "----", ("----------",M=1,3)

    do M = 2,50
    
        call MC_TRIAL(N,M,INTEG,ERR) 
        write(*,"(I4,4X,3(F15.10,4X))") M, 1.0/sqrt(real(M)), INTEG, ERR

    end do

end program CORREL
