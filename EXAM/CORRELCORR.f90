real function F(R1,R2,THETA1,THETA2,PHI1,PHI2)

    real, intent(in)::R1
    real, intent(in)::R2
    real, intent(in)::THETA1
    real, intent(in)::THETA2
    real, intent(in)::PHI1
    real, intent(in)::PHI2

    real :: CBETA,DENOM

    CBETA = cos(THETA1) + cos(THETA2) + sin(THETA1)*sin(THETA2)*cos(PHI1-PHI2)
    DENOM = R1**2 + R2**2 - 2*R1*R2*CBETA
    DENOM = sqrt(DENOM)
    F = exp(-4.0*(R1+R2))
    F = F*(R1**2)
    F = F*(R2**2)
    F = F*sin(THETA1)*sin(THETA2)/DENOM

end function

program CORREL

    implicit none

    real::R1
    real::R2
    real::THETA1
    real::THETA2
    real::PHI1
    real::PHI2

    real, external::F
    integer::I,M,N,J
    real::UL_R
    real :: INTEG,RES
    real,parameter::PI = 4.0*atan(1.0)
    real::T


    UL_R = 10
    
    T = 2*PI*PI*UL_R

    N = 100

    !write(*,*) "Error Threshold (1E): ", -THR
    !write(*,*)
    write(*,1) "#", "M", "INTEGRAL"
    write(*,1) "#", ("----------------", I=1,2)

    !THR = 10**(-THR)

    !do I = 1, 1000

    !    write(*,*) (I-1), F(real(I-1),real(I-1),PI/2,PI/2,0.0,0.0)

    !end do

    do M = 1,10

        RES = 0

        do J = 1,M

            INTEG = 0
            do I=1,N

                do
                call RANDOM_NUMBER(R1)
                call RANDOM_NUMBER(R2)
                call RANDOM_NUMBER(THETA1)
                call RANDOM_NUMBER(THETA2)
                call RANDOM_NUMBER(PHI1)
                call RANDOM_NUMBER(PHI2)

                R1 = R1 * UL_R
                R2 = R2 * UL_R

                THETA1 = THETA1 * PI
                THETA2 = THETA2 * PI

                PHI1 = PHI1 * 2.0 * PI
                PHI2 = PHI2 * 2.0 * PI
                
                if (abs(R1-R2) >= 1E-2) exit

                end do

                INTEG = INTEG + F(R1,R2,THETA1,THETA2,PHI1,PHI2)

                !write(*,*) R1,R2,THETA1,THETA2,PHI1,PHI2, F(R1,R2,THETA1,THETA2,PHI1,PHI2)

            end do

            RES = RES + (T**2)*INTEG/N

        end do

        write(*,2) M, RES/M

    end do

    1 format (1A,1X,2(A15,4X))
    2 format (2X,I15,4X,1(F15.10,4X))
    
end program CORREL
