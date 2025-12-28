subroutine TRAPEZOID(ITER, RAD, EVAL, VAL_T)

    implicit none

    integer, intent(in) :: ITER
    real, intent(in) :: RAD, EVAL 
    real, intent(out) :: VAL_T
    real :: TL, HT, UT, XP

    TL = RAD
    UT = -1.0*RAD
    HT = (TL-UT)/1000000000
    XP = TL - HT
    VAL_T = 0.5*(FN(ITER, RAD, UT) + FN(ITER, RAD, TL))

    do while(XP.gt.UT)
     VAL_T = VAL_T + FN(ITER, RAD, XP)
     XP = XP - HT
    enddo

    VAL_T = HT*VAL_T

    contains

    real function FN(K, PARA, VAR)
            real, intent(in) :: PARA, VAR
            integer, intent(in) :: K
            FN = (PARA**2 - VAR**2)**(0.5*K)
    end function FN

end subroutine TRAPEZOID

program IMPROPER_INTEGRAL

        use FORMULA
        use MONTECARLO

        implicit none

        interface 
                subroutine TRAPEZOID(ITER, RAD, EVAL, VAL_T) 
                implicit none

                real, intent(in) :: RAD, EVAL
                real, intent(out) :: VAL_T
                integer, intent(in) :: ITER
                end subroutine TRAPEZOID
        end interface

        real :: TI, R
        real, allocatable, dimension(:) :: V, V_EXACT, V_ERR
        integer :: N, I

        write(*,*) 'ENTER N'
        read(*,*) N

        write(*,*) 'N = ', N, 'CHECK!'


        R = 1.0

        write(*,*) 'R = ', R, 'CHECK!'

        allocate(V(N))
        allocate(V_EXACT(N))
        allocate(V_ERR(N))
        
        do I = 1,N

                V_EXACT(I) = NBALL_VOLUME(I,R)
                write(*,*) V_EXACT(I)
                write(*,*) 'FOR N = ', I, 'EXACT VOLUME = ', V_EXACT(I)

        if(I.eq.1) then 

                V(I) = 2*R
        else
                call TRAPEZOID(I, R, V_EXACT(I), TI)
                V(I) = V(I-1)*TI
        endif

        V_ERR(I) = V_EXACT(I) - V(I)

        enddo
        
        write(*,*)'----------------------------------------------------------------------------------------------'
        write(*,*) 'VOLUME OF N-SPHERE CALCULATED'
        write(*,*)'----------------------------------------------------------------------------------------------'
        write(*,*) 'INTEGRAL = ', V(N),'EXACT VOLUME = ', V_EXACT(N), 'ERROR(%) = ', 100*abs(V_ERR(N))/V_EXACT(N)
        write(*,*)'----------------------------------------------------------------------------------------------'

        end program IMPROPER_INTEGRAL

