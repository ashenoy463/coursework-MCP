real function FC(X)
    FC = exp(-x**2)
end function FC

program TEST

    use INTEGRATE

    implicit none

    real, external::FC
    real:: T

    call SIMPSON(FC,0.0,1.0,100,T)

    write(*,*) T
    
end program TEST
