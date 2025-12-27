! AYUSH PRAVIN SHENOY
! 24021014
!
! To sample two uniformly distributed random variables and calculate their Box-Muller
! transform, verifying that the result is normally distributed

program BOXMULLER

    use FORMULA

    implicit none

    integer ::  N               ! Number of samples

    integer          :: I
    real             :: M(2,2)
    real,allocatable :: O(:,:)
    real,allocatable :: R(:,:)

    read (*,*) N

    allocate(O(N,2))
    allocate(R(N,2))
    call RANDOM_NUMBER(O)

    R(:,1) = sqrt(-2*log(O(:,1)))*cos(2*PI*O(:,2))
    R(:,2) = sqrt(-2*log(O(:,1)))*sin(2*PI*O(:,2))

    do I=1,N
        write(*,1) O(I,:), R(I,:)
    end do 

    M(1,1) = sum(R(:,1))/N
    M(1,2) = sum(R(:,1)**2)/N - M(1,1)**2

    M(2,1) = sum(R(:,2))/N
    M(2,2) = sum(R(:,2)**2)/N - M(2,1)**2
    
    write(*,2) "#"
    write(*,2) "#", "X", "Y"
    write(*,2) "#", ("--------------", I =1,2)
    write(*,3) "# Mean : ", M(1,1), M(2,1)
    write(*,3) "# Sdev : ", M(2,1), M(2,2)

    1 format(4(F15.10,4X))
    2 format(A1,8X,A15,4X,A15)
    3 format(A9,F15.10,4X,F15.10)
    
end program BOXMULLER
