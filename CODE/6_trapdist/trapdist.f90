program TRAP_DIST

    implicit none


    integer :: N
    real    :: A
    real    :: B
    real    :: C

    real,allocatable    ::  U(:)
    real,allocatable    ::  R(:)
    real                ::  D
    integer             ::  I

    ! Read parameters
    read (*,*) N
    read (*,*) A
    read (*,*) B
    read (*,*) C


    allocate(U(N))
    allocate(R(N))
    call RANDOM_NUMBER(U)

    D = A/(2*B)

    ! Calculate inverse CDF
    do I = 1,N
 
        if ( U(I) >= 0 .and. U(I) <= D ) then
            R(I) = C + sqrt(2.0*A*B*U(I))
        elseif ( U(I) > D .and. U(I) < (1- D)) then
            R(I) = B*U(I) + 0.5*A + C
        elseif ( U(I) > (1-D) .and. U(I) <= 1.0 ) then
            R(I) = A+B+C - sqrt(2*A*B*(1-U(I)))
        endif

        write(*,1) U(I), R(I)

    end do

    1 format(2(F15.10,4X))
    
end program TRAP_DIST
