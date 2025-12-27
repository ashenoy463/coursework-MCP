subroutine PRINT_LATTICE(LATT,U)

    integer, intent(in) :: U            ! Output unit (6 for stdout)
    integer, intent(in) :: LATT(:,:)    ! Lattice to print
    integer             ::  I
    integer             ::  J
    integer             ::  L

    L = size(LATT(1,:))

    do I = 1,L
        write(U,*) (LATT(I,J), J=1,L)
        write(U,*)
    end do
    !write(U,*) ("------",I=1,2*L)
    write(U,*)

end subroutine

subroutine PRINT_OBSERVABLES(TIME,E,E_SQ,M,M_SQ,U)
        ! Write observables to output file

        integer, intent(in) :: TIME
        real,    intent(in) :: E
        real,    intent(in) :: E_SQ
        real,    intent(in) :: M
        real,    intent(in) :: M_SQ
        integer, intent(in) :: U

        105 format (2X,I8,4X,4(F20.10,4X))    ! Observables
        write(U,105) TIME, E ,E_SQ, M, M_SQ

end subroutine

integer function WRAP(X,L)
    ! Wrap co-ordinate around if needed

    integer, intent(in) :: X
    integer, intent(in) :: L

    if (X .eq. 0) then
        WRAP = X + L
    elseif (X .eq. L+1) then
        WRAP = X - L 
    else
        WRAP = X
    end if

end function

program ISING_MODEL

    implicit none
    
    ! Interfaces and other formalities
    interface
    subroutine PRINT_LATTICE(LAT, U)
        integer, dimension(:,:), intent(in) :: LAT
        integer, intent(in)                 :: U
    end subroutine
    end interface
    integer,external :: WRAP

    ! Model Parameters
    ! ---------------------------------------------------------------
    real    :: T_START  ! Starting Temperature          (units of kB)
    real    :: T_STOP   ! Ending Temperature            (units of kB)
    real    :: T_STEP   ! Temperature step              (units of kB)
    real    :: J        ! Spin-Spin coupling strength   (units of kB)
    real    :: G        ! Spin-Field coupling strength  (units of J)
    integer :: L        ! Lattice size

    ! Numerical Parameters
    ! ---------------------------------------------------------------------
    real    :: CONV_THR = 1E-3         ! Equilibriation Threshold
    integer :: EQ_SAMPLES = 100000     ! Number of post-equilibrium samples
    

    ! Observables
    ! ---------------------------------------------
    real :: M    = 0    ! Net magnetization
    real :: M_SQ = 0    ! Net magnetization squared
    real :: E    = 0    ! Total energy
    real :: E_SQ = 0    ! Total energy sq

    ! Observable Averages
    ! ---------------------------------------------
    real :: MEAN_M    = 0
    real :: MEAN_M_SQ = 0
    real :: MEAN_E    = 0
    real :: MEAN_E_SQ = 0 

    ! Response functions
    ! -------------------------------------------
    real :: CHI         ! Magnetic susceptibility
    real :: C_V         ! Specific heat

    ! Other Variables
    ! ------------------------------------------------------------------------------------
    
    ! IO and Dummy
    character (len=4) :: DUMMY
    integer           :: I,P,Q,S
    logical           :: VERBOSE


    ! Simulation
    integer, allocatable    :: LATTICE(:,:)         ! Current state
    real                    :: T                    ! Current temperature
    real                    :: E_DIFF               ! Difference in energy after M-H step
    real ,dimension(2)      :: R_SF                 ! Spin-flip roll
    real                    :: R_CA                 ! Configuration acceptance roll 
    integer                 :: X_F                  ! Spin-flip X-coord
    integer                 :: Y_F                  ! Spin-flip Y-coord
    integer                 :: TIME                 ! Timestep counter


    ! Model Initialization
    ! -------------------------------------------
    
    ! Read input file from standard input
    read(*,*) T_START
    read(*,*) T_STOP
    read(*,*) T_STEP
    read(*,*) J
    read(*,*) G
    read(*,*) L

    ! Open output files
    open(unit=1,file="./ising.out")
    open(unit=2,file="./observables.out") 
    open(unit=3,file="./state.out")

    ! Write inital inputs
    write(1,*) "INPUT PARAMETERS"
    write(1,*) ("-", I = 1,16)
    write(1,103) "T_STRT = ", T_START
    write(1,103) "T_STOP = ", T_STOP
    write(1,103) "T_STOP = ", T_STEP
    write(1,103) "J = ", J
    write(1,103) "G = ", G
    write(1,104) "L = ", L
    write(1,*)

    ! Write output headers
    write(2,101) "#", "TIME", "E", "E_SQ", "M", "M_SQ"
    write(2,101) "#", "--------", ("--------------------", I=1,4)

    ! Initialize system
    ! ----------------------------------------------
    allocate(LATTICE(L,L))

    LATTICE = +1

    write(1,*) "INITIAL CONFIGURATION"
    write(1,*) "---------------------"
    write(*,*)
    call PRINT_LATTICE(LATTICE,1)
    write(*,*)

    ! Calculate Hamiltonian
    E = - G*sum(LATTICE)     ! Spin-Field Interaction Energy
    do P = 1,L                ! Spin-Spin  Interaction Energy
        do Q = 1,L
            E = E - J*LATTICE(P,Q)*(LATTICE(P          , WRAP(Q+1,L))&  ! Right
                                   +LATTICE(WRAP(P+1,L), Q         ))   ! Down
        end do
    end do

    E_SQ = E**2
    M = sum(LATTICE)
    M_SQ = M**2

    !call PRINT_OBSERVABLES(0,E,E_SQ,M,M_SQ,2)

    write(1,103) "Initial E   =", E
    write(1,103) "Initial M   =", M


    T = T_START
    write(*,106) "#", "TEMP", "CHI_S", "C_V_S", "E_S", "M_S"
    write(*,106) "#", ("--------------------", I=1,5)

    ! Run T : [T_START,T_STOP] Experiment
    ! ----------------------------------------------
    do
        
        MEAN_E = 0
        MEAN_M = 0
        MEAN_E_SQ = 0
        MEAN_M_SQ = 0
        
        TIME = 0 

        ! Timesteps
        ! --------

        do TIME = 1, EQ_SAMPLES
        ! Monte Carlo Cycle
        ! -----------------

            ! Metropolis steps
            do I = 1,SIZE(LATTICE)

                ! Propose configuration

                ! Flip a random spin
                call RANDOM_NUMBER(R_SF)         
                R_SF = R_SF*L
                R_SF = int(ceiling(R_SF))
                X_F = R_SF(1)
                Y_F = R_SF(2)
                S = LATTICE(X_F,Y_F)
                LATTICE(X_F,Y_F) = -1*S

                ! Calculate spin-spin energy change (Up,Down Left Right)
                E_DIFF = +2*J*S*(LATTICE(WRAP(X_F-1,L),Y_F) + LATTICE(WRAP(X_F+1,L),Y_F)&
                                +LATTICE(X_F,WRAP(Y_F-1,L)) + LATTICE(X_F,WRAP(Y_F+1,L)))
                ! Calculate spin-field energy change in non-zero field
                if (abs(G) >= 1E-6) then
                    E_DIFF = E_DIFF - SUM(LATTICE)
                end if

                ! Accept/Reject Config
                if (E_DIFF > 0) then
                    ! Reject with 1 - P(Boltzmann)
                    call RANDOM_NUMBER(R_CA)
                    if (exp(-E_DIFF/T) <= R_CA) then
                        LATTICE(X_F,Y_F) = S
                        E_DIFF = 0.0
                    end if
                end if

                ! Update energy
                E = E+E_DIFF

            end do

            ! Calculate and accumulate observables

            E_SQ = E**2
            M = sum(LATTICE)
            M_SQ = M**2

            MEAN_E    = MEAN_E + E
            MEAN_M    = MEAN_M + M
            MEAN_E_SQ = MEAN_E_SQ + E_SQ
            MEAN_M_SQ = MEAN_M_SQ + M_SQ

        end do

        ! Average out observables and normalize to lattice size

        MEAN_E      = MEAN_E    / EQ_SAMPLES
        MEAN_E_SQ   = MEAN_E_SQ / EQ_SAMPLES
        MEAN_M      = MEAN_M    / EQ_SAMPLES
        MEAN_M_SQ   = MEAN_M_SQ / EQ_SAMPLES

        C_V = (MEAN_E_SQ - MEAN_E**2)
        CHI = (MEAN_M_SQ - MEAN_M**2)

        CHI = CHI/size(LATTICE)
        C_V = C_V/size(LATTICE)
        
        MEAN_E = MEAN_E/size(LATTICE)
        MEAN_M = MEAN_M/size(LATTICE)

        ! Output to log and standard out

        write(1,103) "T = ", T
        write(1,103) "CHI_S =", CHI
        write(1,103) "C_V_S =", C_V
        write(1,103) "-----------------------------"

        write(*,102) T, CHI, C_V, MEAN_E, MEAN_M
        !write(*,*) T, CHI,C_V,E/EQ_SAMPLES,M/EQ_SAMPLES

        T = T+T_STEP

        if (T>T_STOP) stop

    end do

    ! Format : Output
    ! -------------------------------------------
    101 format (A1,1X,A8,4X,5(A20,4X))      ! Table header
    106 format (A1,1X,5(A20,4X))            ! Table header
    102 format (2X,5(F20.10,4X))            ! Average Observables
    103 format (A15,F15.10)                 ! Single variables
    104 format (A15,I8)            

end program ISING_MODEL
