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

    use IO

    implicit none
    
    ! Interfaces and other formalities
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
    write(1,203) "T_STRT = ", T_START
    write(1,203) "T_STOP = ", T_STOP
    write(1,203) "T_STOP = ", T_STEP
    write(1,203) "J = ", J
    write(1,203) "G = ", G
    write(1,204) "L = ", L
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
    write(1,*)
    call PRINT_LATTICE(LATTICE,1)
    write(1,*)

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

    write(1,203) "Initial E   =", E
    write(1,203) "Initial M   =", M


    T = T_START
    write(*,102) "#", "TEMP", "CHI_S", "C_V_S", "E_S", "M_S"
    write(*,102) "#", ("--------------------", I=1,5)

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
            do I = 1,(L**2)

                ! Propose configuration

                ! Flip a random spin
                call RANDOM_NUMBER(R_SF)         
                R_SF = R_SF*L
                R_SF = int(ceiling(R_SF))
                X_F = R_SF(1)
                Y_F = R_SF(2)
                S = LATTICE(X_F,Y_F)
                LATTICE(X_F,Y_F) = -1*S

                M = M - 2*S

                !write(3,*) "TEMP", T, "TIME", TIME
                !write(3,*) I, "PROP", X_F, Y_F, S, -S

                ! Calculate spin-spin energy change (Up,Down Left Right)
                E_DIFF = +2*J*S*(LATTICE(WRAP(X_F-1,L), Y_F          ) + LATTICE(WRAP(X_F+1,L), Y_F          )&
                                +LATTICE(X_F,           WRAP(Y_F-1,L)) + LATTICE(X_F,           WRAP(Y_F+1,L)))
                ! Calculate spin-field energy change in non-zero field
                if (abs(G) >= 1E-6) then
                    E_DIFF = E_DIFF + 2*G*S
                end if

                ! Accept/Reject Config
                if (E_DIFF > 0) then
                    ! Reject with 1 - P(Boltzmann)
                    call RANDOM_NUMBER(R_CA)
                    if (exp(-E_DIFF/T) <= R_CA) then
                        !write(3,*) "POS REJECT", E_DIFF,R_CA, exp(-E_DIFF/T)

                        LATTICE(X_F,Y_F) = S
                        E_DIFF = 0.0
                        M = M + 2*S
                    else
                        !write(3,*) "POS ACCEPT", E_DIFF,R_CA, exp(-E_DIFF/T)
                    end if
                else
                        !write(3,*) "NEG ACCEPT", E_DIFF,R_CA, exp(-E_DIFF/T)
                end if

                ! Update energy
                E = E+E_DIFF

                !call PRINT_LATTICE(LATTICE, 3)
                !write(3,*) E
                !call PRINT_OBSERVABLES(TIME,LATTICE,E,3)
                !write(3,*) "-------------------------"
            end do

            ! Calculate and accumulate observables

            E_SQ = E**2
            M = sum(LATTICE)
            M_SQ = M**2

            MEAN_E    = MEAN_E + E
            MEAN_M    = MEAN_M + M
            MEAN_E_SQ = MEAN_E_SQ + E_SQ
            MEAN_M_SQ = MEAN_M_SQ + M_SQ

            ! Data and debug Outputs

            !call PRINT_OBSERVABLES(TIME,E,E_SQ,M,M_SQ,2)
            
            !write(3,*) "-------------------------"
            !write(3,*) "AFTER TIME", TIME
            !write(3,*) 
            !call PRINT_LATTICE(LATTICE,3)
            !call PRINT_OBSERVABLES(TIME,E,E_SQ,M,M_SQ,M_S,3)
            !write(3,*) "-------------------------"

        end do

        ! Average out observables and normalize to lattice size

        MEAN_E      = MEAN_E    / EQ_SAMPLES
        MEAN_E_SQ   = MEAN_E_SQ / EQ_SAMPLES
        MEAN_M      = MEAN_M    / EQ_SAMPLES
        MEAN_M_SQ   = MEAN_M_SQ / EQ_SAMPLES

        C_V = (MEAN_E_SQ - MEAN_E**2)/((L**2)*T**2)
        CHI = MEAN_M_SQ - MEAN_M**2/((L**2)*T)

        MEAN_E = MEAN_E/(L**2)
        MEAN_M = MEAN_M/(L**2)

        ! Output to log and standard out

        write(*,201) T, CHI, C_V, MEAN_E, MEAN_M

        write(3,*) "T = " , T
        write(3,*) "---------------------"
        write(3,*)
        call PRINT_LATTICE(LATTICE,3)

        T = T+T_STEP

        if (T>T_STOP) exit

    end do

    ! Exit cleanly

    do I=1,3
        close(I)
    end do

    ! Format : Output
    ! -------------------------------------------
    101 format (A1,1X,A8,4X,5(A20,4X))      ! Table header
    102 format (A1,1X,5(A20,4X))            ! Table header
    201 format (2X,5(F20.10,4X))            ! Average Observables
!   202 format (2X,I8,4X,4(F20.10,4X))      ! Time Series (COMMON/IO.f90)
    203 format (A15,F15.10)                 ! Single variables
    204 format (A15,I8)            

end program ISING_MODEL
