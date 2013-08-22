! This program computes the mean squared displacement of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

program meanSquaredDisplacement

    implicit none
    character(len("traj.in")) :: trajectoryFileName = "traj.in"
    character(len("msd.out")) :: outputFile = "msd.out"
    integer :: Nat
    integer :: i, nbTimeStepsInTraj, iostat, dt, t, nline, nt
    integer, parameter :: x=1, y=2, z=3
    double precision, dimension(:,:,:), allocatable :: r ! position of site i at timestep t
    double precision :: tmp ! dummy
    double precision, dimension(:), allocatable :: msd
    character(len=300) :: arg, trajectoryFileName
    
    call get_command_argument(1,arg,status=i)
    if( i < 0 ) then
        stop "STOP. The length of the argument is too big for me :( "
    else if ( i > 0 ) then
        stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
    end if
    read(arg,*) Nat
    print*, Nat*2

    call get_command_argument(2,arg,status=i)
    if( i < 0 ) then
        stop "STOP. The length of the argument is too big for me :( "
    else if ( i > 0 ) then
        stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
    end if
    print*, trim(adjustl(arg))

    
STOP    
    
    call opentraj
    ! computes the number of lines in traj.in and deduces the number of timesteps
    nline = -1
    do while (iostat == 0)
        read(10,*,iostat=iostat) tmp
        nline = nline + 1
    end do
    call closetraj
    nbTimeStepsInTraj = nline/Nat

    ! allocate consequently
    allocate( msd(nbTimeStepsInTraj-1), source=0.d0 )
        
    ! read positions of all sites i at all timesteps t
    allocate( r(Nat,nbTimeStepsInTraj,x:z) )
    call opentraj
    do t = 1, nbTimeStepsInTraj
        do i = 1, Nat
            read(10,*) r(i,t,x), r(i,t,y), r(i,t,z)
        end do
    end do
    call closetraj

    open(11,file=outputfile)

    ! compute msd(dt)= <|r_i(t)-r_i(t+dt)|²>_{i,t}
    do dt = 1, nbTimeStepsInTraj-1
        nt = nbTimeStepsInTraj-dt
        do t = 1, nt
            msd(dt) = msd(dt) +sum(  (r(:,t,x)-r(:,t+dt,x))**2 + (r(:,t,x)-r(:,t+dt,y))**2 + (r(:,t,x)-r(:,t+dt,z))**2   )/dble(Nat) ! average over sites i
        end do
        msd(dt) = msd(dt) /dble(nt) ! average over reference time t
        write(11,*) dt, msd(dt)
    end do
    
    contains

    subroutine opentraj
        call inquireFileExistence(trajectoryFileName)
        ! read positions
        open(10, file=trajectoryFileName,status='old',iostat=iostat)
        if (iostat /= 0) then
            write (*,*) 'File open failed for',trajectoryFileName
            write (*,*) 'Error code is ', iostat
            stop
        end if
    end subroutine
    
    subroutine closetraj
        close(10)
    end subroutine

    subroutine inquireFileExistence(fileName)
        character(len=*), intent(in) :: fileName
        integer, parameter :: stderr = 0
        logical :: exist
        inquire(file=fileName, exist=exist)
        if( .not. exist) then
            write(stderr,*) "YOUR ERROR (not mine ;): The file ",fileName," does not exist. It should."
            stop
        end if
    end subroutine

  
end program
