! This program computes the mean squared displacement of a list of sites during a trajectory in time.
! It is written by Maximilien Levesque while in postdoc in the group of Mathieu Salanne
! at UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France.

program meanSquaredDisplacement

    implicit none
    character(len("msd.out")) :: outputFile = "msd.out"
    integer :: Nat
    integer :: i, nbTimeStepsInTraj, iostat, dt, t, nt
    integer, parameter :: x=1, y=2, z=3
    double precision, dimension(:,:,:), allocatable :: r ! position of site i at timestep t
    double precision :: tmp, lx, ly, lz
    double precision, dimension(:), allocatable :: msd
    character(len=300) :: arg, trajectoryFileName

    ! read all arguments that MUST be given at execution
    call readArguments(lx,ly,lz,Nat,trajectoryFileName)

    ! deduce the number of timesteps in the trajectory from the number of lines in the trajectory file
    nbTimeStepsInTraj = NbOfLinesInTraj()/Nat

    print*,'You have' ,nbTimeStepsInTraj,' time steps in your trajectory file ',trim(adjustl(trajectoryFileName))
    print*,'Please be patient. Everything seems fine... Multiorigin effect! ;)'

    ! allocate consequently what will store the MSD(t)
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

    ! MSD(t) will be written in file unit 11
    open(11,file=outputfile)

    ! compute msd(dt)= <min{|r_i(t)-r_i(t+dt)|}_{PBC(r_i(t+dt))}²>_{i,t}
    msd=0.d0
    do dt = 1, nbTimeStepsInTraj-1
        nt = nbTimeStepsInTraj-dt
        do t = 1, nt
            tmp = 0.d0
            do i = 1, Nat
                tmp = tmp + min(    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,ly,lz] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,ly,0.d0] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,ly,-lz] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,0.d0,lz] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,0.d0,0.d0] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,0.d0,-lz])**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,-ly,lz] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,-ly,0.d0] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [lx,-ly,-lz])**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,ly,lz] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,ly,0.d0] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,ly,-lz])**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,0.d0,lz] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,0.d0,0.d0] )**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,0.d0,-lz])**2)   ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,-ly,lz] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,-ly,0.d0] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [0.d0,-ly,-lz])**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,ly,lz] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,ly,0.d0] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,ly,-lz])**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,0.d0,lz] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,0.d0,0.d0] )**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,0.d0,-lz])**2)  ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,-ly,lz] )**2) ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,-ly,0.d0] )**2) ,&
                                    sum(( r(i,t,:) - r(i,t+dt,:) - [-lx,-ly,-lz])**2)  &
                                )
            end do
            msd(dt) = msd(dt) + tmp/dble(Nat)
!~             sum(  (r(:,t,x)-r(:,t+dt,x))**2 + (r(:,t,y)-r(:,t+dt,y))**2 + (r(:,t,z)-r(:,t+dt,z))**2   )/dble(Nat) ! average over sites i
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

    function NbOfLinesInTraj()
        integer :: NbOfLinesInTraj
        call opentraj
        ! computes the number of lines in traj.in and deduces the number of timesteps
        NbOfLinesInTraj = -1
        do while (iostat == 0)
            read(10,*,iostat=iostat) tmp
            NbOfLinesInTraj = NbOfLinesInTraj + 1
        end do
        call closetraj
    end function
    
    subroutine readArguments(lx,ly,lz,Nat,trajectoryFileName)
        double precision, intent(out) :: lx, ly, lz
        integer, intent(out) :: Nat
        character (len=*), intent(out) :: trajectoryFileName
        
        call get_command_argument(1,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
        end if
        read(arg,*) lx
    
        call get_command_argument(2,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
        end if
        read(arg,*) ly
    
        call get_command_argument(3,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
        end if
        read(arg,*) lz
    
        call get_command_argument(4,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
        end if
        read(arg,*) Nat
    
        call get_command_argument(5,arg,status=i)
        if( i < 0 ) then
            stop "STOP. The length of the argument is too big for me :( "
        else if ( i > 0 ) then
            stop "Argument retrieval failed. You should execute the program with the number of atoms as argument, e.g. in ./msd 10 "
        end if
        trajectoryFileName = trim(adjustl(arg))
    
    end subroutine
  
end program
