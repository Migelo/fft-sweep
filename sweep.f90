!------------------------------------------------------------------------------
! Program: fft_test_driver
! Description: Sweeping turbulence generator
!              Creates evolving turbulence by sweeping through k-space with decay
!------------------------------------------------------------------------------
program fft_test_driver
    use, intrinsic :: iso_c_binding
    use mpi_f08
    implicit none

    include '/usr/include/fftw3-mpi.f03'

    ! Grid parameters
    integer, parameter :: nx = 256, ny = 256
    real(C_DOUBLE), parameter :: pi = 3.14159265358979d0

    ! Sweep parameters (command line)
    real(C_DOUBLE) :: k_min, k_max
    real(C_DOUBLE) :: t_decay, dt
    integer :: n_steps
    real(C_DOUBLE) :: amplitude
    character(len=64) :: arg
    real(C_DOUBLE) :: t
    real(C_DOUBLE) :: k_old, dk

    ! Arrays
    complex(C_DOUBLE_COMPLEX), pointer :: F_kspace(:,:)   ! k-space working array
    real(C_DOUBLE), pointer :: F_real_local(:,:)          ! accumulated real-space field (local)
    type(C_PTR) :: cdata_k

    ! FFTW MPI variables
    type(C_PTR) :: plan_backward
    integer(C_INTPTR_T) :: alloc_local, local_n0, local_0_start
    integer(C_INT32_T) :: fftw_comm

    ! MPI variables
    integer :: ierr, pid, nproc

    ! Local variables
    integer :: step
    real(C_DOUBLE) :: k_current, decay_factor
    character(len=128) :: filename
    integer, dimension(12) :: seed

    ! Check for help flag first
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        if (trim(arg) == '-h' .or. trim(arg) == '--help') then
            call print_help()
            stop
        end if
    end if

    ! Parse command line: k_min k_max t_decay dt n_steps [amplitude]
    k_min = 2.0d0
    k_max = 32.0d0
    t_decay = 0.2d0
    dt = 0.075d0
    n_steps = 84

    amplitude = 1.0d0
    t = 0.0d0
    k_old = k_min
    dk = 0.0d0

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) k_min
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) k_max
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *) t_decay
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read(arg, *) dt
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        read(arg, *) n_steps
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg)
        read(arg, *) amplitude
    end if

    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    ! Initialize FFTW MPI
    call fftw_mpi_init()

    if (pid == 0) then
        write(*,'(A)') '=== Sweeping Turbulence Generator ==='
        write(*,'(A,F8.2,A,F8.2)') 'Tend: ', n_steps * dt, ' n_periods: ', real(n_steps * dt / (2.0d0 * pi), C_DOUBLE)
        write(*,'(A,I0,A,I0)') 'Grid: ', nx, ' x ', ny
        write(*,'(A,F8.2,A,F8.2)') 'k_min=', k_min, '  k_max=', k_max
        write(*,'(A,F8.2,A,F8.2,A,I0)') 't_decay=', t_decay, '  dt=', dt, '  n_steps=', n_steps
        write(*,'(A,F8.4)') 'amplitude=', amplitude
        write(*,'(A,I0)') 'nproc=', nproc
    end if

    ! Get FFTW communicator
    fftw_comm = int(MPI_COMM_WORLD%MPI_VAL, C_INT32_T)

    ! Get local data size
    alloc_local = fftw_mpi_local_size_2d(int(nx, C_INTPTR_T), int(ny, C_INTPTR_T), &
                                          fftw_comm, local_n0, local_0_start)

    ! Allocate FFTW-aligned memory for k-space
    cdata_k = fftw_alloc_complex(alloc_local)
    call c_f_pointer(cdata_k, F_kspace, [int(ny, C_INTPTR_T), local_n0])

    ! Allocate local real-space accumulator (just real values)
    allocate(F_real_local(ny, local_n0))
    F_real_local = 0.0d0

    ! Create backward FFT plan
    plan_backward = fftw_mpi_plan_dft_2d(int(nx, C_INTPTR_T), int(ny, C_INTPTR_T), &
                                          F_kspace, F_kspace, fftw_comm, &
                                          FFTW_BACKWARD, FFTW_MEASURE)

    ! Set random seed (different per rank for variety)
    seed = 42 + pid * 1000
    call random_seed(put=seed)

    ! Initialize sweep
    k_current = k_min
    decay_factor = exp(-dt / t_decay)

    if (pid == 0) then
        write(*,'(A,F8.4)') 'Decay factor per step: ', decay_factor
        write(*,'(A)') 'Starting sweep...'
    end if

    ! Main time loop
    do step = 1, n_steps
        ! 1. Decay existing real-space field
        F_real_local = F_real_local * decay_factor
        t = t + dt
        ! 2. Generate new k-space contribution for current band
        call fill_kspace_band(F_kspace, int(local_n0), ny, int(local_0_start), nx, &
                              k_current, k_current+1, amplitude)

        ! 3. FFT to real space
        call fftw_mpi_execute_dft(plan_backward, F_kspace, F_kspace)

        ! 4. Add new contribution (normalized) to accumulated field: |F|²
        F_real_local = F_real_local + (real(F_kspace)**2 + aimag(F_kspace)**2) / real(nx * ny, C_DOUBLE)**2

        ! 5. Update k_current (bounce between k_min and k_max)
        ! k_current = k_current + direction * dk
        k_current = (sin(t) + 1) * (k_max - k_min) / 2 + k_min
        dk = abs(k_current - k_old)
        k_old = k_current
        ! if (k_current + dk > k_max) then
        !     direction = -1
        !     k_current = k_max - dk
        ! else if (k_current < k_min) then
        !     direction = 1
        !     k_current = k_min
        ! end if

        ! Save frame
        write(filename, '(A,I4.4,A)') 'frame_', step, '.bin'
        call gather_and_save_real(F_real_local, int(local_n0), int(local_0_start), &
                                   nx, ny, trim(filename), pid, nproc, step, k_current)
    end do

    if (pid == 0) then
        write(*,'(A,I0,A)') 'Done! Generated ', n_steps, ' frames.'
    end if

    ! Cleanup
    call fftw_destroy_plan(plan_backward)
    call fftw_free(cdata_k)
    deallocate(F_real_local)
    call fftw_mpi_cleanup()
    call MPI_Finalize(ierr)

contains

    !--------------------------------------------------------------------------
    ! Fill k-space for a specific k-band [k_lo, k_hi) with random phases
    !--------------------------------------------------------------------------
    subroutine fill_kspace_band(F, local_n0, ny, local_0_start, nx, k_lo, k_hi, amp)
        implicit none
        integer, intent(in) :: local_n0, ny, local_0_start, nx
        real(C_DOUBLE), intent(in) :: k_lo, k_hi, amp
        complex(C_DOUBLE_COMPLEX), intent(out) :: F(ny, local_n0)

        integer :: i, j, gi
        real(C_DOUBLE) :: kmag, phase

        do i = 1, local_n0
            gi = local_0_start + i  ! global i index (1-based)
            do j = 1, ny
                ! Centered k values
                kmag = sqrt(real((gi - nx/2 - 1)**2 + (j - ny/2 - 1)**2, C_DOUBLE))
                if (kmag >= k_lo .and. kmag < k_hi) then
                    call random_number(phase)
                    F(j, i) = amp * exp(cmplx(0.0d0, 2.0d0 * pi * phase, C_DOUBLE_COMPLEX))
                else
                    F(j, i) = cmplx(0.0d0, 0.0d0, C_DOUBLE_COMPLEX)
                end if
            end do
        end do
    end subroutine fill_kspace_band

    !--------------------------------------------------------------------------
    ! Gather real-space field and save to binary file
    !--------------------------------------------------------------------------
    subroutine gather_and_save_real(F_local, local_n0, local_0_start, nx, ny, &
                                     filename, pid, nproc, step, k_current_val)
        implicit none
        real(C_DOUBLE), intent(in) :: F_local(:,:)
        integer, intent(in) :: local_n0, local_0_start, nx, ny, pid, nproc, step
        real(C_DOUBLE), intent(in) :: k_current_val
        character(len=*), intent(in) :: filename

        real(C_DOUBLE), allocatable :: F_global(:,:), buffer(:,:)
        integer, allocatable :: all_local_n0(:), all_start(:)
        integer :: p, i, j, ierr, src_n0, src_start
        type(MPI_Status) :: status

        allocate(all_local_n0(nproc), all_start(nproc))
        call MPI_Gather(local_n0, 1, MPI_INTEGER, &
                        all_local_n0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Gather(local_0_start, 1, MPI_INTEGER, &
                        all_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if (pid == 0) then
            allocate(F_global(ny, nx))
            F_global = 0.0d0

            ! Copy local data from rank 0
            do i = 1, local_n0
                do j = 1, ny
                    F_global(j, all_start(1) + i) = F_local(j, i)
                end do
            end do

            ! Receive from other ranks
            do p = 1, nproc - 1
                src_n0 = all_local_n0(p + 1)
                src_start = all_start(p + 1)
                allocate(buffer(ny, src_n0))
                call MPI_Recv(buffer, src_n0 * ny, MPI_DOUBLE, &
                              p, 0, MPI_COMM_WORLD, status, ierr)
                do i = 1, src_n0
                    do j = 1, ny
                        F_global(j, src_start + i) = buffer(j, i)
                    end do
                end do
                deallocate(buffer)
            end do

            ! Save binary
            call save_real_bin(F_global, nx, ny, filename, step, k_current_val)
            deallocate(F_global)
        else
            call MPI_Send(F_local, local_n0 * ny, MPI_DOUBLE, &
                          0, 0, MPI_COMM_WORLD, ierr)
        end if

        deallocate(all_local_n0, all_start)
    end subroutine gather_and_save_real

    !--------------------------------------------------------------------------
    ! Save real-valued field to binary file
    !--------------------------------------------------------------------------
    subroutine save_real_bin(F, nx, ny, filename, step, k_current)
        implicit none
        integer, intent(in) :: nx, ny, step
        real(C_DOUBLE), intent(in) :: F(ny, nx)
        real(C_DOUBLE), intent(in) :: k_current
        character(len=*), intent(in) :: filename

        integer :: unit_num
        real(C_DOUBLE) :: minv, maxv

        minv = minval(F)
        maxv = maxval(F)

        ! Progress output every 10 steps
        if (mod(step, 10) == 0 .or. step == 1) then
            write(*,'(A,I4,A,F10.4,A,F10.4,A,F8.2)') '  Step ', step, ': min=', minv, ' max=', maxv, ' k=', k_current
        end if

        unit_num = 20
        open(unit=unit_num, file=filename, status='replace', access='stream', form='unformatted')
        write(unit_num) nx, ny
        write(unit_num) k_current
        write(unit_num) F
        close(unit_num)
    end subroutine save_real_bin

    !--------------------------------------------------------------------------
    ! Print help message
    !--------------------------------------------------------------------------
    subroutine print_help()
        implicit none
        write(*,'(A)') 'Sweeping Turbulence Generator'
        write(*,'(A)') ''
        write(*,'(A)') 'Usage: sweep [OPTIONS] k_min k_max t_decay dt n_steps [amplitude]'
        write(*,'(A)') '       mpirun -np N sweep k_min k_max t_decay dt n_steps [amplitude]'
        write(*,'(A)') ''
        write(*,'(A)') 'Options:'
        write(*,'(A)') '  -h, --help    Show this help message and exit'
        write(*,'(A)') ''
        write(*,'(A)') 'Arguments:'
        write(*,'(A)') '  k_min         Minimum wavenumber (default: 1.0)'
        write(*,'(A)') '  k_max         Maximum wavenumber (default: 32.0)'
        write(*,'(A)') '  t_decay       Exponential decay timescale (default: 5.0)'
        write(*,'(A)') '  dt            Timestep (default: 1.0)'
        write(*,'(A)') '  n_steps       Number of frames to generate (default: 100)'
        write(*,'(A)') '  amplitude     Injection amplitude (default: 1.0)'
        write(*,'(A)') ''
        write(*,'(A)') 'Output:'
        write(*,'(A)') '  Generates frame_NNNN.bin files containing |F|^2 field data'
        write(*,'(A)') ''
        write(*,'(A)') 'Example:'
        write(*,'(A)') '  mpirun -np 4 sweep 4 64 5.0 1.0 250 1.0'
        write(*,'(A)') ''
    end subroutine print_help

end program fft_test_driver
