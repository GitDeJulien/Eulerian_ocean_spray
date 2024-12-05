module structured_mesh_mod

    use precision
    use data_mod
    implicit none

    type, public :: MeshType
        real(pr), dimension(:), allocatable :: x_tab
        real(pr), dimension(:), allocatable :: r_tab
        real(pr), dimension(:), allocatable :: T_tab
        real(pr), dimension(:), allocatable :: m_tab
        real(pr), dimension(:), allocatable :: vx_tab

        real(pr), dimension(:,:,:,:), allocatable :: SOL
    end type MeshType


    private :: logspace
    public  :: init_mesh, free_mesh
    
contains

    subroutine init_mesh(data, mesh)

        !In
        type(DataType), intent(in) :: data

        !InOut
        type(MeshType), intent(inout) :: mesh

        !Local
        integer  :: i

        allocate(mesh%x_tab(data%Nx))
        allocate(mesh%r_tab(data%N_r))
        allocate(mesh%vx_tab(data%N_vx))
        allocate(mesh%m_tab(data%N_r))
        allocate(mesh%T_tab(data%N_T))

        do i=1,data%Nx
            mesh%x_tab(i) = data%dx * (i-1) + data%x_min !Space
        end do

        mesh%r_tab = logspace(data%r_min, data%r_max, data%N_r)

        do i=1,data%N_vx
            mesh%vx_tab(i) = data%vx_min + i*data%dvx !Velocity
        enddo

        mesh%m_tab = 4.0_pr/3.0_pr*pi*data%rho_p*mesh%r_tab !Mass (radius dependant)

        do i=1,data%N_T
            mesh%T_tab(i) = data%T_min + i*data%dT !Temperature
        enddo

    end subroutine init_mesh

    subroutine free_mesh(mesh)

        !In
        type(MeshType), intent(inout) :: mesh

        deallocate(mesh%x_tab)
        deallocate(mesh%r_tab)
        deallocate(mesh%m_tab)
        deallocate(mesh%vx_tab)
        deallocate(mesh%T_tab)

        deallocate(mesh%SOL)

    end subroutine free_mesh

    function logspace(start, stop, num) result(log_vector)

        !In
        real(pr), intent(in) :: start, stop
        integer, intent(in) :: num

        !Out
        real(pr), dimension(num) :: log_vector

        !Local
        integer :: i
        real(pr) :: log_start, log_stop, log_step
    
        ! Calculate the logarithms of the start and stop values
        log_start = log10(start)
        log_stop = log10(stop)
    
        ! Calculate the step size in log space
        log_step = (log_stop - log_start) / (num - 1)
    
        ! Generate the logarithmically spaced vector
        do i = 1, num
          log_vector(i) = 10.0_pr ** (log_start + (i - 1) * log_step)
        end do

    end function logspace
    
end module structured_mesh_mod