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


    public :: init_mesh, free_mesh
    
contains

    subroutine init_mesh(data, mesh)

        !In
        type(DataType), intent(in) :: data

        !InOut
        type(MeshType), intent(inout) :: mesh

        !Local
        integer  :: i, num
        real(pr) :: factor

        num = 10
        factor = (data%r_max/data%r_min)**(1.0/(num-1))

        allocate(mesh%x_tab(data%Nx))
        allocate(mesh%r_tab(data%N_r))
        allocate(mesh%vx_tab(data%N_vx))
        allocate(mesh%m_tab(data%N_r))
        allocate(mesh%T_tab(data%N_T))

        do i=1,data%Nx
            mesh%x_tab(i) = data%dx * (i-1) + data%x_min !Space
        end do

        do i=1,data%N_r
            mesh%r_tab(i) = data%r_min*factor**i !Radius (log10 space)
        enddo

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


    
end module structured_mesh_mod