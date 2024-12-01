module structured_mesh_mod
    use precision
    use data_mod
    implicit none

    type, public :: StructCelleType 

        integer  :: code                              ! code (1,2...,nb of celles) of the celle
        logical  :: boudary                           ! if TRUE boundary celle; else FALSE 
        real(pr) :: center_x_coord, center_y_coord    ! coord center celle
        real(pr) :: x_left, x_right                   ! coord x left and right of edge
        real(pr) :: y_up, y_down                      ! coord y up and down of edge
        real(pr) :: r_p, m_p, T_p, m_sel, vx_p, vy_p  ! center celle radius, mass, Temperature and value
        real(pr), dimension(:,:,:,:), allocatable :: SOL

    end type StructCelleType

    
    private :: mesh1D
    private :: mesh2D

    public :: init_mesh

    public :: meshbuild
    interface meshbuild
        module procedure mesh1D
        module procedure mesh2D
    end interface meshbuild
    
contains

    subroutine init_mesh(data, kcelle)

        !In
        type(DataType), intent(in) :: data

        !InOut
        type(StructCelleType), dimension(:), allocatable, intent(inout) :: kcelle
        
        !Local
        integer :: i, j, cd
        real(pr), dimension(data%Nx) :: x_1D
        real(pr), dimension(data%Nx,data%Ny) :: x_2D, y_2D

        if (data%dim == 1) then
            allocate(kcelle(data%Nx-1))
            call meshbuild(data, x_1D)
            do i=1,data%n_celle
                kcelle(i)%code           = i
                kcelle(i)%center_x_coord = 0.5_pr*(x_1D(i) + x_1D(i+1))
                kcelle(i)%center_y_coord = 0.0
                kcelle(i)%x_left         = x_1D(i)
                kcelle(i)%x_right        = x_1D(i+1)
                kcelle(i)%y_up           = 0.0
                kcelle(i)%y_down         = 0.0

                if (i==1 .or. i==data%n_celle) then
                    kcelle(i)%boudary = .true.
                else
                    kcelle(i)%boudary = .false.
                end if
            end do

        else if (data%dim == 2) then
            allocate(kcelle((data%Nx-1)*(data%Ny-1)))
            call meshbuild(data, x_2D, y_2D)
            cd = 0
            do j=1,data%Ny-1
                do i=1,data%Nx-1
                    print*, "i*j=", i*j
                    cd                        = cd +1
                    kcelle(cd)%code           = cd
                    kcelle(cd)%center_x_coord = (x_2D(i,j) + x_2D(i+1,j))*0.5_pr
                    kcelle(cd)%center_y_coord = (y_2D(i,j) + y_2D(i,j+1))*0.5_pr
                    kcelle(cd)%x_left         = x_2D(i,j)
                    kcelle(cd)%x_right        = x_2D(i+1,j)
                    kcelle(cd)%y_down         = y_2D(i,j)
                    kcelle(cd)%y_up           = y_2D(i,j+1)

                    if (i==1 .or. i==data%Nx-1 .or. j==1 .or. j==data%Ny-1) then
                        kcelle(cd)%boudary = .true.
                    else
                        kcelle(cd)%boudary = .false.
                    end if
                end do
            end do

        else
            print*, "Error: Dimension can't be greater than 2"
            stop

        end if

    end subroutine init_mesh

    subroutine mesh1D(data, x_coord)

        !In
        type(DataType), intent(in) :: data

        !Out
        real(pr), dimension(:), intent(out) :: x_coord

        !Local
        integer :: i

        do i=1,size(x_coord)
            x_coord(i) = data%dx * (i-1) + data%x_min
        end do
        
    end subroutine mesh1D

    subroutine mesh2D(data, x_coord, y_coord)

        !In
        type(DataType), intent(in) :: data

        !Out
        real(pr), dimension(:,:), intent(out) :: x_coord, y_coord

        !Local
        integer :: i,j

        do i=1,data%Nx
            do j=1,data%Ny
                x_coord(i,j) = data%dx * (i-1) + data%x_min
                y_coord(i,j) = data%dy * (j-1) + data%y_min
            end do
        end do

        
    end subroutine mesh2D


    
end module structured_mesh_mod