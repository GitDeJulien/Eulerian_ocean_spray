module init_mod

    use structured_mesh_mod
    use sample_mod

    implicit none

    !private :: initialize1D
    public :: initial_step
    
contains

    subroutine initial_step(data, celles)

        !In
        type(DataType), intent(in) :: data

        !Out
        type(StructCelleType), dimension(:), intent(inout) :: celles

        !Local
        integer :: p

        !celles already allocated

        if (data%dim==1) then
            do p=1,data%n_celle
                allocate(celles(p)%SOL(data%N_r, data%N_vx, data%N_m, data%N_T))
                call initialize1D(data, celles(p))
            enddo
            !TODO > le cas dim==2
        endif
    end subroutine initial_step

    subroutine initialize1D(data, celle)

        !In
        type(DataType), intent(in) :: data

        !Out
        type(StructCelleType), intent(inout) :: celle

        !Local
        integer :: i,j,k,l
        real(pr) :: init_radius(data%N_r)
        real(pr) :: init_velocity(data%N_vx)
        real(pr) :: init_temperature(data%N_T)
        real(pr) :: init_mass(data%N_m)

        celle%code = data%cas

        do i=1,data%N_r
            init_radius(i) = acceptation_rejet() !Radius
            !init_radius(i) = 0.0
        enddo

        do j=1,data%N_vx
            init_velocity(j) = 0.0 !Velocity
        enddo

        
        do k=1,data%N_m
            init_mass = 4.0_pr/3.0_pr*pi*data%rho_p*init_radius !Mass
            !init_mass(k) = 0.0
        enddo
        

        do l=1,data%N_T
            init_temperature(l) = data%T_p_0 !Temperature
        enddo

        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k=1,data%N_m !mass
                    do l=1,data%N_T !Temperature

                        celle%SOL(:,j,k,l) = init_radius
                        celle%SOL(i,:,k,l) = init_velocity
                        celle%SOL(i,j,:,l) = init_mass
                        celle%SOL(i,j,k,:) = init_temperature

                    enddo
                enddo
            enddo
        enddo

    end subroutine initialize1D
    
end module init_mod