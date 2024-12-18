module init_mod

    use functions_mod
    use structured_mesh_mod
    use sample_mod

    implicit none

    private :: find_best_loc
    public :: initialize_sol, initialize_coeff
    
contains


    subroutine initialize_sol(data, mesh)

        !In
        type(DataType), intent(in) :: data

        !Out
        type(MeshType), intent(inout) :: mesh

        !Local
        integer  :: i,j,k
        integer  :: T_loc, v_loc
        real(pr) :: dm

        allocate(mesh%SOL(data%N_r, data%N_vx, data%N_m, data%N_T))

        v_loc = find_best_loc(mesh%vx_tab, 0.0_pr)

        T_loc = find_best_loc(mesh%T_tab, data%T_p_0)

        mesh%SOL = 0.0

        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k=1,data%N_T !Temperature
                    

                    if ( j == v_loc .and. k == T_loc) then
                        dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
                        mesh%SOL(i,j,i,k) = dCd_r(dFdr(mesh%r_tab(i)), Vdp(vp4(mesh%r_tab(i)))) /&
                        (dm*data%dT*data%dvx)
                    end if

                enddo
            enddo
        enddo

    end subroutine initialize_sol

    subroutine initialize_coeff(data, mesh)
        
        !In
        type(DataType), intent(in) :: data
        

        !InOut
        type(MeshType), intent(inout) :: mesh

        !Local
        integer  :: i,j,l
        real(pr) :: dr, dm

        allocate(mesh%m_sel(data%N_r))
        allocate(mesh%V_coeff(data%N_r, data%N_vx, data%N_m, data%N_T))
        allocate(mesh%R_coeff(data%N_r, data%N_vx, data%N_m, data%N_T))
        allocate(mesh%M_coeff(data%N_r, data%N_vx, data%N_m, data%N_T))
        allocate(mesh%T_coeff(data%N_r, data%N_vx, data%N_m, data%N_T))


        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do l=1,data%N_T !Temperature

                    dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                    dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
                
                    mesh%m_sel(i) = (4.0/3.0)*pi*(mesh%r_tab(i)+dr/2._pr)**3*data%rho_p*data%Salinity_p/1000.0
    
                    mesh%V_coeff(i,j,i,l) = F_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%vx_tab(j), mesh%m_tab(i)+dm/2._pr)/mesh%m_tab(i)

                    mesh%R_coeff(i,j,i,l) = R_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%vx_tab(j), mesh%m_tab(i)+dm/2._pr, mesh%m_sel(i), mesh%T_tab(l))

                    mesh%M_coeff(i,j,i,l) = M_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%vx_tab(j), mesh%m_tab(i)+dm/2._pr, mesh%m_sel(i), mesh%T_tab(l))
                    
                    mesh%T_coeff(i,j,i,l) = T_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%vx_tab(j), mesh%m_tab(i)+dm/2._pr, mesh%m_sel(i), mesh%T_tab(l))
    
                enddo
            enddo
        enddo


    end subroutine initialize_coeff

    function find_best_loc(tab, val) result(loc)

        !In
        real(pr), dimension(:), intent(in) :: tab
        real(pr), intent(in)               :: val

        !Out
        integer  :: loc
        
        !Local
        real(pr) :: res
        integer :: N, i
        real(pr) :: eps

        N = size(tab)
        res = 10.0

        do i=1,N
            
            eps = abs(tab(i)-val)
            !print*, "eps=", eps
            if (res > eps) then
                loc = i
                res = eps
            end if
        enddo

    end function find_best_loc
    
end module init_mod