module init_mod

    use structured_mesh_mod
    use sample_mod

    implicit none

    private :: find_best_loc
    public :: initialize_sol
    
contains


    subroutine initialize_sol(data, mesh)

        !In
        type(DataType), intent(in) :: data

        !Out
        type(MeshType), intent(inout) :: mesh

        !Local
        integer :: i,j,k
        integer :: T_loc, v_loc

        allocate(mesh%SOL(data%N_r, data%N_vx, data%N_m, data%N_T))

        v_loc = find_best_loc(mesh%vx_tab, 0.0_pr)

        T_loc = find_best_loc(mesh%T_tab, data%T_p_0)

        mesh%SOL = 0.0

        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k=1,data%N_T !Temperature

                    if ( j == v_loc .and. k == T_loc) then
                        mesh%SOL(i,j,i,k) = dCd_r(dFdr(mesh%r_tab(i)), Vdp(vp4(mesh%r_tab(i))))
                    end if

                enddo
            enddo
        enddo

    end subroutine initialize_sol

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