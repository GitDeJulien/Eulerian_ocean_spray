module init_mod

    use functions_mod
    use structured_mesh_mod
    use sample_mod

    implicit none

    private :: find_best_loc
    public :: initialize_sol!, initialize_coeff
    
contains


    subroutine initialize_sol(data, mesh)

        !In
        type(DataType), intent(in) :: data

        !Out
        type(MeshType), intent(inout) :: mesh

        !Local
        integer  :: i,j,k
        integer  :: v_loc
        real(pr) :: dm

        v_loc = find_best_loc(mesh%vx_tab, 0.0_pr)

        allocate(mesh%n_bar(data%N_r, data%N_T))
        allocate(mesh%u_bar(data%N_r, data%N_T))

        mesh%n_bar = 0.d0
        mesh%u_bar = 0.d0

        do i=1,data%N_r !radius
            dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
            do j=1,data%N_vx !velocity
                do k=1,data%N_T
                    if ( j == v_loc) then
                        mesh%n_bar(i,k) = dCd_r(dFdr(mesh%r_tab(i)), Vdp(vp4(mesh%r_tab(i))))/&
                        (dm*data%dT*data%dvx)
                    endif
                enddo
            enddo
        enddo

        !print*, mesh%n_bar

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
            if (res > eps) then
                loc = i
                res = eps
            end if
        enddo

    end function find_best_loc
    
end module init_mod