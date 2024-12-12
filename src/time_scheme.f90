module time_mod

    use flux_mod
    implicit none

    contains

    subroutine advance (data, mesh, dt)

        !In
        type(DataType), intent(in) :: data

        !Inout
        type(MeshType), intent(inout) :: mesh

        !Out
        real(pr) :: dt

        !Local
        integer  :: i,j,k,l
        real(pr) :: R, M, T, V, som
        real(pr) :: dr, dm

        dt = 10._pr

        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k=1,data%N_m
                    do l=1,data%N_T !Temperature

                        dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                        dm = abs(mesh%m_tab(k+1) - mesh%m_tab(k))
                    
                        V = abs(mesh%V_coeff(i,j,k,l))
                        R = abs(mesh%R_coeff(i,j,k,l))
                        M = abs(mesh%M_coeff(i,j,k,l))
                        T = abs(mesh%T_coeff(i,j,k,l))

                        som = min(data%dvx/V, dr/R, dm/M, data%dT/T)
                        ! print*, "T=",T

                        if (dt > som) dt = som

                    enddo
                enddo
            enddo
        enddo

        SELECT CASE(data%time_scheme_key)
        CASE(1)
            mesh%SOL = mesh%SOL - dt*vector_flux(data, mesh)
        CASE(2)
            print*,"No other time scheme have been implemented for the moment..."

        END SELECT





    end subroutine

end module