subroutine advance (data, mesh, dt)

    use data_mod
    use structured_mesh_mod
    use flux_mod
    implicit none
    

    !In
    type(DataType), intent(in) :: data

    !Inout
    type(MeshType), intent(inout) :: mesh

    !Out
    real(pr) :: dt

    !Local
    integer  :: i,j,l
    real(pr) :: R, M, T, V, som
    real(pr) :: dr, dm

    dt = 10._pr

    do i=1,data%N_r !radius
        do j=1,data%N_vx !velocity
            do l=1,data%N_T !Temperature

                dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
            
                V = abs(mesh%V_coeff(i,j,i,l))
                R = abs(mesh%R_coeff(i,j,i,l))
                M = abs(mesh%M_coeff(i,j,i,l))
                T = abs(mesh%T_coeff(i,j,i,l))

                som = min(data%dvx/V, dr/R, dm/M, data%dT/T)
                print*, "T=",T

                if (dt > som) dt = som

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