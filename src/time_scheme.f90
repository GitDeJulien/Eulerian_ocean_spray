subroutine advance (data, mesh, dt)

    use flux_mod
    implicit none
    

    !In
    type(DataType), intent(in) :: data

    !Inout
    type(MeshType), intent(inout) :: mesh

    !Out
    real(pr) :: dt

    !Local
    integer :: i,j,l
    real(pr) :: R, M, T, V, m_sel, som

    do i=1,data%N_r !radius
        do j=1,data%N_vx !velocity
            do l=1,data%N_T !Temperature
            

                m_sel = (4.0/3.0)*pi*mesh%r_tab(i)**3*data%rho_p*data%Salinity_p/1000.0

                V = abs(F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(i))/mesh%m_tab(i))
                R = abs(R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(i), m_sel, mesh%T_tab(l)))
                M = abs(M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(i), m_sel, mesh%T_tab(l)))
                T = abs(T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(i), m_sel, mesh%T_tab(l)))

                som = data%cfl * (data%dr/R + data%dvx/V + data%dm/M + data%dT/T) !j'en suis pas sur :*(

                print*, som

                if (dt > som) dt = som

            enddo
        enddo
    enddo


    SELECT CASE(data%time_scheme_key)
    CASE(1)
        mesh%SOL = mesh%SOL - dt*vector_flux(data, mesh)
        ! do i=1,data%N_r !radius
        !     do j=1,data%N_vx !velocity
        !         do l=1,data%N_T !Temperature
                
        !             mesh%SOL(i,j,i,l) = mesh%SOL(i,j,i,l) - 
        !               dt*(fonction_flux_de_maïlys(data, ud, um) + fonction_flux_de_maïlys(data, um, ug))
                
        !         enddo
        !     enddo
        ! enddo
    CASE(2)
        print*,"No other time scheme have been implemented for the moment..."

    END SELECT





end subroutine