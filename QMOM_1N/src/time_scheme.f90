module time_mod

    implicit none
    
contains

    subroutine advance (data, mesh, dt)

        use flux_mod
        implicit none
        

        !In
        type(DataType), intent(in) :: data

        !Inout
        type(MeshType), intent(inout) :: mesh

        !Out
        real(pr), intent(out) :: dt

        !Local
        integer  :: i!,j,l
        real(pr) :: R, Fp, Fm, tau!, som, M, T, V, som
        real(pr) :: dr, dm

        dt = 10._pr

        ! do i=1,data%N_r !radius
        !     do j=1,data%N_vx !velocity
        !         do l=1,data%N_T !Temperature

        !             dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
        !             dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
                
        !             V = abs(mesh%V_coeff(i,j,i,l))
        !             R = abs(mesh%R_coeff(i,j,i,l))
        !             M = abs(mesh%M_coeff(i,j,i,l))
        !             T = abs(mesh%T_coeff(i,j,i,l))

        !             som = min(data%dvx/V, dr/R, dm/M, data%dT/T)
        !             print*, "T=",T

        !             if (dt > som) dt = som

        !         enddo
        !     enddo
        ! enddo

        ! -- dt computing
        do i=1,data%N_r !radius
            dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
            dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))

            R = abs(R_function(data, mesh%r_tab(i)+dr/2._pr, &
            mesh%u_bar(i), mesh%m_tab(i)+dm/2._pr, mesh%m_sel(i), data%T_air))

            if (dt > dr/R) dt = data%cfl*dr/R
        enddo


        SELECT CASE(data%time_scheme_key)
        CASE(1)
            !! > TRANSPORT < !!

            !!celle loop
            do i=1,data%N_r !radius

                dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))

                R = R_function(data, mesh%r_tab(i)+dr/2._pr, &
                mesh%u_bar(i), mesh%m_tab(i)+dm/2._pr, mesh%m_sel(i), data%T_air)

                !! -- For n_bar
                if(i==1) then 
                    Fp = upwind_flux(mesh%n_bar(i),mesh%n_bar(i+1),R)
                    Fm = upwind_flux(0.0_pr,mesh%n_bar(i),R)
                else if(i==data%N_r) then
                    Fp = upwind_flux(mesh%n_bar(i),0.0_pr,R)
                    Fm = upwind_flux(mesh%n_bar(i-1),mesh%n_bar(i),R)
                else 
                    Fp = upwind_flux(mesh%n_bar(i),mesh%n_bar(i+1),R)
                    Fm = upwind_flux(mesh%n_bar(i-1),mesh%n_bar(i),R)
                end if


                mesh%n_bar(i) = mesh%n_bar(i) - dt*(Fp-Fm)/dr

                !! -- For ubar
                if(i==1) then 
                    Fp = upwind_flux(mesh%u_bar(i),mesh%u_bar(i+1),R)
                    Fm = upwind_flux(0.0_pr,mesh%u_bar(i),R)
                else if(i==data%N_r) then
                    Fp = upwind_flux(mesh%u_bar(i),0.0_pr,R)
                    Fm = upwind_flux(mesh%u_bar(i-1),mesh%u_bar(i),R)
                else 
                    Fp = upwind_flux(mesh%u_bar(i),mesh%u_bar(i+1),R)
                    Fm = upwind_flux(mesh%u_bar(i-1),mesh%u_bar(i),R)
                end if

                mesh%u_bar(i) = mesh%u_bar(i) - dt*(Fp-Fm)/dr
            enddo

            !! > RELAXATION < !!
            do i=1,data%N_r !radius
                tau = tau_p(data, mesh%r_tab(i), mesh%m_tab(i))
                mesh%u_bar(i) = (1 - exp(-dt/tau))*data%U_air + exp(-dt/tau)*mesh%u_bar(i)
            enddo
        CASE(2)
            print*,"No other time scheme have been implemented for the moment..."

        END SELECT


    end subroutine
    
end module time_mod

