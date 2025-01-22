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
        integer  :: i,j!,l
        real(pr) :: R, Fp, Fm, Ip, Im, tau, T, som
        real(pr) :: dr, dm

        dt = 10._pr

        ! -- dt computing
        do i=1,data%N_r !radius
            dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
            dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))
            do j=1,data%N_T !Temperature

                R = abs(R_function(data, mesh%r_tab(i)+dr/2._pr, &
                mesh%u_bar(i,j), mesh%m_tab(i)+dm/2._pr, &
                mesh%m_sel(i), mesh%T_tab(j)+data%dT/2._pr))

                T = abs(T_function(data, mesh%r_tab(i)+dr/2._pr, &
                mesh%u_bar(i,j), mesh%m_tab(i)+dm/2._pr, &
                mesh%m_sel(i), mesh%T_tab(j)+data%dT/2._pr))

                som = min(dr/R, data%dT/T)

                if (dt > som) dt = data%cfl*som
            enddo
        enddo


        SELECT CASE(data%time_scheme_key)
        CASE(1)

            !! > TRANSPORT < !!

            !!celle loop
            do i=1,data%N_r !radius

                dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                dm = abs(mesh%m_tab(i+1) - mesh%m_tab(i))

                do j=1,data%N_T !Temperature

                    R = abs(R_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%u_bar(i,j), mesh%m_tab(i)+dm/2._pr, &
                    mesh%m_sel(i), mesh%T_tab(j)+data%dT/2._pr))

                    T = abs(T_function(data, mesh%r_tab(i)+dr/2._pr, &
                    mesh%u_bar(i,j), mesh%m_tab(i)+dm/2._pr, &
                    mesh%m_sel(i), mesh%T_tab(j)+data%dT/2._pr))

                    !! -- For n_bar
                    !! - radius
                    if(i==1) then 
                        Fp = upwind_flux(mesh%n_bar(i,j),mesh%n_bar(i+1,j),R)
                        Fm = upwind_flux(0.0_pr,mesh%n_bar(i,j),R)
                    else if(i==data%N_r) then
                        Fp = upwind_flux(mesh%n_bar(i,j),0.0_pr,R)
                        Fm = upwind_flux(mesh%n_bar(i-1,j),mesh%n_bar(i,j),R)
                    else 
                        Fp = upwind_flux(mesh%n_bar(i,j),mesh%n_bar(i+1,j),R)
                        Fm = upwind_flux(mesh%n_bar(i-1,j),mesh%n_bar(i,j),R)
                    end if


                    mesh%n_bar(i,j) = mesh%n_bar(i,j) - dt*(Fp-Fm)/dr

                    !! - Temperature
                    if(j==1) then 
                        Ip = upwind_flux(mesh%n_bar(i,j),mesh%n_bar(i,j+1),R)
                        Im = upwind_flux(0.0_pr,mesh%n_bar(i,j),R)
                    else if(j==data%N_T) then
                        Ip = upwind_flux(mesh%n_bar(i,j),0.0_pr,R)
                        Im = upwind_flux(mesh%n_bar(i,j-1),mesh%n_bar(i,j),R)
                    else 
                        Ip = upwind_flux(mesh%n_bar(i,j),mesh%n_bar(i,j+1),R)
                        Im = upwind_flux(mesh%n_bar(i,j-1),mesh%n_bar(i,j),R)
                    end if

                    mesh%n_bar(i,j) = mesh%n_bar(i,j) - dt*(Ip-Im)/data%dT

                    !! -- For ubar
                    !! - radius
                    if(i==1) then 
                        Fp = upwind_flux(mesh%u_bar(i,j),mesh%u_bar(i+1,j),R)
                        Fm = upwind_flux(0.0_pr,mesh%u_bar(i,j),R)
                    else if(i==data%N_r) then
                        Fp = upwind_flux(mesh%u_bar(i,j),0.0_pr,R)
                        Fm = upwind_flux(mesh%u_bar(i-1,j),mesh%u_bar(i,j),R)
                    else 
                        Fp = upwind_flux(mesh%u_bar(i,j),mesh%u_bar(i+1,j),R)
                        Fm = upwind_flux(mesh%u_bar(i-1,j),mesh%u_bar(i,j),R)
                    end if

                    mesh%u_bar(i,j) = mesh%u_bar(i,j) - dt*(Fp-Fm)/dr

                    !! - Temperature
                    if(j==1) then 
                        Ip = upwind_flux(mesh%u_bar(i,j),mesh%u_bar(i,j+1),R)
                        Im = upwind_flux(0.0_pr,mesh%u_bar(i,j),R)
                    else if(j==data%N_T) then
                        Ip = upwind_flux(mesh%u_bar(i,j),0.0_pr,R)
                        Im = upwind_flux(mesh%u_bar(i,j-1),mesh%u_bar(i,j),R)
                    else 
                        Ip = upwind_flux(mesh%u_bar(i,j),mesh%u_bar(i,j+1),R)
                        Im = upwind_flux(mesh%u_bar(i,j-1),mesh%u_bar(i,j),R)
                    end if

                    mesh%u_bar(i,j) = mesh%u_bar(i,j) - dt*(Ip-Im)/data%dT
                enddo
            enddo

            !! > RELAXATION < !!
            
            do i=1,data%N_r !radius
                tau = tau_p(data, mesh%r_tab(i), mesh%m_tab(i))
                do j=1,data%N_T
                    mesh%u_bar(i,j) = (1 - exp(-dt/tau))*data%U_air + exp(-dt/tau)*mesh%u_bar(i,j)
                enddo
            enddo
        CASE(2)
            print*,"No other time scheme have been implemented for the moment..."

        END SELECT


    end subroutine
    
end module time_mod

