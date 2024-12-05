module flux_mod

    use data_mod
    use structured_mesh_mod
    use functions_mod

    implicit none

    public :: upwind_flux
    public :: vector_flux
    
contains

    function upwind_flux(solg,sold,vel) result(Flux)

        real(pr), intent(in) :: solg, sold, vel
        real(pr)             :: Flux

        real(pr)             :: vel_p, vel_m

        vel_p = MAX(vel, 0.0)
        vel_m = MIN(vel, 0.0)

        Flux = vel_p*solg + vel_m*sold


    end function upwind_flux

    function vector_flux(data, mesh) result(flux)
        
        !In
        type(DataType), intent(in)  :: data
        type(MeshType), intent(in) :: mesh

        !Out
        real(pr), dimension(data%N_r, data%N_vx, data%N_m, data%N_T) :: flux

        !Local
        integer :: i, j, k, l
        real(pr) :: m_sel
        real(pr) :: Fp, Fm, Gp, Gm, Hp, Hm, Ip, Im


        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k = 1,data%N_m !masse
                    do l=1,data%N_T !Temperature
                    
                        m_sel = (4.0/3.0)*pi*mesh%r_tab(i)**3*data%rho_p*data%Salinity_p/1000.0

                        !radius
                        if(i==1) then 
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i+1,j,k,l),&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Fm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else if(i==data%N_r) then
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Fm = upwind_flux(mesh%SOL(i-1,j,k,l),mesh%SOL(i,j,k,l),&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else 
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i+1,j,k,l),&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Fm = upwind_flux(mesh%SOL(i-1,j,k,l),mesh%SOL(i,j,k,l),&
                                R_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        end if

                        !velocity
                        if(j==1) then
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j+1,k,l),&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                            Gm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                        else if(j==data%N_vx) then
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                            Gm = upwind_flux(mesh%SOL(i,j-1,k,l),mesh%SOL(i,j,k,l),&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                        else 
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j+1,k,l),&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                            Gm = upwind_flux(mesh%SOL(i,j-1,k,l),mesh%SOL(i,j,k,l),&
                                F_function(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k))/mesh%m_tab(k))
                        end if

                        !masse
                        if(k==1) then
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k+1,l),&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Hm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else if(k==data%N_m) then
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Hm = upwind_flux(mesh%SOL(i,j,k-1,l),mesh%SOL(i,j,k,l),&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else 
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k+1,l),&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Hm = upwind_flux(mesh%SOL(i,j,k-1,l),mesh%SOL(i,j,k,l),&
                                M_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        end if

                        !Temperature
                        if(l==1) then
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k,l+1),&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Im = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else if (l==data%N_T) then
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Im = upwind_flux(mesh%SOL(i,j,k,l-1),mesh%SOL(i,j,k,l),&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        else
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k,l+1),&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                            Im = upwind_flux(mesh%SOL(i,j,k,l-1),mesh%SOL(i,j,k,l),&
                                T_coeff(data, mesh%r_tab(i), mesh%vx_tab(j), mesh%m_tab(k), m_sel, mesh%T_tab(l)))
                        end if

                        flux(i,j,i,l) = (Fp-Fm)/data%dr + (Gp-Gm)/data%dvx + (Hp-Hm)/data%dm + (Ip-Im)/data%dT 
                    end do
                end do
            end do
        end do

    end function vector_flux
    
end module flux_mod