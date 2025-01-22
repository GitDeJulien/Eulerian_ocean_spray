module flux_mod

    use structured_mesh_mod
    use functions_mod

    implicit none

    public :: upwind_flux
    
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
        integer  :: i, j, k, l
        real(pr) :: V, R, M, T, Fp, Fm, Gp, Gm, Hp, Hm, Ip, Im
        real(pr) :: dr, dm


        do i=1,data%N_r !radius
            do j=1,data%N_vx !velocity
                do k = 1,data%N_m !masse
                    do l=1,data%N_T !Temperature

                        dr = abs(mesh%r_tab(i+1) - mesh%r_tab(i))
                        dm = abs(mesh%m_tab(k+1) - mesh%m_tab(k))
                    
                        V = mesh%V_coeff(i,j,k,l)
                        R = mesh%R_coeff(i,j,k,l)
                        M = mesh%M_coeff(i,j,k,l)
                        T = mesh%T_coeff(i,j,k,l)

                        !radius
                        if(i==1) then 
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i+1,j,k,l),R)
                            Fm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),R)
                        else if(i==data%N_r) then
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,R)
                            Fm = upwind_flux(mesh%SOL(i-1,j,k,l),mesh%SOL(i,j,k,l),R)
                        else 
                            Fp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i+1,j,k,l),R)
                            Fm = upwind_flux(mesh%SOL(i-1,j,k,l),mesh%SOL(i,j,k,l),R)
                        end if

                        !velocity
                        if(j==1) then
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j+1,k,l),V)
                            Gm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),V)
                        else if(j==data%N_vx) then
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,V)
                            Gm = upwind_flux(mesh%SOL(i,j-1,k,l),mesh%SOL(i,j,k,l),V)
                        else 
                            Gp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j+1,k,l),V)
                            Gm = upwind_flux(mesh%SOL(i,j-1,k,l),mesh%SOL(i,j,k,l),V)
                        end if

                        !masse
                        if(k==1) then
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k+1,l),M)
                            Hm = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),M)
                        else if(k==data%N_m) then
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,M)
                            Hm = upwind_flux(mesh%SOL(i,j,k-1,l),mesh%SOL(i,j,k,l),M)
                        else 
                            Hp = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k+1,l),M)
                            Hm = upwind_flux(mesh%SOL(i,j,k-1,l),mesh%SOL(i,j,k,l),M)
                        end if

                        !Temperature
                        if(l==1) then
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k,l+1),T)
                            Im = upwind_flux(0.0_pr,mesh%SOL(i,j,k,l),T)
                        else if (l==data%N_T) then
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),0.0_pr,T)
                            Im = upwind_flux(mesh%SOL(i,j,k,l-1),mesh%SOL(i,j,k,l),T)
                        else
                            Ip = upwind_flux(mesh%SOL(i,j,k,l),mesh%SOL(i,j,k,l+1),T)
                            Im = upwind_flux(mesh%SOL(i,j,k,l-1),mesh%SOL(i,j,k,l),T)
                        end if

                        flux(i,j,k,l) = (Fp-Fm)/dr + (Gp-Gm)/data%dvx + (Hp-Hm)/dm + (Ip-Im)/data%dT 
                    end do
                end do
            end do
        end do

    end function vector_flux
    
end module flux_mod