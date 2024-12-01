module flux_mod

    use structured_mesh_mod
    use functions_mod

    implicit none

    public :: upwind_flux1D
    
contains

    function upwind_flux1D(data, celle) result(flux)
        
        !In
        type(DataType), intent(in)  :: data
        type(StructCelleType), intent(in) :: celle

        !Out
        real(pr), dimension(data%N_r, data%N_vx, data%N_m, data%N_T) :: flux

        !Local
        integer :: i, j, k, l
        real(pr) :: r_mean, vx_mean, m_mean, T_mean, m_sel
        real(pr) :: R_plus, M_plus, T_plus, V_plus
        real(pr) :: R_moins, M_moins, T_moins, V_moins
        real(pr) :: Fiplus1demi, Fimoins1demi


        do i=2,data%N_r-1 !radius
            do j=2,data%N_vx-1 !velocity
                do k=2,data%N_m-1 !mass
                    do l=2,data%N_T-1 !Temperature
                        
                        r_mean = SUM(celle%SOL(:,j,k,l)) / data%N_r
                        vx_mean = SUM(celle%SOL(i,:,k,l)) / data%N_vx
                        m_mean = SUM(celle%SOL(i,j,:,l)) / data%N_m
                        T_mean = SUM(celle%SOL(i,j,k,:)) / data%N_T

                        m_sel = (4.0/3.0)*pi*r_mean**3*data%rho_p*data%Salinity_p/1000.0

                        V_plus = MAX(F_function(data, r_mean, vx_mean, m_mean)/m_mean, 0.0)
                        R_plus = MAX(R_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)
                        M_plus = MAX(M_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)
                        T_plus = MAX(T_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)
                
                        V_moins = MIN(F_function(data, r_mean, vx_mean, m_mean)/m_mean, 0.0)
                        R_moins = MIN(R_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)
                        M_moins = MIN(M_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)
                        T_moins = MIN(T_coeff(data, r_mean, vx_mean, m_mean, m_sel, T_mean), 0.0)

                        Fiplus1demi = R_plus*(celle%SOL(i,j,k,l) - celle%SOL(i-1,j,k,l))/data%dr +&
                                    V_plus*(celle%SOL(i,j,k,l) - celle%SOL(i,j-1,k,l))/data%dvx + &
                                    M_plus*(celle%SOL(i,j,k,l) - celle%SOL(i,j,k-1,l))/data%dm + &
                                    T_plus*(celle%SOL(i,j,k,l) - celle%SOL(i,j,k,l-1))/data%dT
                                    

                        Fimoins1demi = R_moins*(celle%SOL(i+1,j,k,l) - celle%SOL(i,j,k,l))/data%dr + &
                                    V_moins*(celle%SOL(i,j+1,k,l) - celle%SOL(i,j,k,l))/data%dvx + &
                                    M_moins*(celle%SOL(i,j,k+1,l) - celle%SOL(i,j,k,l))/data%dm + &
                                    T_moins*(celle%SOL(i,j,k,l+1) - celle%SOL(i,j,k,l))/data%dT

                        flux(i,j,k,l) = Fiplus1demi + Fimoins1demi


                    end do
                end do
            end do
        end do

        !CL
        flux(1,:,:,:) = flux(2,:,:,:)
        flux(:,1,:,:) = flux(:,2,:,:)
        flux(:,:,1,:) = flux(:,:,2,:)
        flux(:,:,:,1) = flux(:,:,:,2)

        flux(data%N_r,:,:,:) = flux(data%N_r-1,:,:,:)
        flux(:,data%N_vx,:,:) = flux(:,data%N_vx-1,:,:)
        flux(:,:,data%N_m,:) = flux(:,:,data%N_m-1,:)
        flux(:,:,:,data%N_T) = flux(:,:,:,data%N_T-1)



    end function upwind_flux1D
    
end module flux_mod