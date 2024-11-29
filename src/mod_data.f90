module data_mod

    use precision
    use toml_parser
    implicit none

    public :: config_data

    type, public :: DataType

        ![USEFULL PATH]
        character(len=256) :: output_path

        ![PHYSICAL PARMAETERS]
        real(pr) :: U_air
        real(pr) :: T_air_celcius
        real(pr) :: T_p_0_celcius
        real(pr) :: Q_RH
        real(pr) :: M_w
        real(pr) :: M_s
        real(pr) :: M_a
        real(pr) :: R_g
        real(pr) :: p_0
        real(pr) :: Salinity_w
        real(pr) :: Salinity_p
        real(pr) :: Delta_v
        real(pr) :: alpha_c
        real(pr) :: Delta_T
        real(pr) :: alpha_T
        real(pr) :: I

        ![PHYSICAL PARMAETERS COMPUTED]
        real(pr) :: k_a
        real(pr) :: T_air
        real(pr) :: T_p_0
        real(pr) :: qs0
        real(pr) :: rs0
        real(pr) :: r0
        real(pr) :: q0
        real(pr) :: Tv0
        real(pr) :: rho_air
        real(pr) :: nu_air
        real(pr) :: mu_air
        real(pr) :: p_air
        real(pr) :: pv_sat_T_air
        real(pr) :: L_v
        real(pr) :: D_v
        real(pr) :: Gamma_p
        real(pr) :: Phi_s
        real(pr) :: rho_w
        real(pr) :: rho_p
        real(pr) :: ms
        real(pr) :: c_p_s
        real(pr) :: c_p_air_sec
        real(pr) :: rho_air_sec
        
        ![NUMERICAL PARMAETERS]
        integer  :: dim
        integer  :: cas
        real(pr) :: Lx
        integer  :: Nx
        integer  :: n_celle
        real(pr) :: dx

    end type DataType

contains

    subroutine config_data(data, filename)
        character(len=*), intent(in) :: filename
        type(DataType), intent(inout) :: data

        call parse_toml(filename, "output_path", data%output_path)
        call parse_toml(filename, "U_air", data%U_air)
        call parse_toml(filename, "T_air_celcius", data%T_air_celcius)
        call parse_toml(filename, "T_p_0_celcius", data%T_p_0_celcius)
        call parse_toml(filename, "Q_RH", data%Q_RH)
        call parse_toml(filename, "M_w", data%M_w)
        call parse_toml(filename, "M_s", data%M_s)
        call parse_toml(filename, "M_a", data%M_a)
        call parse_toml(filename, "R_g", data%R_g)
        call parse_toml(filename, "p_0", data%p_0)
        call parse_toml(filename, "Salinity_w", data%Salinity_w)
        call parse_toml(filename, "Salinity_p", data%Salinity_p)
        call parse_toml(filename, "Delta_v", data%Delta_v)
        call parse_toml(filename, "alpha_c", data%alpha_c)
        call parse_toml(filename, "Delta_T", data%Delta_T)
        call parse_toml(filename, "alpha_T", data%alpha_T)
        call parse_toml(filename, "I", data%I)

        call parse_toml(filename, "dim", data%dim)
        call parse_toml(filename, "cas", data%cas)
        call parse_toml(filename, "Lx", data%Lx)
        call parse_toml(filename, "Nx", data%Nx)

        data%k_a = 0.02411*(1+3.309e-3*(data%T_p_0_celcius) -1.441e-6*(data%T_p_0_celcius)**2)
        data%T_air = data%T_air_celcius + 273.15
        data%T_p_0 = data%T_p_0_celcius + 273.15
        data%qs0 = 0.62197*611.2/data%p_0*EXP((17.67*(data%T_p_0_celcius))/(data%T_p_0-29.66))
        data%rs0 = data%qs0/(1-data%qs0)
        data%r0 = 98.0/100.0*data%rs0
        data%q0 = data%r0/(1+data%r0)
        data%Tv0 = data%T_p_0*(1+0.6078*data%q0)
        data%rho_air = 1.2929*273.15/data%Tv0
        data%nu_air = 1.33e-5 +(0.0084*(data%Tv0-273.15))*1e-5
        data%mu_air = data%nu_air*data%rho_air
        data%p_air = data%p_0-data%rho_air*gravity*10.0
        data%pv_sat_T_air = 2.33e3
        data%L_v = (25.00895-0.02274*(data%T_air-273.15))*1e5
        data%D_v = 2.11e-5*(data%T_air-273.15)**1.94*(1013.25/(data%p_air/100.0))
        data%Gamma_p = (75.63 - 0.144*data%T_p_0_celcius + 0.221*data%Salinity_w)*1e-3
        data%Phi_s = 0.91154614+1.7317496707e-4*data%Salinity_p+4.7616058412e-6*data%Salinity_p**2-&
        & 9.2541509027e-9*data%Salinity_p**3+7.3475024678e-12*data%Salinity_p**4
        data%rho_w = (999.842594+6.793952e-2*data%T_p_0_celcius-9.095290e-3*data%T_p_0_celcius**2+&
        & 1.001685e-4*data%T_p_0_celcius**3-1.120083e-6*data%T_p_0_celcius**4+6.536332e-9*data%T_p_0_celcius**5)
        data%rho_p = data%rho_w + data%Salinity_p*(0.824493-4.0899e-3*data%T_p_0_celcius + &
        & 7.6438e-5*data%T_p_0_celcius**2-8.2467e-7*data%T_p_0_celcius**3+5.3875e-9*data%T_p_0_celcius**4)+&
        & data%Salinity_p**(3./2.)*(-5.72466e-3+1.0227e-4*data%T_p_0_celcius&
        & -1.6546e-6*data%T_p_0_celcius**2) +4.8314e-4*data%Salinity_p**2
        data%c_p_s = 4217.4 -3.720283*(data%T_air-273.15)+0.1412855*(data%T_air-273.15)**2&
        & -2.654387e-3*(data%T_air-273.15)**3 +2.093236e-5*(data%T_air-273.15)**4 +&
        & data%Salinity_p*(-7.6444+0.107276*(data%T_air-273.15)-1.3839e-3*data%T_p_0_celcius**2)+&
        & data%Salinity_p**(3./2)*(0.17709-4.0772e-3*(data%T_air-273.15)+5.3539e-5*data%T_p_0_celcius**2)
        data%c_p_air_sec = 1.9327e-10*data%T_p_0**4-7.9999e-7*data%T_p_0**3+&
        & 1.1407e-3*data%T_p_0**2-4.4890e-1*data%T_p_0+1.0575e+3
        data%rho_air_sec = 1.2929*273.13/data%T_air


        if (data%dim == 1) then
            data%dx = data%Lx / (data%Nx)
        end if

        if (data%dim == 1) then
            data%n_celle = data%Nx-1
        else
            print*, "Error: Dimension can't be greater than 1 for the moment"
            stop
        end if

    end subroutine config_data
    
end module data_mod