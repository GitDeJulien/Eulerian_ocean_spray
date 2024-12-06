module functions_mod

    use precision
    use data_mod

    implicit none

    private :: f_v, Dv_star_function, ka_star_function
    public  :: Re_p, tau_p, tau_t, F_function, M_function, R_function, T_function, b_function

contains

function Re_p(data, r_p, v_p) result(res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: res

    !Local
    real(pr) :: v_s

    v_s = v_p - data%U_air
    res = 2.0*r_p*abs(v_s)/data%nu_air

end function Re_p

function tau_p(data, r_p, m_p) result(res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p
    real(pr), intent(in)         :: m_p 

    !Out
    real(pr)     :: res

    res = (2.0*r_p**2*(m_p/(4./3.*pi*r_p**3)))/(9.0*data%mu_air);

end function tau_p

function F_function(data, r_p, v_p, m_p) result(res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p, m_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: res

    !Local
    real(pr) :: tau_D
    real(pr) :: fct

    tau_D = tau_p(data, r_p, m_p)
    fct = (m_p/tau_D)*(data%U_air - v_p)
    res = fct


end function F_function

function f_v(data, r_p, v_p) result(res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: res

    res = 1.0 + sqrt(Re_p(data, r_p,v_p))/4.0

end function f_v



function Dv_star_function(data, r_p, v_p) result(Dv_star)

    !> Function computing the Dv* expression
    !> Dv* is the modified diffusivity for water vapour

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: Dv_star

    !Local
    real(pr) :: b, c
    real(pr) :: a


    a = f_v(data, r_p,v_p)*data%D_v
    b = r_p/(r_p + data%Delta_v)
    c = (data%D_v/(r_p*data%alpha_c))*sqrt((2.0*pi*data%M_w)/(data%R_g*data%T_air))

    Dv_star = a/(b+c)

end function Dv_star_function

function ka_star_function(data, r_p, v_p) result(ka_star)
    !> Function computing the ka* expression

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)       :: r_p
    real(pr), intent(in)       :: v_p 

    !Out
    real(pr) :: ka_star

    !Local
    real(pr) :: b, c
    real(pr) :: a


    a = f_v(data, r_p,v_p)*data%k_a
    b = r_p/(r_p + data%Delta_T)
    c = (data%k_a/(r_p*data%alpha_T*data%rho_air_sec*data%c_p_air_sec))&
    *sqrt((2.0*pi*data%M_a)/(data%R_g*data%T_air))

    ka_star = a/(b+c)

end function ka_star_function


function M_function(data, r_p, v_p, m_p, m_sel, T_p) result(M_res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)         :: r_p, m_p, m_sel, T_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: M_res

    !Local
    real(pr) :: a, b

    a = (4*pi*r_p*Dv_star_function(data, r_p, v_p)*data%M_w*data%pv_sat_T_air)/&
    (data%R_g*data%T_air)
    b = data%Q_RH - (data%T_air/T_p)*EXP(((data%L_v*data%M_w)/data%R_g)&
    *((1.0_pr/data%T_air) - 1.0_pr/T_p) + (2.0_pr*data%M_w*data%Gamma_p)/&
    (data%R_g*data%rho_w*r_p*T_p) - &
    (data%I*data%Phi_s*m_sel*(data%M_w/data%M_s))/(m_p-m_sel))

    M_res = a*b

end function M_function

function R_function(data, r_p, v_p, m_p, m_sel, T_p) result(R_res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)       :: r_p, m_p, m_sel, T_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: R_res

    R_res = M_function(data,r_p, v_p, m_p, m_sel, T_p) / &
    (4.0_pr*pi*r_p**2*data%rho_w)

end function R_function

function T_function(data, r_p, v_p, m_p, m_sel, T_p) result(T_res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)       :: r_p, m_p, m_sel, T_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: T_res

    !Local
    real(pr) :: a

    a = (4.0_pr*pi*ka_star_function(data, r_p, v_p)*&
    (data%T_air - T_p))/(m_p*data%c_p_s)

    T_res = a + (data%L_v/(m_p*data%c_p_s))*&
    M_function(data,r_p, v_p, m_p, m_sel, T_p)

end function T_function

function b_function(data, r_p, v_p, m_p, m_sel, T_p) result(b_res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)       :: r_p, m_p, m_sel, T_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: b_res

    b_res = data%T_air + (data%L_v*M_function(data, r_p, v_p, m_p, m_sel, T_p)/(4*pi*r_p*ka_star_function(data, r_p, v_p)))

end function b_function

function tau_t(data, r_p, v_p, m_p) result(tau_t_res)

    !In
    type(DataType), intent(in) :: data
    real(pr), intent(in)       :: r_p, m_p
    real(pr), intent(in) :: v_p 

    !Out
    real(pr) :: tau_t_res

    tau_t_res = m_p*data%c_p_s/(4*pi*r_p*ka_star_function(data, r_p, v_p))

end function tau_t



end module