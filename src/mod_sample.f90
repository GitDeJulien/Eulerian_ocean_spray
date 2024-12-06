module sample_mod

    use precision
    implicit none

    !private :: normalized_N_r!, acceptation_rejet
    public :: dCd_r, dFdr, Vdp!, normalized_N_r
    
contains
    
    function dFdr(R) result(dFdr_F)

        !In
        real(pr), intent(in) :: R

        !Out
        real(pr) :: dFdr_F

        !Local
        real(pr) :: r_80, Uwind10, a
        real(pr), parameter :: NAN = -huge(0.0)  ! Define NAN as a large negative number

        Uwind10 = 15.0
        r_80 = 0.518 * (R * 1.0e6)**0.976

        if (r_80 > 0.9 .and. r_80 < 15) then
            dFdr_F = 75.7672 / 2.0 * 3.8e-6 * Uwind10**3.4 * r_80**(-0.024) * &
                    10**(4.405 - 2.646 * log10(r_80) - 3.156 * log10(r_80)**2 + &
                        8.902 * log10(r_80)**3 - 4.482 * log10(r_80)**4)
        else if (r_80 >= 15 .and. r_80 <= 37.5) then
            dFdr_F = 75.7672 / 2.0 * 3.8e-6 * Uwind10**3.4 * r_80**(-0.024) * &
                    1.02e4 * r_80**(-1)
        else if (r_80 >= 37.5 .and. r_80 <= 100) then
            dFdr_F = 75.7672 / 2.0 * 3.8e-6 * Uwind10**3.4 * r_80**(-0.024) * &
                    6.95e6 * r_80**(-2.8)
        else if (r_80 > 100) then
            dFdr_F = 75.7672 / 2.0 * 3.8e-6 * Uwind10**3.4 * r_80**(-0.024) * &
                    1.75e17 * r_80**(-8)
        else
            dFdr_F = NAN
        end if

        if (r_80 > 0.5 .and. r_80 < 0.9) then
            a = 34623
            dFdr_F = a * r_80**(-3)
        end if

        dFdr_F = dFdr_F * 1.0e6

    end function dFdr


    function vp4(R) result(Vt_p)

        !In
        real(pr), intent(in) :: R

        !Out
        real(pr) :: Vt_p

        !Local

        real(pr), parameter :: lambda = 68e-9
        real(pr), parameter :: rhop = 1024.0
        real(pr), parameter :: rhoa = 1.19448
        real(pr), parameter :: nua = 1.51905e-05
        real(pr) :: Cd_p, RE, CU, vp2, f
        integer :: iter

        CU = 1 + lambda / R * (1.257 + 0.4 * exp(-0.55 * 2 * R / lambda))

        vp2 = 0.01 * gravity * (R * 2)**2 * (rhop - rhoa) / (18 * nua * rhoa)

        do iter = 1, 200
            RE = 2 * R * vp2 / nua

            if (RE < 500) then
                Cd_p = 24.0 / RE * (1 + 0.15 * RE**0.687 + &
                                    0.175 * (1 + 4.25e4 * RE**(-1.16))**(-1))
            else
                Cd_p = 0.85 - 9.76e-4 * RE + 1.091e-6 * RE**2 - &
                    6.84e-10 * RE**3 + 2.72e-13 * RE**4 - &
                    6.68e-17 * RE**5 + 9.83e-21 * RE**6 - &
                    7.96e-25 * RE**7 + 2.73e-29 * RE**8
            end if

            Cd_p = Cd_p / CU

            f = Cd_p * RE / 24.0
            vp2 = 0.5 * (gravity * (R * 2)**2 * (rhop - rhoa) / (18 * nua * rhoa * f)) + 0.5 * vp2

            if (2 * R > 1e-2) then
                vp2 = 9.36
            end if
        end do

        Vt_p = vp2

    end function vp4

    function Vdp(Vt_p)

        real(pr) :: Vdp
        real(pr), intent(in) :: Vt_p

        !Local
        real(pr) :: Uwind10, res

        Uwind10 = 15.0

        res = Vt_p/(1-EXP(-Vt_p/1e-3/Uwind10))

        Vdp = res

    end function

    function dCd_r(dFdr_var, Vdp_r) result(dCdr)

        !In
        real(pr), intent(in) :: dFdr_var, Vdp_r

        !Out
        real(pr) :: dCdr

        dCdr = dFdr_var/Vdp_r

    end function dCd_r

    ! function normalized_N_r(dCdr) result(Nr_normed)

    !     !In
    !     real(pr), intent(in) :: dCdr

    !     !Out
    !     real(pr) :: Nr_normed

    !     !Local
    !     real(pr), parameter :: K = 8.55422e+08

    !     Nr_normed = dCdr/K

    ! end function normalized_N_r

    ! function acceptation_rejet(seed) result(y)
    
    !     integer, intent(in), optional :: seed
    
    !     real(pr) :: max, a, b, c, x, y, random_x, random_y
    !     integer  :: seed_array(8)
    
    !     ! Initialize constants
    !     max = 19412.8
    !     a = 1e-6
    !     b = 1e-3
    !     c = (b - a) * max
    
    !     ! Initialize random number generator
    !     if (present(seed)) then
    !         seed_array(:) = seed
    !         call random_seed(PUT=seed_array) ! Optional: Initialize with a seed if needed
    !     else
    !         call random_seed()
    !     end if

    !     do
    !         ! Generate random numbers x in [0, 1) and y in [a, b)
    !         call random_number(random_x)
    !         x = random_x  ! Scale is not needed as random_number generates in [0, 1)
    
    !         call random_number(random_y)
    !         y = a + random_y * (b - a)  ! Scale to range [a, b)
    
    !         ! Accept or reject
    !         if (x <= (b - a) * normalized_N_r(dCd_r(dFdr(y), Vdp(vp4(y)))) / c) then
    !             exit
    !         end if
    !     end do
    
    ! end function acceptation_rejet

end module sample_mod