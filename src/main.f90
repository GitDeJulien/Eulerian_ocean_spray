program EulerianOceanSpray

    use init_mod
    use time_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df
    type(MeshType)     :: me
    integer :: i,j,k,l
    real(pr)           :: tn, dt
    integer            :: nt

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    !Initialize the mesh
    call init_mesh(df, me)

    !Initialize coefficient
    call initialize_coeff(df, me)
    
    !Initialize de solution
    call initialize_sol(df, me)

    ! do i=1,df%N_T+1 
    !     print*, me%T_tab(i) 
    ! end do

    ! print*, R_function(df, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0)

    ! print*, df%Dv_star_function(df, 1.0, 1.0)

    ! do i=1,df%N_r !radius
    !     ! print*, (me%r_tab(i+1) - me%r_tab(i))/2.d0
    !     ! print*, me%m_sel(i)
    !     do j=1,df%N_vx !velocity
    !         ! print*, me%vx_tab(j)
    !         do k=1,df%N_m
    !             do l=1,df%N_T !Temperature
    !                 ! print*, me%T_tab(l)
    !                 ! if (me%SOL(i,j,i,l) > 1e-8) print*, me%SOL(i,j,i,l)
    !                 print*, me%V_coeff(i,j,k,l)
    !                 ! print*, me%r_tab(i) + (me%r_tab(i+1) - me%r_tab(i))/2.d0
    !             enddo
    !         enddo
    !     enddo
    ! enddo

    ! print*, me%T_coeff(df%N_r,df%N_vx,df%N_m,df%N_T)

    ! print*,"r_init"
    ! print*,me%r_tab
    ! print*,"vx_init"
    ! print*,me%vx_tab
    ! print*,"m_init"
    ! print*,me%m_tab
    ! print*,"T_init"
    ! print*,me%T_tab
    ! print*, "SOL!=0.0"

    ! do i=1,df%N_r !radius
    !     do j=1,df%N_vx !velocity
    !         do l=1,df%N_T !Temperature
            
    !             if (me%SOL(i,j,i,l) > 1e-8) print*, me%SOL(i,j,i,l)

    !         enddo
    !     enddo
    ! enddo

    !TODO > Boucle en temps + calcul moyenne
    tn = df%t0
    do nt=1,df%ntime
        call advance(df, me, dt)
        tn = tn + dt
        print*, "tn = ", tn
    enddo

    do i=1,df%N_r !radius
        do j=1,df%N_vx !velocity
            do k=1,df%N_m
                do l=1,df%N_T !Temperature
                    if (me%SOL(i,j,k,l) > 1e-8) print*, me%SOL(i,j,k,l)
                enddo
            enddo
        enddo
    enddo

    call free_mesh(me)
    
end program EulerianOceanSpray