program EulerianOceanSpray

    use data_mod
    use structured_mesh_mod
    use init_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df
    type(MeshType)     :: me
    ! integer :: i,j,l
    ! real(pr)           :: tn, dt
    ! integer            :: nt

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    !Initialize the mesh
    call init_mesh(df, me)

    !Initialize coefficient
    call initialize_coeff(df, me)

    print*, me%R_coeff
    
    !Initialize de solution
    call initialize_sol(df, me)

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

    ! !TODO > Boucle en temps + calcul moyenne
    ! tn = df%t0
    ! do nt=1,df%ntime
    !     call advance(df, me, dt)
    !     print*, dt
    !     tn = tn + dt
    ! enddo

    call free_mesh(me)
    
end program EulerianOceanSpray