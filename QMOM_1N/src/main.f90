program EulerianOceanSpray

    use init_mod
    use time_mod
    use save_out_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df
    type(MeshType)     :: mesh
    ! integer :: i,j,l
    real(pr)           :: tn, dt
    integer            :: nt, io, ios

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    !Initialize the mesh
    call init_mesh(df, mesh)
    
    !Initialize de solution
    call initialize_sol(df, mesh)
    call save_approx_sol(df, mesh%r_tab, mesh%T_tab, mesh%n_bar, mesh%u_bar, 0)
    call save_mesh(mesh%r_tab, "radius")
    call save_mesh(mesh%T_tab, "Temperature")

    open(newunit=io, file="output/other/tn.dat", status="replace", action="write", iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', " output/sol/time.dat"
        stop
    end if
    write(io,*) 0.d0

    !!Boucle en temps
    tn = df%t0
    do nt=1,df%ntime
        call advance(df, mesh, dt)
        !print*, "diff =", abs(sum(mesh%u_bar(:)/df%N_r) - df%U_air)
        tn = tn + dt
        print*, "tn=", tn
        write(io,*) tn
        call save_approx_sol(df, mesh%r_tab, mesh%T_tab, mesh%n_bar, mesh%u_bar, nt)
    enddo

    print*, 'tfinal =', tn
    close(io)

    call free_mesh(mesh)
    
end program EulerianOceanSpray