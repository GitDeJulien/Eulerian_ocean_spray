program EulerianOceanSpray

    use init_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df
    type(MeshType)     :: me

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    !Initialize the mesh
    call init_mesh(df, me)
    
    !Initialize de solution
    call initialize_sol(df, me)

    !!TODO > Boucle en temps + calcul moyenne

    call free_mesh(me)
    
end program EulerianOceanSpray