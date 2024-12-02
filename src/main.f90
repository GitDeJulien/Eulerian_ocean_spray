program EulerianOceanSpray

    use init_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df
    type(StructCelleType), dimension(:), allocatable :: celles

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    ! Initialize the mesh
    call init_mesh(df, celles)
    
    ! Initialize de solution
    call initial_step(df, celles)

    print*, celles(1)%SOL(1,1,1,:)

    deallocate(celles)
    
end program EulerianOceanSpray