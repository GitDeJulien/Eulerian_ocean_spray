program EulerianOceanSpray

    use data_mod
    implicit none

    character(len=256) :: filepath
    type(DataType)     :: df

    filepath = 'data/data.toml'

    ! Display TOML data file
    call display_toml_file(filepath)

    ! Read and stock all data from TOML data file
    call config_data(df, filepath)

    
end program EulerianOceanSpray