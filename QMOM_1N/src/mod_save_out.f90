module save_out_mod

    use structured_mesh_mod
    implicit none
    
contains

    subroutine save_approx_sol(df, r_tab, T_tab, n_bar, u_bar, iter)

        !In
        type(DataType), intent(in)           :: df
        real(pr), dimension(:), intent(in)   :: r_tab, T_tab
        real(pr), dimension(:,:), intent(in) :: n_bar, u_bar
        integer, intent(in)                  :: iter

        !Local
        integer :: i, j, io, ios
        character(len=125) :: ch_iter

        write(ch_iter, '(I5)') iter

        open(newunit=io, file="output/sol/sol."//trim(adjustl(ch_iter))//".dat",&
                status="replace", action="write", iostat=ios)

            if (ios /= 0) then
                print *, 'Error opening file: ', " output/sol/sol.*"
                stop
            end if

            write(io, *) "## r(i)  ", "  T(j)  ", "  n_bar(i,j)  ", "  u_bar(i,j)"
            do i=1,df%N_r
                do j=1,df%N_T
                    write(io,*) r_tab(i), T_tab(j)-273.15d0, n_bar(i,j), u_bar(i,j)
                enddo
            enddo

        close(io)

    end subroutine save_approx_sol

    subroutine save_mesh(mesh_tab, name)

        !In
        real(pr), dimension(:), intent(in) :: mesh_tab
        character(len=*), intent(in)       :: name

        !Local
        integer :: i, io, ios

        open(newunit=io, file="output/other/"//trim(adjustl(name))//".mesh.dat",&
                status="replace", action="write", iostat=ios)

            if (ios /= 0) then
                print *, 'Error opening file: ', " output/sol/sol.*"
                stop
            end if

            write(io, *) "## "//trim(adjustl(name))
            do i=1,size(mesh_tab,dim=1)
                write(io,*) mesh_tab(i)
            enddo

        close(io)

    end subroutine save_mesh
    
end module save_out_mod