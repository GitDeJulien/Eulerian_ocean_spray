module save_out_mod

    use structured_mesh_mod
    implicit none
    
contains

    subroutine save_approx_sol(df, r_tab, n_bar, u_bar, iter)

        !In
        type(DataType), intent(in)           :: df
        real(pr), dimension(:), intent(in)   :: r_tab, n_bar, u_bar
        integer, intent(in)                  :: iter

        !Local
        integer :: i, io, ios
        character(len=125) :: ch_iter

        write(ch_iter, '(I5)') iter

        open(newunit=io, file="output/sol/sol."//trim(adjustl(ch_iter))//".dat",&
                status="replace", action="write", iostat=ios)

            if (ios /= 0) then
                print *, 'Error opening file: ', " output/sol/sol.*"
                stop
            end if

            write(io, *) "## r(i)  ", "  n_bar(i)  ", "  u_bar(i)"
            do i=1,df%N_r
                write(io,*) r_tab(i), n_bar(i), u_bar(i)
            enddo

        close(io)

    end subroutine save_approx_sol
    
end module save_out_mod