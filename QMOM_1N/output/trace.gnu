#set zrange [0.0:0.2]
#set xrange [-21.0:21.0]

set xlabel "radius"

#pour faire un gif
set term gif animate 
set output "sol_u.gif"
set logscale x

do for [i = 0:2000:10] {
    set title "iter = ".sprintf("%d", i)
    show title

    #-> 1D
    plot "sol/sol.".i.".dat" u 1:3 title "u(t,r)" w l lc 8 lw 2
    #plot "sol/sol.".i.".dat" u 1:2 title "n(t,r)" w l lc 8


}

unset term
unset view