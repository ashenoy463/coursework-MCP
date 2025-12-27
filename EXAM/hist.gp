# https://psy.swansea.ac.uk/staff/carter/gnuplot/gnuplot_frequency.htm

clear
reset
set key off
#set border 3

# Each bar is half the (visual) width of its x-range.

set boxwidth 0.05 absolute
set style fill solid 1.0 noborder

bin_width = 0.1;

#set xrange [-8:8]

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * ( bin_number(x) + 0.5 )

plot 'DHIST.dat' using (rounded($1)):(1) smooth frequency with boxes
