set title dirname
set logscale y
set datafile separator ","
set term pdf enhanced size 4in, 3in
set output "gnuplot-pdf.pdf"
plot "QR_1.0" using 3:4 w l lw 3
replot "OQDS" using 3:4 w l lw 3
replot "One-Sided_Jacobi" using 3:4 w l lw 3

