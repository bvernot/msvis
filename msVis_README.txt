msVis is designed to take ms command line arguments exactly as ms would.  So for example, instead of:
ms 15 1000 -t 2.0 -I 3 10 4 1 5.0
You can run:
python2 msVis.py 15 1000 -t 2.0 -I 3 10 4 1 5.0

For options, run:
python2 msVis.py -h

msVis requires matplotlib to draw plots - the easiest way to install matplotlib is to install the basic version of the EPD Python distribution (click on the appropriate link for your operating system, then run the installation program):
https://www.enthought.com/repo/free/

I will probably add a flag later that allows you to run it without matplotlib, and it will just check your command line arguments.  It also currently prints out a lot of debug information, that will be hidden in future versions.

There are some commands that do not work yet: specifically, -f, -eM or -ema commands.  It does work with all of the examples used in the ms documentation.  Because I would like this to also serve as a way to check for proper ms input, there are some commands (-t, -s, -c, -r, etc) that don't affect the output, but are accepted, or even required.  I haven't yet implemented all of the command line checking, so if you notice something 'wrong' that doesn't raise a legible error, please let me know.  It definitely doesn't check for errors in the order of commands (e.g., -eg .1 1 10 -en .1 1 1), but it will in the future.  A number of non-ms commands are provided to tweak plotting options (use -h for details).

=============
Sample plots

Most of the sample plots are taken from the examples in msdoc.pdf.

The plot ea_aa.pdf is based on a demographic model of Africans and Europeans, and is generated with the command:

python msVis.py 36 1 -t 3.362600e+01 -r 14.6197076 50000 -I 2 18 18 0  -n 1 58.002735978 -n 2 70.041039672 -eg 0 1 482.46 -eg 0 2 570.18 -em 0 1 2 0.7310 -em 0 2 1 0.7310 -eg 0.006997264 1 0 -eg 0.006997264 2 89.7668 -en 0.006997264 1 1.98002736 -en 0.031463748 2 0.141176471 -en 0.03146375 2 0.254582763 -em 0.03146375 1 2 4.386 -em 0.03146375 2 1 4.386 -em 0.069767442 1 2 0 -em 0.069767442 2 1 0 -ej 1.751026e-01 2 1 -en 1.751026e-01 1 1.98002736  -en 0.20246238 1 1 -o ea_aa.pdf -pops Afr Eur -years -N0 7310