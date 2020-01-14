C resources
http://anythingbutrbitrary.blogspot.com/2014/02/?m=1
http://www.biostat.jhsph.edu/~rpeng/docs/interface.pdf
https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Handling-R-objects-in-C
http://adv-r.had.co.nz/C-interface.html
http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_11.html

Compiling:
R CMD SHLIB <jsij>.c  -L/usr/local/lib/ -lgsl -lgslcblas\n

