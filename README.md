GP here
http://krasserm.github.io/2018/03/19/gaussian-processes/
Better here
https://katbailey.github.io/post/gaussian-processes-for-dummies/

Stuff to read about C
http://anythingbutrbitrary.blogspot.com/2014/02/?m=1
http://www.biostat.jhsph.edu/~rpeng/docs/interface.pdf
https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Handling-R-objects-in-C
http://adv-r.had.co.nz/C-interface.html
http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_11.html

Compiling
R CMD SHLIB <jsij>.c  -I/usr/local/lib/ -lgsl -lgslcblas
