MEX= mex
# You may have to specify explicitly the path of the MATLAB mex program on
# your computer. Ask a linux user if you do not understand what this means.

MFLAG= -g -compatibleArrayDims
# The option -compatibleArrayDims is for the compatability under the
# latest releases of MATLAB (i.e., R2018a). If your MATLAB does not
# recognize this option, remove "-compatibleArrayDims".

FSRC=./Fsrc

all: evalfun setuptest

evalfun: evalfun.f evalfun_mexgateway.F
	$(MEX) -output evalfun $(MFLAG) evalfun.f evalfun_mexgateway.F

setuptest: setuptest.f setuptest_mexgateway.F
	$(MEX) -output setuptest $(MFLAG) setuptest.f setuptest_mexgateway.F

clean:
	rm -f *.o *.mod evalfun.mex* setuptest.mex* fort.6
