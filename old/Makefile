all: testfun setup 

testfun: testfun.f testfun_mexgateway.F
	mex -output testfun testfun.f testfun_mexgateway.F

setup: setup.f setup_mexgateway.F
	mex -output setup setup.f setup_mexgateway.F

clean: 
	rm -f *.o *.mod testfun.mex* setup.mex*
