FSRC=./Fsrc

newuoa: $(FSRC)/*.f newuoa_mexgateway.F
	mex -output newuoa $(FSRC)/newuoa.f $(FSRC)/newuob.f $(FSRC)/bigden.f $(FSRC)/biglag.f $(FSRC)/trsapp.f $(FSRC)/update.f newuoa_mexgateway.F
	rm global_var.mod

clean: 
	rm -f *.o *.mod newuoa.mex*
