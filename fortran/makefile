FORT=gfortran
FFLAGS=-Wall -Wextra -pedantic -g -rdynamic
#FORT=g95
#FORT=f77
#FFLAGS=-Wall -Wextra -pedantic -g

test: bigden.o biglag.o calfun.o main.o newuob1.o newuob2.o newuob3.o newuob4.o trsubp.o update.o newuoa.o setup.o aux.o modules.o subnewuoa.o subnewuob.o
	$(FORT) $(FFLAGS) $(LIBS) -o test bigden.o biglag.o calfun.o main.o newuob1.o newuob2.o newuob3.o newuob4.o trsubp.o update.o  newuoa.o setup.o aux.o modules.o subnewuoa.o subnewuob.o

modules.o: modules.f
	$(FORT) $(FFLAGS) -c modules.f
prob.mod: modules.f
	$(FORT) $(FFLAGS) -c modules.f
permute.mod: modules.f
	$(FORT) $(FFLAGS) -c modules.f
bigden.o: bigden.f
	$(FORT) $(FFLAGS) -c bigden.f
biglag.o: biglag.f
	$(FORT) $(FFLAGS) -c biglag.f
calfun.o: calfun.f modules.o 
	$(FORT) $(FFLAGS) -c calfun.f 
main.o: main.f modules.o
	$(FORT) $(FFLAGS) -c main.f 
trsubp.o: trsubp.f
	$(FORT) $(FFLAGS) -c trsubp.f
newuob1.o: newuob1.f
	$(FORT) $(FFLAGS) -c newuob1.f
newuob2.o: newuob2.f
	$(FORT) $(FFLAGS) -c newuob2.f
newuob3.o: newuob3.f
	$(FORT) $(FFLAGS) -c newuob3.f
newuob4.o: newuob4.f
	$(FORT) $(FFLAGS) -c newuob4.f
update.o: update.f
	$(FORT) $(FFLAGS) -c update.f
newuoa.o: newuoa.f
	$(FORT) $(FFLAGS) -c newuoa.f
setup.o: setup.f
	$(FORT) $(FFLAGS) -c setup.f
aux.o: aux.f
	$(FORT) $(FFLAGS) -c aux.f
subnewuoa.o: subnewuoa.f
	$(FORT) $(FFLAGS) -c subnewuoa.f
subnewuob.o: subnewuob.f
	$(FORT) $(FFLAGS) -c subnewuob.f
clean:
	rm -f *.o; rm -f test; rm -f *.mod
