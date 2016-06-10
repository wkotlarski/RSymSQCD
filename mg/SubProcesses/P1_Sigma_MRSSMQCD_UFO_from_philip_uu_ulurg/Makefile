LIBDIR=../../lib
INCDIR=../../src
MODELLIB=model_MRSSMQCD_UFO_from_philip
CXXFLAGS= -O -I$(INCDIR) -I. -std=c++11
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)
CXX= g++#-mp-4.9
main=check

all: $(main)

objects=Process_uu_ulurg.o check_sa.o

$(main): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -o $@ $(objects) $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	cd ../../src
	make

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

