LIBDIR=../lib
CXXFLAGS= -O -I. -std=c++11
CXX= g++
target=$(LIBDIR)/libmodel_MRSSMQCD_UFO_from_philip.a

all: $(target)

objects=HelAmps_MRSSMQCD_UFO_from_philip.o Parameters_MRSSMQCD_UFO_from_philip.o rambo.o read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

