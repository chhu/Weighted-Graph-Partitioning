CPPFLAGS=-O2 -march=native 
LDFLAGS=
LDLIBS=
CXX=g++

all: gpart mapinfo 

gpart: gpart.o
	$(CXX) $(LDFLAGS) -o gpart gpart.o $(LDLIBS)
 
mapinfo: mapinfo.o
	$(CXX) $(LDFLAGS) -o mapinfo mapinfo.o $(LDLIBS)
 
clean:
	$(RM) gpart gpart.o mapinfo mapinfo.o
