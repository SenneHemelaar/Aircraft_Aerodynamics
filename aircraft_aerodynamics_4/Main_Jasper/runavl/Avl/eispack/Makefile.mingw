# Makefile for MinGW (with gfortran)

FC = gfortran
FLG = -O
AR = ar r
RAN = ranlib

##---------------------------------------
##  Uncomment for MingGW with gfortran, links to Win32
### GNU gfortran
FC = gfortran
# add -m32 or -m64 flag to generate 32 bit or 64 bit for your environment  
#FC = gfortran -m64  

FLG = -O2 -fomit-frame-pointer
# Debug flags (symbols, array bounds)
#FLG = -g -O0 -C

### single precision
SRC = eispack.f
BIN = eispack.o
DP =
EISLIB = eispack_gSP.a

### double precision
#SRC = eispackDP.f
#BIN = eispackDP.o
#DP = -r8
#EISLIB = eispack_gDP.a


$(EISLIB): $(SRC)
	$(FC) -c $(FLG) $(DP) $(SRC)
	$(AR) $(EISLIB) $(BIN)
	$(RAN) $(EISLIB)

clean:
	/bin/rm *.o *.a


