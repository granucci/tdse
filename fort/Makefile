# Makefile created by mkmf.perl $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include ../make.defs

.DEFAULT:
	-touch $@
all: quant_ev rd_resu
diag.o: ./diag.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./diag.f
lapack_zheev.o: ./lapack_zheev.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./lapack_zheev.f
quant_efield.o: ./quant_efield.f90 quant_util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quant_efield.f90
quant_ev.o: ./quant_ev.f90 quant_util.o quant_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quant_ev.f90
quant_mod.o: ./quant_mod.f90 quant_util.o quant_sys.o quant_efield.o singleton.o fft_drv.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quant_mod.f90
quant_sys.o: ./quant_sys.f90 quant_util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quant_sys.f90
quant_util.o: ./quant_util.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quant_util.f90
rd_resu.o: ./rd_resu.f90 quant_util.o singleton.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./rd_resu.f90
singleton.o: ./singleton.f90 quant_util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./singleton.f90
fft_drv.o: ./fft_drv.F90 quant_util.o singleton.o
	$(FC) $(FFLAGS) $(CPPDEFS) $(OTHERFLAGS) -c	./fft_drv.F90
SRC = ./quant_efield.f90 ./lapack_zheev.f ./quant_sys.f90 ./diag.f ./quant_util.f90 ./quant_mod.f90 ./singleton.f90 ./fft_drv.F90 ./quant_ev.f90
OBJ = quant_efield.o lapack_zheev.o quant_sys.o diag.o quant_util.o singleton.o fft_drv.o quant_mod.o quant_ev.o
clean: neat
	-rm -f .cppdefs $(OBJ) quant_ev rd_resu
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
quant_ev: $(OBJ)
	$(LD) $(OBJ) -o quant_ev $(LDFLAGS)
rd_resu: rd_resu.o quant_util.o singleton.o
	$(LD) rd_resu.o quant_util.o singleton.o -o rd_resu $(LDFLAGS)
