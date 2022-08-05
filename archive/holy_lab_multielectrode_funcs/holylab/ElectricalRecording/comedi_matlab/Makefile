include version
LIBVERSION_MAJOR=0
THISDIR=`pwd`
MEXEXT=mexglx

include Config

CFLAGS += -I . -I $(MEXHDIR)
LDFLAGS += -lcomedi -lm 
LDFLAGS += -lcomedipersist -L$(INSTALLDIR_LIB) -I$(INSTALLDIR_LIBH)

lib:
	g++ -fPIC -c comedipersist.cc
	g++ -shared -o libcomedipersist.so.$(libversion) comedipersist.o

lib_install:
	@if [ ! -d $(INSTALLDIR_LIB) ]; then echo "Library install directory does not exist"; exit 1; fi
	if [ ! -d $(INSTALLDIR_LIBH) ]; then install -d $(INSTALLDIR_LIBH); fi
	install -m 644 comedipersist.h $(INSTALLDIR_LIBH)
	install libcomedipersist.so.$(libversion) $(INSTALLDIR_LIB)
	(cd $(INSTALLDIR_LIB); ln -sf libcomedipersist.so.$(libversion) libcomedipersist.so)

commono = commat_util.o
commonh = commat.h

$(commono): $(commonh)
	@if [ ! -f $(MEXHDIR)/mex.h ]; then echo "Path to mex.h incorrect; edit Config file"; exit 1; fi
	$(CC) $(CFLAGS) -c $(@:.o=.cc)

matfuncs = comedi_open comedi_isopen comedi_close comedi_data_read comedibuf_read comedi_find_range comedi_get_range comedi_command_test comedi_command comedi_get_n_ranges comedi_get_max_buffer_size comedi_set_buffer_size comedi_cancel_and_flush comedi_get_subdevice_type comedi_find_subdevice_by_type comedi_dio_config comedi_dio_bitfield comedi_get_n_channels comedibuf_write comedi_data_write comedi_do_insnlist comedi_internal_trigger 

$(matfuncs): % : $(commonh) $(commono:.o=.cc)
	mex $(LDFLAGS) $@.cc $(commono)

allmex: commat_util.o $(matfuncs)

test: % : $(commonh) $(commono:.o=.cc)
	mex $(LDFLAGS) $@.cc $(commono)


.PHONY: clean
clean:
	-rm *.o lib*.so* $(addsuffix .$(MEXEXT),$(matfuncs))
