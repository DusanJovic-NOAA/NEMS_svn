
include       $(ESMFMKFILE)
ESMF_INC    = $(ESMF_F90COMPILEPATHS)
ESMF_LIB    = $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)

gfs=gsm
gfsdir=atmos/$(gfs)

NETCDF_LIB  =   $(NWPROD_LIB)/netcdf/lib/libnetcdf.a

NEMSIO_INC  = -I$(NWPROD_LIB)/incmod/nemsio
NEMSIO_LIB  = -L$(NWPROD_LIB) -lnemsio

W3_LIB      = -L$(NWPROD_LIB) -lw3nco_4 -lw3emc_4
BACIO_LIB   = -L$(NWPROD_LIB) -lbacio_4
SP_LIB      = -L$(NWPROD_LIB) -lsp_4
SYS_LIB     =

EXTLIBS     = $(NEMSIO_LIB) \
              $(W3_LIB) \
              $(BACIO_LIB) \
              $(SP_LIB) \
              $(ESMF_LIB) \
              $(NETCDF_LIB) \
              $(SYS_LIB)

FC          = mpif90 -fc=gfortran
FREE        = -ffree-form
FIXED       = -ffixed-form
R8          = -fdefault-real-8

FINCS       = $(ESMF_INC) $(NEMSIO_INC)
TRAPS       =

FFLAGS      = $(TRAPS) $(FINCS) -fconvert=big-endian -fno-range-check -ffree-line-length-256

OPTS_NMM    = -O2
OPTS_GFS    = -O2
OPTS_GEN    = -O2
OPTS_FIM    = -O2

FFLAGS_NMM  = $(OPTS_NMM) $(FFLAGS)
FFLAGS_GFS  = $(OPTS_GFS) $(FFLAGS) $(FREE)
FFLAGS_GFSF = $(OPTS_GFS) $(FFLAGS) $(FIXED)
FFLAGS_GEN  = $(OPTS_GEN) $(FFLAGS)
FFLAGS_FIM  = $(OPTS_FIM) $(FFLAGS)

CPP         = /lib/cpp -P -traditional
CPPFLAGS    =

AR          = ar
ARFLAGS     = -r

RM          = rm
