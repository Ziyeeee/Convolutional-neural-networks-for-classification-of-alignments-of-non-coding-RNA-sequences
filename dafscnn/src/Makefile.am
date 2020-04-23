SUBDIRS = contrafold probconsRNA contralign
bin_PROGRAMS = dafs

dafs_SOURCES = dafs.cpp typedefs.h fa.cpp fa.h \
	align.h align.cpp needleman_wunsch.h needleman_wunsch.cpp \
	fold.h fold.cpp nussinov.h nussinov.cpp ipknot.h ipknot.cpp \
	alifold.cpp alifold.h boltzmann_param.c boltzmann_param.h \
	ip.cpp ip.h cmdline.c cmdline.h
dafs_LDADD = contrafold/libcontrafold.a contralign/libcontralign.a probconsRNA/libprobcons.a
EXTRA_DIST = gnugetopt.h getopt.c getopt1.c

#dafs.cpp: cmdline.h
#cmdline.c cmdline.h: dafs.ggo
#	gengetopt --unamed-opts < dafs.ggo
