noinst_LIBRARIES = libpartalign.a
bin_PROGRAMS=pa_train

libpartalign_a_SOURCES = partalign.cpp partalign.h log_value.h

pa_train_SOURCES = pa_train.cpp optimizer.cpp optimizer.h lbfgsb.c lbfgsb.h ../fa.cpp ../fa.h ../fold.cpp ../fold.h ../boltzmann_param.c ../boltzmann_param.h
pa_train_LDADD = libpartalign.a ../contrafold/libcontrafold.a
