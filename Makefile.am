# Keep everything that is generated.
.SECONDARY:
# Remove works-in-progress.
.DELETE_ON_ERROR:

mp5_CPPFLAGS = -I$(top_srcdir)/src/include/ -DINTELSWP -DINLINE

mp5_CFLAGS = $(FFTW3_CFLAGS)

mp5_LDADD = $(FFTW3_LIBS)

bin_PROGRAMS=mp5
mp5_SOURCES = src/atom.c \
	src/cmd.c \
	src/dic.c \
	src/io_mp5.c \
	src/main.c \
	src/matrix.c \
	src/mmp.c \
	src/mp5.c \
	src/queue.c \
	src/r250.c \
	src/smp.c \
	src/stringTools.c \
	src/tools.c \
	src/vector.c \
	src/random.c \
	src/mtwist.c \
	src/include/atom.h \
	src/include/cmd.h \
	src/include/def.h \
	src/include/dic.h \
	src/include/io_mp5.h \
	src/include/matrix.h \
	src/include/mmp.h \
	src/include/mp5.h \
	src/include/queue.h \
	src/include/r250.h \
	src/include/smp.h \
	src/include/stringTools.h \
	src/include/tools.h \
	src/include/types.h \
	src/include/vector.h \
	src/include/random.h \
	src/include/mtwist.h

install-exec-local:
	test -z "$(FFTW3DLL)" || $(INSTALL) "$(FFTW3DLL)" "$(DESTDIR)$(bindir)"
