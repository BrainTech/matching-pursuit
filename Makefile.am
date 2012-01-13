# Keep everything that is generated.
.SECONDARY:
# Remove works-in-progress.
.DELETE_ON_ERROR:

mp5_CPPFLAGS = -I$(top_srcdir)/src/include/ -DINTELSWP -DINLINE

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
	src/vector.c
