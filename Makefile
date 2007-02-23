CC      = gcc
CCFLAGS = -O3 -Wall -DINTELSWP

INCPATH = src/include
OBJDIR  = obj
BINDIR  = bin
TARGET  = mp5

#          $(OBJDIR)/mmp1.o        \
#          $(OBJDIR)/mmp2.o        \
#          $(OBJDIR)/mmp3.o        \

#          $(OBJDIR)/mmp3.o        \

OBJECTS = $(OBJDIR)/cmd.o         \
          $(OBJDIR)/dic.o         \
          $(OBJDIR)/gabor.o       \
          $(OBJDIR)/io.o          \
          $(OBJDIR)/main.o        \
          $(OBJDIR)/mmp3.o        \
          $(OBJDIR)/mp5.o         \
          $(OBJDIR)/matrix.o      \
          $(OBJDIR)/new_io.o      \
          $(OBJDIR)/queue.o       \
          $(OBJDIR)/r250.o        \
          $(OBJDIR)/smp.o         \
          $(OBJDIR)/stringTools.o \
          $(OBJDIR)/tools.o       \
          $(OBJDIR)/vector.o

.PHONY: all

all: mp5

$(TARGET): $(OBJECTS)
	$(CC) $^ $(LIBS) -lm -o $(BINDIR)/$@

.PHONY: clean

clean:
	rm -f obj/*.o

# DEPENDIECES AND COMPILATION #

$(OBJDIR)/cmd.o: src/cmd.c                 \
                 src/include/cmd.h         \
                 src/include/def.h         \
                 src/include/queue.h       \
                 src/include/stringTools.h \
                 src/include/types.h       \
                 src/include/vector.h
	$(CC) -c -I $(INCLPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/dic.o: src/dic.c            \
                 src/include/def.h    \
                 src/include/dic.h \
                 src/include/gabor.h  \
                 src/include/r250.h   \
                 src/include/types.h  \
                 src/include/vector.h  
	
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/gabor.o: src/gabor.c         \
                   src/include/gabor.h \
                   src/include/types.h \
                   src/include/vector.h 
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/io.o: src/io.c                      \
                src/include/def.h             \
                src/include/gabor.h           \
                src/include/io.h              \
                src/include/matrix.h          \
                src/include/new_io.h          \
                src/include/queue.h           \
                src/include/types.h           \
                src/include/vector.h
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<
      
$(OBJDIR)/main.o: src/main.c                    \
                  src/include/cmd.h             \
                  src/include/def.h             \
                  src/include/dic.h             \
                  src/include/io.h              \
                  src/include/mmp3.h            \
                  src/include/mp5.h             \
                  src/include/smp.h             \
                  src/include/stringTools.h     \
                  src/include/types.h           
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/matrix.o: src/matrix.c                \
                    src/include/matrix.h
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

#$(OBJDIR)/mmp1.o: src/mmp1.c                    \
#                  src/include/def.h             \
#                  src/include/gabor.h           \
#                  src/include/mp5.h             \
#                  src/include/queue.h           \
#                  src/include/tools.h           \
#		  src/include/types.h           \
#		  src/include/vector.h           
#	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

#$(OBJDIR)/mmp2.o: src/mmp2.c          \
#                  src/include/mp5.h   \
#                  src/include/types.h \
#                  src/include/def.h
#	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/mmp3.o: src/mmp3.c                    \
                  src/include/def.h             \
                  src/include/gabor.h           \
                  src/include/mmp3.h            \
                  src/include/mp5.h             \
                  src/include/queue.h           \
                  src/include/tools.h           \
                  src/include/types.h           \
                  src/include/vector.h           
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/mp5.o: src/mp5.c                     \
                 src/include/def.h             \
                 src/include/dic.h             \
                 src/include/gabor.h           \
                 src/include/matrix.h          \
                 src/include/mp5.h             \
                 src/include/queue.h           \
                 src/include/tools.h           \
                 src/include/types.h           \
                 src/include/vector.h
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/new_io.o: src/new_io.c         \
                    src/include/new_io.h
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<


$(OBJDIR)/queue.o: src/queue.c         \
                   src/include/def.h   \
		   src/include/types.h
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/r250.o: src/r250.c           \
                  src/include/r250.h   
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/smp.o: src/smp.c                     \
                 src/include/def.h             \
                 src/include/gabor.h           \
                 src/include/mp5.h             \
                 src/include/queue.h           \
                 src/include/smp.h             \
                 src/include/tools.h           \
                 src/include/types.h        
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/stringTools.o: src/stringTools.c         \
                         src/include/def.h         \
                         src/include/stringTools.h 
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/tools.o: src/tools.c         \
                   src/include/def.h   \
                   src/include/tools.h 
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<

$(OBJDIR)/vector.o: src/vector.c          \
                    src/include/vector.h 
	$(CC) -c -I $(INCPATH) $(CCFLAGS) -o $@ $<
