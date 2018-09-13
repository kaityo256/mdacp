#------------------------------------------------------
# Makefile for MDACP
#------------------------------------------------------

VERSION=2.20
TARGET=mdacp
RELEASEDIR=mdacp-$(VERSION)

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=mpic++
CPPFLAGS=-O3 -fopenmp -std=c++11
LDFLAGS=

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Source Files
#------------------------------------------------------

.SUFFIXES: .c .cc .h. .o
.PHONY: clean dep

SRC=$(shell ls *.cc)
HED=$(shell ls *.h)
OBJ=$(SRC:.cc=.o)

#------------------------------------------------------
# Rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ)

.cc.o:
	$(CC) $(CPPFLAGS) -c $< 

dep:
	g++ -MM -MG $(SRC) >makefile.depend

makefile.depend: 
	g++ -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(OBJ) gmon.*.out gmon.out

tar:
	tar cvzf $(TARGET).tar.gz *.cfg $(SRC) $(HED) makefile

release:
	mkdir $(RELEASEDIR)
	cp $(SRC) $(HED) LICENSE makefile *.cfg README makefile.opt.* $(RELEASEDIR)
	tar cvzf $(RELEASEDIR).tar.gz $(RELEASEDIR)

#--------------------------------------------------
-include makefile.depend
