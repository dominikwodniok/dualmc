# build dual marching cubes example app
ROOTDIR := ../..
TARGET := $(ROOTDIR)/dmc
include ${ROOTDIR}/Makefile.inc

CXXFLAGS += -I${ROOTDIR}/include

SOURCES := $(wildcard [^_]*.cpp)
${TARGET}: ${SOURCES:.cpp=.o}
	$(LINK) $^ $(LDLIBS) -o $@

clean:
	${RM} ${TARGET} *.o Makefile.dep

.PHONY: clean