# build gentables app
ROOTDIR := ../..
TARGET := $(ROOTDIR)/gentables
include ${ROOTDIR}/Makefile.inc

CXXFLAGS += -I${ROOTDIR}/include

SOURCES := $(wildcard [^_]*.cpp)
${TARGET}: ${SOURCES:.cpp=.o}
	$(LINK) $^ $(LDLIBS) -o $@

clean:
	${RM} ${TARGET} *.o Makefile.dep

.PHONY: clean