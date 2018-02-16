# shared compiler settings and rules
CXX ?= g++
CXXWARNINGS ?= -Wall -Wextra -Wundef -pedantic
CXXFLAGS ?= ${CXXWARNINGS} -g -O3

COMPILE = ${CXX} ${CXXFLAGS} ${CPPFLAGS} -c
LINK = ${CXX} ${LDFLAGS}
.DEFAULT_GOAL := ${TARGET}

# link executables from object files
%: %.o
	${LINK} -o $@ $^ ${LDLIBS}

# object files from cpp files
%.o: %.cpp
	${COMPILE} -o $@ $<

# Dependency determination.
Makefile.dep: [^_]*.cpp
	${COMPILE} -MM $^ >$@

-include Makefile.dep