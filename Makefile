CXX = g++
CXX_FLAGS = -Wall -O3
INPUTDIR = src
OBJDIR = bin

EXAMPLEBIN = dmc
GENTABLEBIN = gentable
EXAMPLEOBJ = $(OBJDIR)/main.o $(OBJDIR)/example.o
GENTABLEOBJ = $(OBJDIR)/gentable.o
HEADERS = $(INPUTDIR)/dualmc.h $(INPUTDIR)/dualmc.tpp $(INPUTDIR)/dualmc_table.tpp

all: $(OBJDIR) $(EXAMPLEOBJ) $(GENTABLEOBJ) $(HEADERS)
	$(CXX) $(CXX_FLAGS) -o $(EXAMPLEBIN) $(EXAMPLEOBJ)
	$(CXX) $(CXX_FLAGS) -o $(GENTABLEBIN) $(GENTABLEOBJ)

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: $(INPUTDIR)/%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

.PHONY: clean all
clean:
	rm -rf $(EXAMPLEBIN) $(GENTABLEBIN) $(EXAMPLEOBJ) $(GENTABLEOBJ)
