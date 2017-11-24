CXX = g++
CXX_FLAGS = -Wall -O3
INPUTDIR = src
OBJDIR = bin

BIN = dmc
OBJ = $(OBJDIR)/main.o $(OBJDIR)/example.o
HEADERS = $(INPUTDIR)/dualmc.h $(INPUTDIR)/dualmc.tpp $(INPUTDIR)/dualmc_table.tpp

all: $(OBJDIR) $(OBJ) $(HEADERS)
	$(CXX) $(CXX_FLAGS) -o $(BIN) $(OBJ)

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: $(INPUTDIR)/%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

.PHONY: clean all
clean:
	rm -rf $(BIN) $(OBJ)
