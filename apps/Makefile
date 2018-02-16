# recursively builds apps contained in this folder
all:
	$(MAKE) -C example
	$(MAKE) -C gentables

clean:
	$(MAKE) -C example $@
	$(MAKE) -C gentables $@

.PHONY: all clean