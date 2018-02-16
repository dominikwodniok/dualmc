# build apps in apps subdir
all: $(BINDIR)
	$(MAKE) -C apps

clean:
	$(MAKE) -C apps $@

.PHONY: all clean