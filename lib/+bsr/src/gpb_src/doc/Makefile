# constants
empty := 
space := $(empty) $(empty)
comma := ,

# utility functions - lists
wordlist_tail = $(wordlist 2, $(words $(1)), $(1))
wordlist_last = $(word $(words $(1)), $(1))

# utility functions - directories
dirlist_tail = $(subst $(space),/, $(call wordlist_tail, $(subst /, , $(1))))
dirlist_last = $(call wordlist_last, $(subst /, , $(1)))

# compute paths and include global options
path_to_root_dir := 
path_from_root_dir := 
path_to_working_dir := $(CURDIR)/
include Rules.make
path_in_tree := $(subst //,/,$(call dirlist_tail, $(path_from_root_dir))/)

# default target
.PHONY: default
default: all

# rule for building all documentation
.PHONY: all
all: | clean
	$(DOXYGEN) Doxyfile

# rule for cleaning all documentation
.PHONY: clean
clean:
	rm -rf html/

.PHONY: clean_all
clean_all: | clean
