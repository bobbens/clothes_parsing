NUM_THREADS ?= 
ifneq ($(strip $(NUM_THREADS)),)
	NUM_THREADS_OPTION := -DNUM_THREADS=$(NUM_THREADS)
endif

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

# utility functions - find subdirectories of the specified directory
subdir_ignore_file := $(wildcard .ignore_dirs)
subdirs_ignore := $(foreach ignore_file, $(subdir_ignore_file), $(shell cat $(ignore_file)))
find_subdirs_containing_file = $(filter-out $(subdirs_ignore), $(dir $(wildcard $(addprefix $(addprefix $(1), */), $(2)))))
find_subdirs = $(call find_subdirs_containing_file, $(1), .)

# compute paths and include global options
path_to_root_dir := 
path_from_root_dir := 
path_to_working_dir := $(CURDIR)/
include Rules.make
path_in_matlab := $(subst //,/,$(call dirlist_tail, $(path_from_root_dir))/)

# find a list of object files to link by examining a depend object (.do) file
o_link_list_from_depend_o = $(patsubst \,,$(patsubst %:,,$(1)))

# compute subdirectories
subdirs := $(call find_subdirs, )

# compute subdirectories containing makerules/makefiles
subdirs_with_makerules := $(call find_subdirs_containing_file, , Rules.make)
subdirs_with_makefiles := $(call find_subdirs_containing_file, , Makefile)

# compute makefiles that should be found in each subdirectory
makerules := $(addsuffix Rules.make, $(subdirs))
makefiles := $(addsuffix Makefile, $(subdirs))

# compute compilation targets
target_names   := $(patsubst %.cc,%, $(wildcard *.cc))
depend_files   := $(addsuffix .d, $(target_names))
depend_o_files := $(addsuffix .do, $(target_names))
object_files   := $(addsuffix .o, $(target_names))
mex_files      := $(addsuffix .$(MEXEXT), $(target_names))

# default target - just rebuild current branch of tree
.PHONY: default
default: local

# rule for building everything in current branch of tree and lower
.PHONY: all
all: subdir_makerules subdir_makefiles local subcomponents

# debug this branch of the tree
.PHONY: debug
debug:
	@echo "current dir     =" $(CURDIR)
	@echo "path from root  =" $(path_from_root_dir)
	@echo "path to root    =" $(path_to_root_dir)
	@echo "path in matlab  =" $(path_in_matlab)
	@echo "ignored subdirs =" $(subdirs_ignore)
	@echo "subdirs         =" $(subdirs)
	@echo "sub makerules   =" $(makerules)
	@echo "sub makefiles   =" $(makefiles)
	@echo "target names    =" $(target_names)
	@echo "depend files    =" $(depend_files)
	@echo "depend o files  =" $(depend_o_files)
	@echo "object files    =" $(object_files)
	@echo "mex files       =" $(mex_files)

# rules for creating makefiles for subcomponents
.PHONY: subdir_makerules
.PHONY: subdir_makefiles
subdir_makerules: $(makerules)
subdir_makefiles: $(makefiles)
$(makerules):
	ln -s ../Rules.make $@
$(makefiles):
	ln -s ../Makefile $@

# rules for building local components
.PHONY: depend
.PHONY: local
.SECONDARY: $(depend_files) $(depend_o_files)
depend: $(depend_files) $(depend_o_files)
local: $(depend_files) $(depend_o_files) $(object_files) $(mex_files)
%.d: %.cc
	$(CXX) $(MEX_FLAGS) -iquote "$(path_to_root_dir)$(include_dir)" -I"$(MEX_INCLUDE_PATH)" $(NUM_THREADS_OPTION) -MM -MT $@ -o $@ $<
%.do: %.d
	python2 $(path_to_root_dir)$(build_scripts_dir)/depend_o.py $*.$(MEXEXT) $*.d $(path_to_root_dir) $(bin_dir) $(depend_dir) $(include_dir) $(src_dir)
%.o: %.d
	$(CXX) $(MEX_FLAGS) -iquote "$(path_to_root_dir)$(include_dir)" -I"$(MEX_INCLUDE_PATH)" $(NUM_THREADS_OPTION) -c -o $@ $*.cc
%.$(MEXEXT): %.do
	$(CXX) -o $@ $(MEX_LDFLAGS) $(call o_link_list_from_depend_o, $(shell cat $*.do))
ifneq ($(depend_files),)
-include $(depend_files)
endif
ifneq ($(depend_o_files),)
-include $(depend_o_files)
endif
        
# rules for building subcomponents
.PHONY: subcomponents $(subdirs)
subcomponents: $(subdirs)
$(subdirs): | subdir_makerules subdir_makefiles
	@$(MAKE) -r -C $@ all

# clean local branch of the tree
.PHONY: clean
clean:
	rm -f $(depend_files)
	rm -f $(depend_o_files)
	rm -f $(object_files)
	rm -f $(mex_files)
        
# clean lower branches of the tree
.PHONY: clean_subdirs
clean_subdirs: 
	@for dir in $(subdirs_with_makefiles); do \
		$(MAKE) -r -C $$dir clean_all; \
	done

# clean local branch of the tree and lower
.PHONY: clean_all
clean_all: | clean clean_subdirs
	rm -f $(makerules)
	rm -f $(makefiles)
