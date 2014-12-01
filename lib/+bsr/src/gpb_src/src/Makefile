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
path_in_tree := $(subst //,/,$(call dirlist_tail, $(path_from_root_dir))/)

# compute subdirectories for storage of dependency and binary output
subdirs        := $(call find_subdirs, )
depend_subdirs := $(addprefix $(path_to_root_dir)$(depend_dir)$(path_in_tree), $(subdirs))
bin_subdirs    := $(addprefix $(path_to_root_dir)$(bin_dir)$(path_in_tree), $(subdirs))

# compute subdirectories containing makerules/makefiles
subdirs_with_makerules := $(call find_subdirs_containing_file, , Rules.make)
subdirs_with_makefiles := $(call find_subdirs_containing_file, , Makefile)

# compute makefiles that should be found in each subdirectory
makerules := $(addsuffix Rules.make, $(subdirs))
makefiles := $(addsuffix Makefile, $(subdirs))

# compute compilation targets
target_names := $(patsubst %.cc,%, $(wildcard *.cc))
depend_files := $(addprefix $(path_to_root_dir)$(depend_dir)$(path_in_tree), $(addsuffix .d, $(target_names)))
object_files := $(addprefix $(path_to_root_dir)$(bin_dir)$(path_in_tree), $(addsuffix .o, $(target_names)))

# default target - just rebuild current branch of tree
.PHONY: default
default: local

# rule for building everything in current branch of tree and lower
.PHONY: all
all: depend_directories bin_directories subdir_makerules subdir_makefiles local subcomponents

# debug this branch of the tree
.PHONY: debug
debug:
	@echo "current dir     =" $(CURDIR)
	@echo "path from root  =" $(path_from_root_dir)
	@echo "path to root    =" $(path_to_root_dir)
	@echo "path in tree    =" $(path_in_tree)
	@echo "ignored subdirs =" $(subdirs_ignore)
	@echo "subdirs         =" $(subdirs)
	@echo "depend subdirs  =" $(depend_subdirs)
	@echo "bin subdirs     =" $(bin_subdirs)
	@echo "sub makerules   =" $(makerules)
	@echo "sub makefiles   =" $(makefiles)
	@echo "target names    =" $(target_names)
	@echo "depend files    =" $(depend_files)
	@echo "object files    =" $(object_files)
        
# rules for creating depend and bin directories
.PHONY: depend_directories
.PHONY: bin_directories
depend_directories: $(depend_subdirs)
bin_directories: $(bin_subdirs)
$(depend_subdirs) $(bin_subdirs):
	mkdir $@

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
depend: $(depend_files)
local: $(depend_files) $(object_files)
$(path_to_root_dir)$(depend_dir)$(path_in_tree)%.d: %.cc
	$(CXX) $(CXX_FLAGS) -iquote "$(path_to_root_dir)$(include_dir)" -MM -MT $@ -o $@ $<
$(path_to_root_dir)$(bin_dir)$(path_in_tree)%.o: $(path_to_root_dir)$(depend_dir)$(path_in_tree)%.d
	$(CXX) $(CXX_FLAGS) -iquote "$(path_to_root_dir)$(include_dir)" -c -o $@ $*.cc
ifneq ($(depend_files),)
-include $(depend_files)
endif
        
# rules for building subcomponents
.PHONY: subcomponents $(subdirs)
subcomponents: $(subdirs)
$(subdirs): | depend_directories bin_directories subdir_makerules subdir_makefiles
	@$(MAKE) -r -C $@ all

# clean local branch of the tree
.PHONY: clean
clean:
	rm -f $(depend_files)
	rm -f $(object_files)
        
# clean lower branches of the tree
.PHONY: clean_subdirs
clean_subdirs: 
	@for dir in $(subdirs_with_makefiles); do \
		$(MAKE) -r -C $$dir clean_all; \
	done

# clean local branch of the tree and lower
.PHONY: clean_all
clean_all: | clean clean_subdirs
	rm -rf $(depend_subdirs)
	rm -rf $(bin_subdirs)
	rm -f $(makerules)
	rm -f $(makefiles)
