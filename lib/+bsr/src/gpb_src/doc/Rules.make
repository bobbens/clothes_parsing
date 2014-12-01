path_to_root_dir := ../$(path_to_root_dir)
path_from_root_dir := $(call dirlist_last, $(path_to_working_dir))/$(path_from_root_dir)
path_to_working_dir := $(dir $(subst //, , $(path_to_working_dir)/))
include $(path_to_root_dir)Rules.make
