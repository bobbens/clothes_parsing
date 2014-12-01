include Rules.make

.PHONY: default
default: all
 
.PHONY: all
all: src test

.PHONY: world
world: all check_world

.PHONY: src
src:
	@$(MAKE) -r -C $(src_dir) all

.PHONY: test
test: src
	@$(MAKE) -r -C $(test_dir) all

.PHONY: test_results
test_results:
	@$(MAKE) -r -C $(test_dir) test_results_all

.PHONY: test_results_world
test_results_world: test
	@$(MAKE) -r -C $(test_dir) test_results_all

.PHONY: check
check:
	@$(MAKE) -r -C $(test_dir) check_all

.PHONY: check_world
check_world: test_results_world
	@$(MAKE) -r -C $(test_dir) check_all

.PHONY: matlab
matlab:
	@$(MAKE) -r -C $(matlab_dir) all

.PHONY: doc
doc:
	@$(MAKE) -r -C $(doc_dir) all

.PHONY: doc_clean
doc_clean:
	@$(MAKE) -r -C $(doc_dir) clean_all

.PHONY: clean
clean:
	@$(MAKE) -r -C $(src_dir) clean_all
	@$(MAKE) -r -C $(test_dir) clean_all
	@$(MAKE) -r -C $(matlab_dir) clean_all
