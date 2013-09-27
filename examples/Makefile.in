## Allows for making of all examples with single command
## Assumes example directories contain Makefile
## and are located two levels below

exampledirs = $(dir $(wildcard */*/Makefile))

all:
##      Stops on first error
	@$(foreach dir,$(exampledirs),$(MAKE) -C $(dir) all &&) echo ""

clean:
	@$(foreach dir,$(exampledirs),$(MAKE) -C $(dir) clean;)

print:
	@echo $(exampledirs) | tr " " "\n"
