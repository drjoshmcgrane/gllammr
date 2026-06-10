# Custom Makefile for GLLAMMR
# Uses pre-compiled individual TMB shared libraries instead of linking into one

# The .so files are already compiled via TMB::compile()
# Just need to ensure they're available

all: $(SHLIB)

# Override the default rule - don't try to link all .o files together
# Instead, just touch a dummy file to satisfy R's build system
$(SHLIB):
	@echo "Using pre-compiled TMB shared libraries"
	@touch $(SHLIB)

clean:
	rm -f *.o *.so

.PHONY: all clean
