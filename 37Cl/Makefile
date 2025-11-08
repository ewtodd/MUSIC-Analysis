# Check if nm-toolkit is available via pkg-config
ifneq ($(shell pkg-config --exists nm-toolkit && echo yes),yes)
    $(error nm-toolkit not found. Run 'nix develop' first)
endif

TOOLKIT_CFLAGS = $(shell pkg-config --cflags nm-toolkit)
TOOLKIT_LIBDIR = $(shell pkg-config --variable=libdir nm-toolkit)
TOOLKIT_LIBPATH = $(TOOLKIT_LIBDIR)/lib-nm-toolkit.so

# Helper function to run ROOT macros with toolkit loaded
define run_root_macro
	@cd macros && \
	 echo "Loading Nuclear Measurement Toolkit from: $(TOOLKIT_LIBDIR)" && \
	 echo "Executing $(1) macro..." && \
	 root -l -b -q -e 'gSystem->AddIncludePath("$(TOOLKIT_CFLAGS)"); \
	                    gSystem->Load("$(TOOLKIT_LIBPATH)");' \
			    $(2).cc
endef

initial:
	$(call run_root_macro,Si diode initial processing,SiCalibration)

clean:
	@echo "Nothing to clean (macros run interpreted)"

.PHONY: initial clean

