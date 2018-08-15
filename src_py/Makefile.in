#
# Makefile for generating uvspec_lex.l from src_py dictionary
#

TOPDIR = ..

include $(TOPDIR)/Makeconf

PYTHONFILES = writeLex.py.in \
                GUI_definition.py \
		lex_starter.l \
		option_definition.py \
		cloud_options.py \
		profile_options.py \
		aerosol_options.py \
		mc_options.py \
		spectral_options.py \
		geometry_options.py \
		general_atmosphere_options.py \
		molecular_options.py \
		surface_options.py \
		special_options.py \
		solver_options.py \
		output_options.py \
                translate.py

DISTFILES = Makefile.in $(PYTHONFILES) uvspec_lex.l writeLex.py

MYSTICDISTFILES = $(DISTFILES)

all: uvspec_lex.l
.PHONY: all	

uvspec_lex.l: $(PYTHONFILES) writeLex.py
	$(PYTHON) writeLex.py

clean:
.PHONY: clean

distclean: clean
	rm -f Makefile
.PHONY: distclean

realclean: distclean
.PHONY: realclean

install: 
.PHONY: install

uninstall:
.PHONY: uninstall

dist:
	ln $(DISTFILES) ../`cat ../.fname`/src_py
.PHONY: dist

mystic-dist:
	ln $(MYSTICDISTFILES) ../`cat ../.fname`/src_py
.PHONY: mystic-dist

Makefile: Makefile.in ../config.status
	cd ..; sh config.status

writeLex.py: writeLex.py.in
	cd ..; sh config.status
