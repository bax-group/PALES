# Stand-alone GNU Makefile for PALES 2.1 on current Linux systems.
#
# The small compatibility sources in include/ and src/standalone.c replace the
# historical NMRPipe/DYNAMO development library.  Only a C compiler, GNU make,
# and the system math library are required.

SHELL := /bin/sh
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

PROGRAM   ?= PALES
BUILD_DIR ?= build
PREFIX    ?= /usr/local
bindir    ?= $(PREFIX)/bin

CC      ?= cc
INSTALL ?= install

# PALES is pre-ANSI/prototype-era C. GNU89 preserves its original behavior;
# -fcommon supports tentative globals in the legacy files on GCC 10 and newer.
OPTFLAGS    ?= -O2
DEBUGFLAGS  ?= -g
LEGACYFLAGS ?= -std=gnu89 -fcommon
WARNFLAGS   ?= -Wall -Wextra -Wno-comment -Wno-parentheses \
               -Wno-unused-parameter -Wno-unused-variable \
               -Wno-unused-but-set-variable

CPPFLAGS ?=
CFLAGS   ?= $(OPTFLAGS) $(DEBUGFLAGS)
LDFLAGS  ?=
LDLIBS   ?= -lm

PROJECT_CPPFLAGS := -Isrc -Iinclude

SOURCES := \
	src/pales.c \
	src/ssia.c \
	src/analyzeSaupe.c \
	src/analyzeMatrix3D.c \
	src/analyzeEigen.c \
	src/jacobimz.c \
	src/inv_euler.c \
	src/get_euler.c \
	src/do_nsc.c \
	src/do_nscS.c \
	src/nsc.c \
	src/nmrlibMod.c \
	src/pearsn.c \
	src/betacf.c \
	src/betai.c \
	src/erfcc.c \
	src/gammln.c \
	src/fit.c \
	src/gammq.c \
	src/gcf.c \
	src/gser.c \
	src/gauss.c \
	src/genrand.c \
	src/nrutil.c \
	src/standalone.c \
	include/cmndargs.c \
	include/memory.c \
	include/strmatch.c \
	include/token.c

OBJECTS  := $(SOURCES:%.c=$(BUILD_DIR)/%.o)
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: all pales check clean debug sanitize install print-config help

all: $(PROGRAM)

pales: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(BUILD_DIR)/%.o: %.c
	@mkdir -p "$(dir $@)"
	$(CC) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CFLAGS) $(LEGACYFLAGS) \
		$(WARNFLAGS) -MMD -MP -c $< -o $@

# A deterministic six-vector RDC fit exercises argument parsing, PDB/table
# input, atom mapping, SVD fitting, and table output.
check: all
	@mkdir -p "$(BUILD_DIR)"
	./$(PROGRAM) -bestFit -pdb tests/bestfit.pdb -inD tests/bestfit.tbl \
		-outD "$(BUILD_DIR)/bestfit.out"
	grep -q "DATA SAUPE" "$(BUILD_DIR)/bestfit.out"

debug:
	$(MAKE) clean
	$(MAKE) OPTFLAGS=-O0 DEBUGFLAGS=-g3 all

sanitize:
	$(MAKE) clean
	$(MAKE) CFLAGS='-O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer' \
		LDFLAGS='-fsanitize=address,undefined' all

install: $(PROGRAM)
	$(INSTALL) -d "$(DESTDIR)$(bindir)"
	$(INSTALL) -m 0755 $(PROGRAM) "$(DESTDIR)$(bindir)/$(PROGRAM)"

clean:
	rm -rf -- "$(BUILD_DIR)" "$(PROGRAM)"

print-config:
	@echo "CC=$(CC)"
	@echo "CPPFLAGS=$(CPPFLAGS) $(PROJECT_CPPFLAGS)"
	@echo "CFLAGS=$(CFLAGS) $(LEGACYFLAGS) $(WARNFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "LDLIBS=$(LDLIBS)"

help:
	@echo "Build PALES 2.1 on Linux with GNU make and GCC (or CC=clang)."
	@echo
	@echo "No NMRPipe or DYNAMO installation is required."
	@echo "Targets: all, check, debug, sanitize, install, clean, print-config"
	@echo "Useful overrides: CC, CFLAGS, CPPFLAGS, LDFLAGS, PREFIX, DESTDIR"

-include $(DEPFILES)
