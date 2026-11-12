PROGRAM := EtaleUpperBound
SOURCE := EtaleUpperBound.c
BUILD_DIR := build

CC ?= cc
CPPFLAGS ?=
CFLAGS ?= -O3 -DNDEBUG -flto -std=c11 -Wall -Wextra -Wpedantic
LDFLAGS ?= -flto
LDLIBS ?=

PREFIX ?= /usr/local
DESTDIR ?=

# Set PARI_PREFIX when PARI/GP is installed outside the compiler's default
# search path, for example PARI_PREFIX="$(shell brew --prefix pari)" on macOS.
PARI_PREFIX ?=
PKG_CONFIG ?= pkg-config
PARI_PKG_CFLAGS := $(shell $(PKG_CONFIG) --cflags pari 2>/dev/null)
PARI_PKG_LIBS := $(shell $(PKG_CONFIG) --libs pari 2>/dev/null)

ifneq ($(strip $(PARI_PREFIX)),)
CPPFLAGS += -I$(PARI_PREFIX)/include
LDFLAGS += -L$(PARI_PREFIX)/lib
endif
CPPFLAGS += $(PARI_PKG_CFLAGS)
ifeq ($(strip $(PARI_PKG_LIBS)),)
LDLIBS += -lpari
else
LDLIBS += $(PARI_PKG_LIBS)
endif
LDLIBS += -lm

# By default, build the Cliquer sources placed in vendor/cliquer.  On systems
# providing libcliquer and cliquer.h, use: make USE_SYSTEM_CLIQUER=1
USE_SYSTEM_CLIQUER ?= 0
CLIQUER_DIR ?= vendor/cliquer
CLIQUER_VERSION := 1.21
CLIQUER_URL := https://users.aalto.fi/~pat/cliquer/cliquer-$(CLIQUER_VERSION).tar.gz

PROJECT_OBJECTS := $(BUILD_DIR)/EtaleUpperBound.o \
                   $(BUILD_DIR)/etale_bitset.o

ifeq ($(USE_SYSTEM_CLIQUER),1)
OBJECTS := $(PROJECT_OBJECTS)
LDLIBS += -lcliquer
else
CPPFLAGS += -I$(CLIQUER_DIR)
CLIQUER_OBJECTS := $(BUILD_DIR)/cliquer.o \
                   $(BUILD_DIR)/graph.o \
                   $(BUILD_DIR)/reorder.o
OBJECTS := $(PROJECT_OBJECTS) $(CLIQUER_OBJECTS)
endif

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(BUILD_DIR):
	mkdir -p $@

ifeq ($(USE_SYSTEM_CLIQUER),1)
$(BUILD_DIR)/EtaleUpperBound.o: $(SOURCE) etale_bitset.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
else
$(BUILD_DIR)/EtaleUpperBound.o: $(SOURCE) etale_bitset.h | $(BUILD_DIR) prepare-cliquer
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
endif

$(BUILD_DIR)/etale_bitset.o: etale_bitset.c etale_bitset.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

ifeq ($(USE_SYSTEM_CLIQUER),0)
# This order-only prerequisite is deliberately shared by the main source and
# all Cliquer objects.  It prevents a parallel or clean build from trying to
# include cliquer.h before the vendored sources have been downloaded.
.PHONY: prepare-cliquer
prepare-cliquer: | $(BUILD_DIR)
	@if ! test -f "$(CLIQUER_DIR)/cliquer.c" || \
	    ! test -f "$(CLIQUER_DIR)/graph.c" || \
	    ! test -f "$(CLIQUER_DIR)/reorder.c" || \
	    ! test -f "$(CLIQUER_DIR)/cliquer.h"; then \
		echo "Cliquer sources not found; downloading Cliquer $(CLIQUER_VERSION)..."; \
		$(MAKE) --no-print-directory fetch-cliquer; \
	fi

$(BUILD_DIR)/cliquer.o: | $(BUILD_DIR) prepare-cliquer
	$(CC) $(CPPFLAGS) $(CFLAGS) -c "$(CLIQUER_DIR)/cliquer.c" -o $@

$(BUILD_DIR)/graph.o: | $(BUILD_DIR) prepare-cliquer
	$(CC) $(CPPFLAGS) $(CFLAGS) -c "$(CLIQUER_DIR)/graph.c" -o $@

$(BUILD_DIR)/reorder.o: | $(BUILD_DIR) prepare-cliquer
	$(CC) $(CPPFLAGS) $(CFLAGS) -c "$(CLIQUER_DIR)/reorder.c" -o $@
endif

.PHONY: fetch-cliquer
fetch-cliquer:
	@set -eu; \
	mkdir -p vendor $(BUILD_DIR); \
	archive="$(BUILD_DIR)/cliquer-$(CLIQUER_VERSION).tar.gz"; \
	if command -v curl >/dev/null 2>&1; then \
		curl -L --fail "$(CLIQUER_URL)" -o "$$archive"; \
	elif command -v wget >/dev/null 2>&1; then \
		wget -O "$$archive" "$(CLIQUER_URL)"; \
	else \
		echo "Neither curl nor wget is available." >&2; \
		exit 1; \
	fi; \
	rm -rf "$(CLIQUER_DIR)" "vendor/cliquer-$(CLIQUER_VERSION)"; \
	tar -xzf "$$archive" -C vendor; \
	mv "vendor/cliquer-$(CLIQUER_VERSION)" "$(CLIQUER_DIR)"; \
	echo "Cliquer $(CLIQUER_VERSION) installed in $(CLIQUER_DIR)."

.PHONY: test
test: $(PROGRAM)
	@set -eu; \
	tmp=$$(mktemp -d); \
	trap 'rm -rf "$$tmp"' EXIT INT TERM; \
	(cd "$$tmp" && "$(abspath $(PROGRAM))" 3 1 -cpus 1); \
	grep -q '^Orbit computations: 1/1 successfully completed$$' "$$tmp/p=3_d=1/result.out"; \
	grep -q '^Clique number: 3$$' "$$tmp/p=3_d=1/result.out"; \
	(cd "$$tmp" && "$(abspath $(PROGRAM))" 7 1 1 1 -cpus 2); \
	grep -q '^Orbit computations: 3/3 successfully completed$$' "$$tmp/p=7_d=1x1x1/result.out"; \
	for file in "$$tmp"/p=7_d=1x1x1/local/local_*.out; do \
		grep -q '^Best: ' "$$file"; \
		grep -q '^Status: complete$$' "$$file"; \
	done; \
	echo "Smoke tests passed."

.PHONY: debug
debug:
	$(MAKE) clean
	$(MAKE) CFLAGS='-O0 -g3 -std=c11 -Wall -Wextra -Wpedantic' LDFLAGS=''

.PHONY: sanitize
sanitize:
	$(MAKE) clean
	$(MAKE) \
	  CFLAGS='-O1 -g3 -std=c11 -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer' \
	  LDFLAGS='-fsanitize=address,undefined'

.PHONY: install
install: $(PROGRAM)
	mkdir -p "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 $(PROGRAM) "$(DESTDIR)$(PREFIX)/bin/$(PROGRAM)"

.PHONY: uninstall
uninstall:
	rm -f "$(DESTDIR)$(PREFIX)/bin/$(PROGRAM)"

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(PROGRAM)

.PHONY: help
help:
	@printf '%s\n' \
	  'make                     Download Cliquer if needed, then build' \
	  'make fetch-cliquer       Force a fresh Cliquer 1.21 download' \
	  'make USE_SYSTEM_CLIQUER=1  Link against an installed libcliquer' \
	  'make test                Run sequential and parallel smoke tests' \
	  'make debug               Build an unoptimized debug binary' \
	  'make sanitize            Build with ASan and UBSan' \
	  'make clean               Remove generated build files'

.PHONY: all
