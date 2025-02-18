#!/bin/bash
set -e
set -x

# build hisat-3n
make -C hisat-3n -j 16

# build htslib
cd htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make -j 16
cd ../

# build samtools
cd samtools
autoheader										# Build config.h.in (this may generate a warning about
															# AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax					# Generate the configure script
./configure --without-curses	# Needed for choosing optional functionality
make -j 16
cd ..

