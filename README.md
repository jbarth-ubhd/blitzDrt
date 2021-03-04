# blitzDrt

> Tool to correct perspective distortion

* Does not correct verticals (column separators, table lines) yet.
* Uses Blitz++ and ImageMagick++ library.
* "drt" comes from "discrete radon transformation", the intermediate result when rotating the image by -x° to x° to find optimal deskewing angle for top and bottom of the page.

## Installation

Requires ImageMagick++, libfft3, blitz and boost_program_options.

ImageMagick++, libfftw3 and boost_program_options can be installed with `make deps-ubuntu`.

If blitz is not available in your distro, you can build from source with `make install-blitz` (might require `make install-cmake` first for an up-to-date CMake).

Run `make install` to build and install the tool, it will be available as `productiveDrt` in your `$PATH`.

Uninstall with `make uninstall`

## Usage

Example: `./productiveDrt --if input.png --of output.png --doTrans`

Help: `./productiveDrt --help`

Tries to detect if text is on page - to prevent deskewing of e. g. architectural plans or hatchings.

*Warning:* output files are *not* suitable for archival. Metadata could be lost (imagemagick).
