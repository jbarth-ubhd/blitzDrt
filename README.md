# blitzDrt
Tool to correct perspective distortion

* Does not correct verticals (column separators, table lines) yet.
* Uses Blitz++ and ImageMagick++ library.
* "drt" comes from "discrete radon transformation", the intermediate result when rotating the image by -x° to x° to find optimal deskewing angle for top and bottom of the page.

Compile with `./doMake`. Compiles to `productiveDrt`.

Example: `./productiveDrt --if input.png --of output.png --doTrans`

Help: `./productiveDrt --help`

Tries to detect if text is on page - to prevent deskewing of e. g. architectural plans or hatchings.
