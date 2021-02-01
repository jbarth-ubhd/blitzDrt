# blitzDrt
Tool to correct perspective distortion. Does not correct verticals (column separators, table lines) yet. Uses Blitz++ and ImageMagick++ library.

Compile with `./doMake`. Compiles to `productiveDrt`.

Example: `./productiveDrt --if input.png --of output.png --doTrans`

Help: `./productiveDrt --help`

Tries to detect if text is on page - to prevent deskewing of e. g. architectural plans or hatchings.
