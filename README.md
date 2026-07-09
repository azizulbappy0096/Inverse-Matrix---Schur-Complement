# Schur Complement - Matrix Inverse

## Overview

This is a sequential implementation of the Schur Complement algorithm for computing matrix inverses.

## Implementation Status

- **Sequential Program**: This is a sequential implementation at this phase

## Supported Matrix Sizes

Currently, the algorithm only works with square matrices of power-of-2 dimensions:

- 2x2
- 4x4
- 8x8
- 16x16
- 32x32
- .....
- 512x512
- And larger power-of-2 sizes

## How To Run

### Using Code::Blocks

1. Open the project file `schur_complement.cbp` in Code::Blocks.
2. Make sure your compiler toolchain is configured in Code::Blocks.
3. Select the `Debug` or `Release` target from the build dropdown.
4. Build and Run the program with `Build and Run` or `Ctrl+F10`.

The executable is a console application, so it will open in a terminal window and print:

- the inverse matrix
- a verification message showing whether `A * A^-1 = I`
