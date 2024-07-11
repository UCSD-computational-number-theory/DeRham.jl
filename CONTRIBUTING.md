# Contribution Guides

This file contains instuctions on how to do some things in Julia, 
since contributors might not be familiar.

### How to use this package on the REPL while you are developing it

The simplest way is to `cd src` and run
`include("DeRham.jl")` in the Julia REPL.
Then you can access all the functions and methods in the package
with the prefix `DeRham.[----]` where `[----]` is the function
you want to call.

When you modify a file, you can re-run `include("DeRham.jl")`.

This isn't the most elegant, but it gets the job done.

### Running tests

1. Make sure you're in the top-level project folder.
  * you can do this by running `julia` in the project folder
  * or you can type `;` to enter julia's shell mode and use `cd` to get 
    into the project folder. You can use backspace to get out of shell mode.
2. Enter package mode by typing `]`
3. `activate .` (still in package mode)
4. `test` (still in package mode)

You can also try `include("runtests.jl")` while in the test folder, 
but this is not what the continuous integration will do.

### Adding a dependency

This is a little bit convoluted, since we actually have to add each dependency
both to the usual project and to the tests. 
Hopefully soon Julia will have a better solution, see
[here](https://discourse.julialang.org/t/adding-test-specific-dependencies/101100/4).

1. Enter package mode with `]`. While in package mode:
2. `activate .`
3. `add [packagename]` where `[packagename]` the name and/or github url of the package
4. `actiavte test`
5. `add [packagename]` (for a second time)
6. Now that your package is added in the package manager, you can go to `DeRham.jl` and
   add `using [packagename]` at the top
7. Now you're ready to use your package in any of the source/test files!

