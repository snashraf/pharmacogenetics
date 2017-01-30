This directory contains subdirectories, which correspond to individual TeX
formats. These subdirectories contain TeX source code templates, which
can be recognized by the `.tex.m4` suffix, and a file named `COMMANDS.m4`.

The TeX source code templates are documents that will be typeset as a part of
every test. Before being typeset, the TeX source code templates will be
pre-processed using the m4 macro preprocessor with the following macro
definitions:

 1. `TEST_SETUP_FILENAME` corresponds to the name of the file containing the
    test setup TeX source code.
 2. `TEST_INPUT_FILENAME` corresponds to the name of the file containing the
    test markdown source code.

The `COMMANDS` file contains a newline-separated list of commands that will
be used to typeset the pre-processed TeX source code templates. Before use,
the commands will preprocessed with the following macro definitions:

 1. `TEST_FILENAME` corresponds to the name of the pre-processed TeX source
    code template that is being typeset.
