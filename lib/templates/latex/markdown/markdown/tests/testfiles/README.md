This directory contains subdirectories with test files, which can be recognized
by the `.test` suffix. A test file has the following syntax:

    The test setup TeX source code
    <<<
    The test markdown source code
    >>>
    The expected test output

The test setup TeX source code can be used to configure the Markdown package
through its plain TeX interface before the test markdown source code is
processed.

The test markdown source code is the markdown code that will be processed
during the test. The majority of markdown tokens are configured by the support
files to produce output to the log file. This output will be compared against
the expected test output.

The `<<<` and `>>>` markers may be surrounded by optional whitespaces. If the
last section beginning with `>>>` is not present, it will be automatically
generated during the testing and appended to the test file.
