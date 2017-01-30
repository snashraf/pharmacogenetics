% Load the package.
\input markdown\relax
% Load the support files.
\input setup\relax
\input plain-setup\relax
% Load the test-specific setup.
\input TEST_SETUP_FILENAME\relax
% Perform the test.
\markdownBegin
undivert(TEST_INPUT_FILENAME)dnl
\markdownEnd
\bye
