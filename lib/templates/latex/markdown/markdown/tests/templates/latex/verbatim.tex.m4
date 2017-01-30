\documentclass{minimal}
% Load the package.
\usepackage{markdown}
% Load the support files.
\input setup\relax
\usepackage{catchfile}
\CatchFileDef{\LATEXOPTIONS}{latex-setup}\relax
\expandafter\markdownSetup\expandafter{\LATEXOPTIONS}
% Load the test-specific setup.
\input TEST_SETUP_FILENAME\relax
\begin{document}
% Perform the test.
\begin{markdown}
undivert(TEST_INPUT_FILENAME)dnl
\end{markdown}
\end{document}
