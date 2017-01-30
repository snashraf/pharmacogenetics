\documentclass{minimal}
% Load the package.
\usepackage{markdown}
% Load the support files.
\input setup\relax
\usepackage{catchfile}
\CatchFileDef{\LATEXOPTIONS}{latex-setup}\relax
% Load the test-specific setup.
\input TEST_SETUP_FILENAME\relax
\begin{document}
% Perform the test.
\expandafter\markdownInput\expandafter[\LATEXOPTIONS]{TEST_INPUT_FILENAME}%
\end{document}
