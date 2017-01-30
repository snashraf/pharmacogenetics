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
\let\ea\expandafter
\ea\begin\ea{\ea m\ea a\ea r\ea k\ea d\ea o\ea w\ea n\ea *\ea}\ea{\LATEXOPTIONS}
undivert(TEST_INPUT_FILENAME)dnl
\end{markdown*}
\end{document}
