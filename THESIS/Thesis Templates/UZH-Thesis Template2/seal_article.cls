% SEAL ARTICLE DOCUMENT CLASS -- version 1.0 (11-Feb-2005)
% Beat Fluri LaTeX2e support for Special Articles
%
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{seal_article}[2005/02/11 v1.0
^^J LaTeX document class for seal style articles]

\LoadClass[10pt,a4paper]{article}
\RequirePackage{graphicx}
\RequirePackage{listings}

\renewcommand{\rmdefault}{ppl}
\renewcommand{\sfdefault}{phv}
\renewcommand{\ttdefault}{pcr}

% Ragged bottom for the actual page (from LNCS)
\def\thisbottomragged{\def\@textbottom{\vskip\z@ plus.0001fil
\global\let\@textbottom\relax}}

\setlength{\textwidth}{15cm}
\setlength{\textheight}{21.3cm}

\setlength\oddsidemargin   {11\p@}
\setlength\evensidemargin  {11\p@}
\setlength\marginparwidth  {90\p@}

\setlength\headsep   {16\p@}

\setlength\footnotesep{7.7\p@}
\setlength\textfloatsep{8mm\@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {8mm\@plus 2\p@ \@minus 2\p@}

\setcounter{secnumdepth}{2}

\def\ps@sealheadings{\let\@mkboth\@gobbletwo
   \def\@oddfoot{\hfil\fontfamily{pfr}\selectfont\small\bfseries\thepage}
   \def\@evenfoot{\hfil\fontfamily{pfr}\selectfont\small\bfseries\thepage}
   \let\@evenhead\@empty\let\@oddhead\@empty}

\ps@sealheadings

\renewcommand\section{\@startsection{section}{1}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {12\p@ \@plus 4\p@ \@minus 4\p@}%
                       {\fontfamily{pfr}\selectfont\Large\bfseries\boldmath
                        \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {8\p@ \@plus 4\p@ \@minus 4\p@}%
                       {\fontfamily{pfr}\selectfont\bfseries\boldmath
                        \rightskip=\z@ \@plus 8em\pretolerance=10000 }}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                       {-18\p@ \@plus -4\p@ \@minus -4\p@}%
                       {6\p@ \@plus 4\p@ \@minus 4\p@}%
                       {\fontfamily{pfr}\selectfont\fontsize{9}{0}\bfseries\boldmath}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                       {-12\p@ \@plus -4\p@ \@minus -4\p@}%
                       {-0.5em \@plus -0.22em \@minus -0.1em}%
                       {\fontfamily{pfr}\selectfont\itshape}}
\renewcommand\subparagraph[1]{\typeout{LLNCS warning: You should not use
                  \string\subparagraph\space with this class}\vskip0.5cm
You should not use \verb|\subparagraph| with this class.\vskip0.5cm}

% list, itemize, and enumerate margin definitions (from LNCS)
\setlength\leftmargini  {17\p@}
\setlength\leftmargin    {\leftmargini}
\setlength\leftmarginii  {\leftmargini}
\setlength\leftmarginiii {\leftmargini}
\setlength\leftmarginiv  {\leftmargini}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}

\def\@listI{\leftmargin\leftmargini
            \parsep 0\p@ \@plus1\p@ \@minus\p@
            \topsep 8\p@ \@plus2\p@ \@minus4\p@
            \itemsep0\p@}
\let\@listi\@listI
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    0\p@ \@plus2\p@ \@minus\p@}
\def\@listiii{\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep
              \topsep    0\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@}

% smaller dot for first item in itemize
\renewcommand\labelitemi{\tiny\raisebox{1.5pt}{$\m@th\bullet$}}

%\def\fnum@figure{{\bf\fontfamily{pfr}\selectfont \figurename\nobreakspace\thefigure}}
%\def\fnum@table{{\fontfamily{pfr}\selectfont\bf \tablename\nobreakspace\thetable}}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\fontfamily{pfr}\selectfont{#1}: \fontfamily{pfr}\selectfont{#2}}%
  \ifdim \wd\@tempboxa >\hsize
    \fontfamily{pfr}\selectfont{#1}: \fontfamily{pfr}\selectfont{#2}\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\box\@tempboxa}%
  \fi
  \vskip\belowcaptionskip}

\def\subtitle#1{\gdef\@subtitle{#1}}
\def\@subtitle{\@latex@error{No \noexpand\subtitle given}\@ehc}
\renewcommand\maketitle{%
  \newpage
  \vspace*{-3cm}
  \begingroup
  
  \parindent \z@
  \includegraphics[width=7cm]{logo_inf}\hfill\includegraphics[width=3.5cm]{seal_bw}
  \vskip 48\p@
  \vskip 3cm
  \begin{flushright}
  {\fontfamily{pfr}\selectfont\bfseries\Huge\@title}
  \end{flushright}
  \vskip -10\p@
  \vbox{\hrule height0.25pt width\textwidth}
  \begin{flushright}
  \vskip -8\p@
  {\fontfamily{pfr}\selectfont\Large\@subtitle}\par\nobreak
  \end{flushright}
  \vskip 40\p@
  \endgroup
}
%
\def\@maketitle{\newpage
 \let\newline\\
  \begingroup
  \parindent \z@
  \includegraphics[width=7cm]{logo_inf}\hfill\includegraphics[width=3.5cm]{seal_bw}
  \vskip 48\p@
  \vskip 3cm
  \begin{flushright}
 {\fontfamily{pfr}\selectfont\bfseries\Huge\@title}
  \end{flushright}
  \vskip -10\p@
  \vbox{\hrule height0.25pt width\textwidth}
  \begin{flushright}
  \vskip -8\p@
  {\fontfamily{pfr}\selectfont\Large\@subtitle}\par\nobreak
  \end{flushright}
  \vskip 40\p@
  \endgroup
}

\endinput

%end of seal_article.cls