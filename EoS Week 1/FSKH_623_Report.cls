%%
%% This is file `FSKH_623_Report.cls', based on 'Thesis.cls', My ``Katlego Moloto'' MSc
%%
%% Created by Katlego D. Moloto,


%% ----------------------------------------------------------------
%% Creating the document class
%% ---------------------------------------------------------------
\ProvidesClass{FSKH_623_Report}[2016/06/14 v2.0 LaTeX document class]
\def\baseclass{book}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\baseclass}}
\def\@checkoptions#1#2{
  \edef\@curroptions{\@ptionlist{\@currname.\@currext}}
  \@tempswafalse
  \@tfor\@this:=#2\do{
    \@expandtwoargs\in@{,\@this,}{,\@curroptions,}
    \ifin@ \@tempswatrue \@break@tfor \fi}
  \let\@this\@empty
  \if@tempswa \else \PassOptionsToClass{#1}{\baseclass}\fi
}
\@checkoptions{11pt}{{10pt}{11pt}{12pt}}
\PassOptionsToClass{a4paper, oneside}{\baseclass}
\ProcessOptions\relax
\LoadClass{\baseclass}

%% ----------------------------------------------------------------
%% Setting up the title page
%% ----------------------------------------------------------------
\newcommand*{\studentnumber}[1]{\def\studentnumber{#1}}
\newcommand*{\studyunit}[1]{\def\studyunit{#1}}
\newcommand*{\problemnumber}[1]{\def\problemnumber{#1}}
\renewcommand{\maketitle}{%
\makeatletter
\let\thetitle\@title
\let\theauthor\@author
\let\thedate\@date
\makeatother
\begin{titlepage}
    \centering
%     \vspace*{0.5 cm}
%     \includegraphics[scale = 0.75]{NWU_LOGO}\\[1.0 cm]   % University Logo
    \includegraphics[clip, trim=0cm 0cm 0.0cm 0cm]{./Figures/nwu-logo.jpg}\\[1.0 cm]   % University Logo
    \textsc{\LARGE Centre For Space Physics}\\[2.0 cm]                                     % Centre Name
    \textsc{\Large Honours}\\[0.5 cm]                                         % Course Number
%    \textsc{\large Study Unit: \studyunit}\\[0.5 cm]                                       % Study Unit Number
    %\textsc{\large \emph{Week:} \problemnumber}\\[0.5 cm]                              % problem number
    \rule{\linewidth}{0.2 mm} \\[0.4 cm]
    { \huge \bfseries \thetitle}\\
    \rule{\linewidth}{0.2 mm} \\[1.5 cm]
    
    \begin{minipage}{0.4\textwidth}
        \begin{flushleft} \large
            \emph{Author:}\\
            \theauthor
            \end{flushleft}
            \end{minipage}~
            \begin{minipage}{0.4\textwidth}
            \begin{flushright} \large
            \emph{Student Number:} \\
            \studentnumber                                   % Your Student Number
        \end{flushright}
    \end{minipage}\\[2 cm]
    
    {\large \thedate}\\[2 cm]
 
    \vfill
    
\end{titlepage}

}

%% ----------------------------------------------------------------
%% Setting up margins and paragraph styles
%% ----------------------------------------------------------------
\usepackage{setspace}
\onehalfspacing
\setlength{\parindent}{0pt}                                        %% sets the new paragraph indent
\setlength{\parskip}{2.0ex plus0.5ex minus0.2ex}                   %% sets the paragraph spacing
\setstretch{1.3}                                                   %% sets the line spacing to 1.3
\usepackage{vmargin}
\setmarginsrb           { 1.3in}  % left margin
                        { 0.6in}  % top margin
                        { 0.8in}  % right margin
                        { 0.8in}  % bottom margin
                        {  12pt}  % head height
                        {0.25in}  % head sep
                        {   9pt}  % foot height
                        { 0.3in}  % foot sep
\raggedbottom
\setlength{\topskip}{1\topskip \@plus 5\p@}
\doublehyphendemerits=10000       % No consecutive line hyphens.
\brokenpenalty=10000              % No broken words across columns/pages.
\widowpenalty=9999                % Almost no widows at bottom of page.
\clubpenalty=9999                 % Almost no orphans at top of page.
\interfootnotelinepenalty=9999    % Almost never break footnotes.

%% ----------------------------------------------------------------
%% Setting up the fancy pages headers and footer
%% ----------------------------------------------------------------
\usepackage{fancyhdr}
\pagestyle{fancy}                                                      %% use fancy page style
\rhead[\fancyplain{}{\sl{\leftmark}}]{\rm\thepage}                     %% put the page number on the right header
\lhead[\fancyplain{}{\sl{\leftmark}}]{\rm\theauthor}                   %% put the author on the left header
\chead{}\lfoot{}\rfoot{}\cfoot{}                                       %% leave all footers empty
\renewcommand{\thesection}{\arabic{section}}                           %% set the section number to arabic

%% ----------------------------------------------------------------
%% Setting hyperlinks for graphs and references
%% ----------------------------------------------------------------
\usepackage[pdfpagemode={UseOutlines},bookmarks=true,bookmarksopen=true,
   bookmarksopenlevel=0,bookmarksnumbered=true,hypertexnames=false,
   colorlinks,linkcolor={cyan},citecolor={red},urlcolor={blue},
   pdfstartview={FitV},unicode,breaklinks=true]{hyperref}
\hypersetup{urlcolor=blue, colorlinks=true}  % Colours hyperlinks in blue, but this can be distracting if there are many links.

%% ----------------------------------------------------------------
%% Loading all other required packages
%% ----------------------------------------------------------------
\usepackage{amsmath,amsfonts,amssymb,amscd,amsthm,mathtools,xspace,siunitx}   % all the math packages
\usepackage{graphicx}                                                         % for the graphs
\usepackage{epstopdf}
\usepackage{calc}
\usepackage{overpic}
\usepackage{multirow}
\usepackage{natbib}                                                           % Use the "Natbib" style for the references in the Bibliography
\usepackage{verbatim}                                                         % Needed for the "comment" environment to make LaTeX comments
\usepackage{booktabs}
\usepackage{rotating}

%% ----------------------------------------------------------------
%% Setting the graph caption 
%% ----------------------------------------------------------------
\usepackage{subcaption}
\pdfoptionpdfminorversion=7
\usepackage[font={small,it},margin=20pt,labelfont=bf]{caption}
\usepackage[Sonny]{fncychap}
\usepackage{lettrine}
\RequirePackage{fix-cm} % to auto correct the missing fontsize warning.

%% ----------------------------------------------------------------
%% Setting programing Code listings 
%% ----------------------------------------------------------------
\usepackage{listings}
\usepackage{lstpatch}
        
\usepackage{color}
 
\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
 
\lstdefinestyle{mystyle}{
    backgroundcolor=\color{backcolour},
    frame=tb,
    commentstyle=\color{codegreen}\em,
    keywordstyle=\color{magenta},
    numberstyle=\tiny\color{codegray},
    stringstyle=\color{codepurple},
    basicstyle=\footnotesize, %% \scriptsize\ttfamily,
    breakatwhitespace=false,         
    breaklines=true,                 
    captionpos=b,                    
    keepspaces=true,                 
    numbers=left,                    
    numbersep=5pt,                  
    showspaces=false,                
    showstringspaces=false,
    showtabs=false,                  
    tabsize=2
}
 
\lstset{style=mystyle}      

%% End of file `FSKH_623_Report.cls'.
