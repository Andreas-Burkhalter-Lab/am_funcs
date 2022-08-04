REM This is an example of how to build SON for DOS,
REM but it can also be built for small and medium
REM models, and you could play with optimisation ...
REM
cl -W4 -Oselg -AL -c son.c > son.err
type son.err

