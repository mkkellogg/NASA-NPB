c CLASS = B
c  
c  
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
c  
        integer problem_size, niter_default
        parameter (problem_size=102, niter_default=400)
        double precision dt_default
        parameter (dt_default = 0.001d0)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='14 Dec 2011')
        character npbversion*5
        parameter (npbversion='3.3.1')
        character cs1*5
        parameter (cs1='pgf90')
        character cs2*6
        parameter (cs2='$(F77)')
        character cs3*6
        parameter (cs3='(none)')
        character cs4*14
        parameter (cs4='-YI,${INCLUDE}')
        character cs5*28
        parameter (cs5='-O2 -ta=nvidia -Minfo -Mcuda')
        character cs6*46
        parameter (cs6='-O2 -ta=nvidia -Minfo -Mcuda -YL,${LD_LIBRA...')
        character cs7*6
        parameter (cs7='randi8')
