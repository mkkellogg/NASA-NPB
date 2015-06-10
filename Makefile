SHELL=/bin/sh
BENCHMARK=sp
BENCHMARKU=SP

include ../config/make.def


OBJS = sp.o initialize.o exact_solution.o exact_rhs.o \
       set_constants.o adi.o rhs.o      \
       x_solve.o ninvr.o y_solve.o pinvr.o    \
       z_solve.o tzetar.o add.o txinvr.o error.o verify.o glob.o \
       glob_setup.o rhs_kernels.f txinvr_kernels.o lhs_kernels.o \
	x_solve_kernels.o ninvr_kernels.o y_solve_kernels.o pinvr_kernels.o \
	z_solve_kernels.o tzetar_kernels.o add_kernels.o util_kernels.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o ${COMMON}/wtime.o

include ../sys/make.common

# npbparams.h is included by header.h
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to 
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect). 

# header.h: npbparams.h

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

util_kernels.o:    util_kernels.f header.h npbparams.h
add_kernels.o:     add_kernels.f header.h npbparams.h
tzetar_kernels.o:  tzetar_kernels.f header.h npbparams.h
z_solve_kernels.o: z_solve_kernels.f header.h npbparams.h
pinvr_kernels.o:   pinvr_kernels.f header.h npbparams.h
y_solve_kernels.o: y_solve_kernels.f header.h npbparams.h
ninvr_kernels.o:   ninvr_kernels.f header.h npbparams.h
x_solve_kernels.o: x_solve_kernels.f header.h npbparams.h
lhs_kernels.o:     lhs_kernels.f header.h npbparams.h
txinvr_kernels.o:  txinvr_kernels.f header.h npbparams.h
rhs_kernels.o:     rhs_kernels.f header.h npbparams.h
glob_setup.o:      glob_setup.f header.h npbparams.h
glob.o:            glob.f header.h npbparams.h
sp.o:              glob.o glob_setup.o sp.f  header.h npbparams.h
initialize.o:      util_kernels.o lhs_kernels.o initialize.f  header.h npbparams.h
exact_solution.o:  exact_solution.f  header.h npbparams.h
exact_rhs.o:       exact_rhs.f  header.h npbparams.h
set_constants.o:   set_constants.f  header.h npbparams.h
adi.o:             adi.f  header.h npbparams.h
rhs.o:             rhs_kernels.o rhs.f  header.h npbparams.h
#lhsx.o:           lhsx.f  header.h npbparams.h
#lhsy.o:           lhsy.f  header.h npbparams.h
#lhsz.o:           lhsz.f  header.h npbparams.h
x_solve.o:         x_solve_kernels.o ninvr.o lhs_kernels.o util_kernels.o x_solve.f  header.h npbparams.h
ninvr.o:           ninvr_kernels.o ninvr.f  header.h npbparams.h
y_solve.o:         util_kernels.o y_solve_kernels.o lhs_kernels.o y_solve.f  header.h npbparams.h
pinvr.o:           pinvr_kernels.o pinvr.f  header.h npbparams.h
z_solve.o:         util_kernels.o z_solve_kernels.o z_solve.f  header.h npbparams.h
tzetar.o:          tzetar_kernels.o tzetar.f  header.h npbparams.h
add.o:             util_kernels.o add_kernels.o add.f  header.h npbparams.h
txinvr.o:          txinvr_kernels.o txinvr.f  header.h npbparams.h
error.o:           error.f  header.h npbparams.h
verify.o:          util_kernels.o verify.f  header.h npbparams.h

clean:
	- rm -f *.o *~ mputil*
	- rm -f npbparams.h core
