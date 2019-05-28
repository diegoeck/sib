clc
mex sib_oe_c.c sib_basic.c sib_optimize.c   -I. -lmwlapack -lut -lm
mex sib_armax_c.c sib_basic.c sib_optimize.c   -I./ -lmwlapack -lut -lm
mex sib_bj_c.c sib_basic.c sib_optimize.c   -I./ -lmwlapack -lut -lm