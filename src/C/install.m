clc
mex sib_oe_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I. -lmwlapack -lut
mex sib_armax_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I./ -lmwlapack -lut
mex sib_bj_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I./ -lmwlapack -lut