clc
mex sib_oe_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I. -lmwlapack -lut -lmat
mex sib_armax_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I./ -lmwlapack -lut -lmat
mex sib_bj_c.c sib_basic.c sib_optimize.c CFLAGS='\$CFLAGS -std=c99'  -I./ -lmwlapack -lut -lmat