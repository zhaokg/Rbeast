#include "abc_001_config.h"
#include "abc_datatype.h"
#include "abc_common.h"
#include "beast_common.h"
extern int check_options(Options * _restrict opt,char *missing);
#if R_INTERFACE==1
extern int R_check_input(SEXP Y,SEXP opt);
extern int R_read_input(Options * _restrict pOpt,SEXP Y,SEXP opt,char *missing);
extern SEXP R_allocate_output(RESULT * _restrict matOutput,Options * _restrict opt);
#elif M_INTERFACE==1
extern int     M_check_input(int nrhs,const mxArray * _restrict prhs[]);
extern int     M_read_input(Options * _restrict pOpt,int nrhs,const mxArray * _restrict prhs[],char *missing);
extern mxArray *M_allocate_output(RESULT * _restrict matOutput,Options * _restrict opt);
#endif
