/*
 * MATLAB Compiler: 4.9 (R2008b)
 * Date: Fri Sep 30 12:05:52 2011
 * Arguments: "-B" "macro_default" "-W" "lib:libcpd" "-T" "link:lib"
 * "cpd_register.m" 
 */

#ifndef __libcpd_h
#define __libcpd_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libcpd
#define PUBLIC_libcpd_C_API __global
#else
#define PUBLIC_libcpd_C_API /* No import statement needed. */
#endif

#define LIB_libcpd_C_API PUBLIC_libcpd_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libcpd
#define PUBLIC_libcpd_C_API __declspec(dllexport)
#else
#define PUBLIC_libcpd_C_API __declspec(dllimport)
#endif

#define LIB_libcpd_C_API PUBLIC_libcpd_C_API


#else

#define LIB_libcpd_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libcpd_C_API 
#define LIB_libcpd_C_API /* No special import/export declaration */
#endif

extern LIB_libcpd_C_API 
bool MW_CALL_CONV libcpdInitializeWithHandlers(mclOutputHandlerFcn error_handler,
                                               mclOutputHandlerFcn print_handler);

extern LIB_libcpd_C_API 
bool MW_CALL_CONV libcpdInitialize(void);

extern LIB_libcpd_C_API 
void MW_CALL_CONV libcpdTerminate(void);



extern LIB_libcpd_C_API 
void MW_CALL_CONV libcpdPrintStackTrace(void);


extern LIB_libcpd_C_API 
bool MW_CALL_CONV mlxCpd_register(int nlhs, mxArray *plhs[],
                                  int nrhs, mxArray *prhs[]);

extern LIB_libcpd_C_API 
long MW_CALL_CONV libcpdGetMcrID() ;



extern LIB_libcpd_C_API bool MW_CALL_CONV mlfCpd_register(int nargout
                                                          , mxArray** Transform
                                                          , mxArray** C
                                                          , mxArray** sigma2
                                                          , mxArray* X
                                                          , mxArray* Y
                                                          , mxArray* opt);

#ifdef __cplusplus
}
#endif

#endif
