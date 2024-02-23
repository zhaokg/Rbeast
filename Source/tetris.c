#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if (defined(OS_WIN64) || defined(OS_WIN32)) && R_INTERFACE==1
#include "abc_001_config.h"
#include "abc_ide_util.h"

#include "windows.h"
 
#define IDT_TIMER 432432434

  static CRITICAL_SECTION cs;
  static int  isCsSet    = 0;
  static int  isDrawing  = 0;
  static WNDPROC prevWndProc = NULL;
  static HWND    hwnd        = 0;
  static int     isTimerSet    = 0;
  
  static int     counter = 0;
  static SEXP    env     = NULL;
  static LRESULT CALLBACK MyNewWndProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
      switch (uMsg)     {
      case WM_MOUSEMOVE:
          break;
      case WM_LBUTTONDOWN:
          break;
      case WM_DESTROY:
      {
          if (hwnd) {
              WNDPROC curWndProc = GetWindowLongPtr(hwnd, -4);
              if (curWndProc == MyNewWndProc && prevWndProc != NULL) {
                  WNDPROC tmp = SetWindowLongPtr(hwnd, (-4), prevWndProc);
                  #if !(defined(R_RELEASE) || defined(M_RELEASE))
                  r_printf(" Restore winproc: %p  %p \n", tmp, MyNewWndProc);
                  #endif
              }
              if (isTimerSet == 1) {
                  KillTimer(hwnd, IDT_TIMER);
                  isTimerSet = 0;
              }
          }
          if (isCsSet) {
              LeaveCriticalSection(&cs);
              DeleteCriticalSection(&cs);
              isCsSet = 0;
          }
          prevWndProc = NULL;
          hwnd        = 0;
          isTimerSet  = 0;
          env         = NULL;
          break;
      }

      case WM_TIMER: {
        
          if (env) {   
              if (isDrawing) { break; }
              isDrawing = 1;
              SEXP e, tmp;
              int  status;
              PROTECT(e = lang2(install("timer"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("prob"));  
              R_tryEval(e, env, NULL); //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
          break;
      }
      case WM_KEYDOWN:

          switch (wParam)
          {
              SEXP e, tmp;
              int  status;
          case VK_LEFT: {
              if (isDrawing) { break; }
              if (!env) { break; }
              isDrawing = 1;
              PROTECT(e = lang2(install("keyleft"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("key"));  
              R_tryEval(e, env, NULL);  //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
               break;
          case VK_RIGHT: {
              if (isDrawing) { break; }
              if (!env) { break; }
              isDrawing = 1;
              PROTECT(e = lang2(install("keyright"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("key"));  
              R_tryEval(e, env, NULL); //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
              break;
          case VK_DOWN: {
              if (isDrawing) { break; }
              if (!env) { break; }
              isDrawing = 1;
              PROTECT(e = lang2(install("keydwon"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("key"));  
              R_tryEval(e, env, NULL); //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
               break;
          case VK_UP: {
              if (isDrawing) { break; }
              if (!env) { break; }
              isDrawing = 1;
              PROTECT(e = lang2(install("keyup"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("key"));  
              R_tryEval(e, env, NULL); //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
               break;
          case VK_SPACE: {
              if (isDrawing) { break; }
              if (!env) { break; }
              isDrawing = 1;
              PROTECT(e = lang2(install("keyspace"), ScalarReal(0.6)));
              SET_TAG(CDR(e), install("key"));  
              R_tryEval(e, env, NULL); //R_GlobalEnv         
              UNPROTECT(1);
              isDrawing = 0;
          }
               break;

          }
          break;
      }

      return CallWindowProc(prevWndProc, hwnd, uMsg, wParam, lParam);
  }

  SEXP DllExport TetrisSetTimer(SEXP action, SEXP seconds, SEXP envior) {
      int  x = asInteger(action);
      int  msec = asReal(seconds) * 1000.;
      env = envior;
      if (x == 10) {
          EnterCriticalSection(&cs);
          return R_NilValue;
      }
      if (x == 20) {
          LeaveCriticalSection(&cs);
          return R_NilValue;
      }

      //https://stackoverflow.com/questions/2463437/r-from-c-simplest-possible-helloworld
      //https://stackoverflow.com/questions/7457635/calling-r-function-from-c

      hwnd = FindWindowA(NULL, "Rtetris-Rbeast");
      if (hwnd) {
          WNDPROC curProc = GetWindowLongPtr(hwnd, -4);
          //Rprintf(" curProc %d  \\n ",curProc != MyNewWndProc);
          if (prevWndProc == NULL || curProc != MyNewWndProc) {
              //https://stackoverflow.com/questions/31648180/c-changing-hwnd-window-procedure-in-runtime
              prevWndProc = SetWindowLongPtr(hwnd, (-4), &MyNewWndProc);
          }

          if (x == 1) {
              if (isTimerSet == 0) {
                  SetTimer(hwnd, IDT_TIMER, msec, NULL);
                  isTimerSet = 1;
              }
              return R_NilValue;
          }

          if (x == 0) {
              if (isTimerSet == 1) {
                  KillTimer(hwnd, IDT_TIMER);
                  isTimerSet = 0;
              }
              return R_NilValue;
          }

      }

      if (isCsSet == 0) {
          InitializeCriticalSection(&cs);
          isCsSet = 1;
      }
      return R_NilValue;
  }
#else
static char achar='c';
#endif

#include "abc_000_warning.h"