#include "abc_000_warning.h"
#include "abc_000_macro.h"
 

volatile char ctrl_C_Pressed = 0;

#if ( defined(OS_WIN64) || defined (OS_WIN32)) && defined (COMPILER_MSVC)
//#ifdef COMPILER_MSVC

#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#include <signal.h>
#include <stdlib.h>
#include "abc_ide_util.h"



//signal vs signalaction
//windows: SetConsoleCtrlHandler
//https://stackoverflow.com/questions/16826097/equivalent-to-sigint-posix-signal-for-catching-ctrlc-under-windows-mingw
//https://www.tutorialspoint.com/cplusplus/cpp_signal_handling.htm
//https://stackoverflow.com/questions/4217037/catch-ctrl-c-in-c/54267342
//https://www.geeksforgeeks.org/signals-c-language/
//https://stackoverflow.com/questions/17766550/ctrl-c-interrupt-event-handling-in-linux

// Define the function to be called when ctrl-c (SIGINT) is sent to process
static void signal_callback_handler(int signum) {
	ctrl_C_Pressed = 1;
	printf("CTRL+C is pressed!\n");	
}


//https://github.com/rstudio/rstudio/issues/1970
//This isn't working for Rgui/Rstudio because the signal won't be sent directly to the
//the process. 

//https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-December/001502.html
void RegisterCtrlCHandler() {
	_crt_signal_t prev=signal(SIGINT, signal_callback_handler); 
}

#endif

#include "abc_000_warning.h"