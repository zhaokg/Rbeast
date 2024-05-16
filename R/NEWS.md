# CHANGES IN Rbeast 1.0.0

* Fix the warning message regarding "format ‘%x’ expects argument of type ‘unsigned int’, but argument 2 has type ‘VOID_PTR’ {aka ‘void *’} [-Wformat=]"

* Fix the note message regarding checkRd: (-1) beast.Rd:155: Lost braces; missing escapes or markup?

* Add the tseg.leftmargin and tseg.rightmargin arguments upon the request of  Racine Nassau


# CHANGES IN Rbeast 0.9.8/9

* Fix the warning message regarding "long long unsigned int[0]' [-Warray-bounds]" and the Wstrict-prototypes warning.


# CHANGES IN Rbeast 0.9.7

* Fix the C23 and clang 16 issue upon the request of Dr. Brian D. Ripley (Jan 20,2023):
  feenableexcept not defined for clang16 with -std=gnu2x. No implicit function declarations
  are allowed. 'false','true', and 'bool' are keywords in C23.

# CHANGES IN Rbeast 0.9.6

* Fix the C23 issue upon the request of Dr. Brian D. Ripley (Jan 09,2023)
* Also remove the freq arg in favor of period to specify the period of the
 the seasonal component (Jan 09,2023)


# CHANGES IN Rbeast 0.9.5

* Fix multiple bugs and add the citation info, thanks to many users.
  Special thanks go to Dr. Zhanmang Liao. (Aug 09,2022)

# CHANGES IN Rbeast 0.9.4

* Fix multiple bugs and add new features, thanks to many users. Special 
  thanks go to Wenpeng Zhao, Josué M. POLANCO-MARTINEZ, Adam Canning,and 
  Shelby McNeill, among others. (May 15,2022)
  
# CHANGES IN Rbeast 0.9.3  

* Mar 19,2021:  Fix multiple bugs for edge cases
* Jan 19,2021:  Fix multiple bugs for edge cases

# CHANGES IN Rbeast 0.9.2  

* Dec 19,2021:  Fix multiple donttest errors reported at CRAN

# CHANGES IN Rbeast 0.9.1  

* Nov 22,2021: Fix the compiler error for MacOS-ARM64

# CHANGES IN Rbeast 0.9.0  

* Nov 11,2021: A new api interface is provided.

# CHANGES IN Rbeast 0.2.2 
 
* Nov 19,2019: Add two contributors for the use of LinPACK as per the CRAN policies.

# CHANGES IN Rbeast 0.2.1  

* Jul 26,2019: Bug fixes: Fixed the compilation errors in the OS X Linux system.
* Jul 20,2019: Bug fixes: Fixed the compilation errors in the Solaris Linux system.

# CHANGES IN Rbeast 0.1

* May 21,2019: Rbeast v0.1 released!

