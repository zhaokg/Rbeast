.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  requireNamespace("utils",quitely=TRUE)
  if (interactive())
  {
    v = utils::packageVersion("Rbeast")
    packageStartupMessage("Rbeast v", 
	v, 
	". For help or bug reporting, contact Kaiguang Zhao at zhao.1423@osu.edu. ",
	"Major functions available:\n",
	"  ?beast       - for a single regular time series\n",
	"  ?beast.irreg - for a single irregular time series\n",
	"  ?beast123    - for one, multiple, or 3D array of regular/irregular time series \n",
	"               - (e.g., satellite images): support parallel computing internally\n",
	"  ?minesweeper - a poor man\'s implementation of the minesweeper game\n",
	"  ?tetris      - a poor man\'s implementation of the tetris game (Windows only)\n"
	)
  }
   
   
}

.onLoad <- function(libname, pkgname) {
   #library.dynam("beast", pkgname, libname )
   #utils::data(simdata, package=pkgname,         envir=parent.env(environment())) 
   #utils::data(modis_ohio, package=pkgname,      envir=parent.env(environment())) 
}

.onUnload <- function(libpath) {
  library.dynam.unload("Rbeast", libpath)
}

 