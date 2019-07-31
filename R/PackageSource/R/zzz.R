.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  requireNamespace("utils",quitely=TRUE)
  if (interactive())
  {
    v = utils::packageVersion("Rbeast")
    packageStartupMessage("Rbeast v", v, ". For help or bug reporting, please contact Kaiguang Zhao at lidar.rs@gmail.com.\n")
  }
   
   
}

.onLoad <- function(libname, pkgname) {
   #library.dynam("beast", pkgname, libname )
   utils::data(simdata, package=pkgname,         envir=parent.env(environment())) 
   utils::data(modis_ohio, package=pkgname,      envir=parent.env(environment())) 
   #utils::data(simAnnualData01, package=pkgname,    envir=parent.env(environment())) 
}

.onUnload <- function(libpath) {
  library.dynam.unload("Rbeast", libpath)
}

 