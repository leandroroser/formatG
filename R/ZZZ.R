.onUnload <- function (libpath) {
  library.dynam.unload("formatter", libpath)
}