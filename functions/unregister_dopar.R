### Unregister any old open parallel backends
### Fall 2021
### adam.milton.morgan@gmail.com

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
