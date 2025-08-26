.onAttach <- function(libname, pkgname) {
    current_version <- packageVersion(pkgname)
    required_version <- "1.11"
    
    if (current_version <= package_version(required_version)) {
        packageStartupMessage(
            paste0("WARNING: ", pkgname, " version ", current_version, 
                   " is outdated. Please update to a version > ", required_version, ".")
        )
    }
}