SarimStartupMessage <- function()
{
    # Startup message obtained as
    # > figlet -f slant MCLUST
    msg <- c(paste0(
        "                _
   ________ ___(_) ___
  / __/ _  | __| |/   \\
  \\__ \\(_) | | | | Y Y |
  /___/\\_|_|_| |_|_|_|_|  version ",
        packageVersion("Sarim")),
        "\n  Structured Additive Regression Models using Iterative Methods.")
    return(msg)
}

.onAttach <- function(lib, pkg)
{
    # unlock .mclust variable allowing its modification
    ## unlockBinding("sarim", asNamespace("sarim"))
    # startup message
    msg <- SarimStartupMessage()
    if (!interactive())
        msg[1] <- paste("Package 'Sarim' version", utils::packageVersion("Sarim"))
    packageStartupMessage(msg)
    invisible()
}