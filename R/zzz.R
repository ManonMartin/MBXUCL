# WelcomeMessage

.onLoad <- function(...) {
  suppressPackageStartupMessages(require("modeest"))
  suppressPackageStartupMessages(require("clusterSim"))
  packageStartupMessage("Package MBXUCL was correclty loaded\n")
  # message("message from .onLoad via message")
  # packageStartupMessage("message from .onLoad via
                        # packageStartupMessage\n")
}
