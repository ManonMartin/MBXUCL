# WelcomeMessage

.onLoad <- function(...) {
  suppressPackageStartupMessages(require(modeest))
  suppressPackageStartupMessages(library(modeest))
  cat("Package MBXUCL was correclty loaded\n")
  # message("message from .onLoad via message")
  # packageStartupMessage("message from .onLoad via
                        # packageStartupMessage\n")
}
