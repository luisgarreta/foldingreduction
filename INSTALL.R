#
# Create foldred profile according to current directory
# 
#

# Create foldred profile
message ("Creating foldred profile...")

# Write into .profile
profileFile = paste0 (path.expand ("~"), "/.bashrc")
sink (profileFile, append=T)
writeLines ("")
writeLines ("#-------- foldred tool for reduction of protein folding trajectories ---------")
writeLines (paste0 ("export FOLDING_REDUCTION_HOME=", getwd()))
writeLines ("export PATH=$PATH:$FOLDING_REDUCTION_HOME")
sink ()

message ("")
message ("Folding Reduction is ready to use, right after installed!")
message ("")

# Install R libraries
#libpath = paste0 (FOLDING_REDUCTION_HOME, "/opt/Rlibs")
message ("\n\nInstalling R libraries...\n\n")

#.libPaths (libpath)

#if (!require("pacman")) install.packages('pacman', lib=libpath, repos='http://cran.us.r-project.org')
if (!require("pacman")) install.packages('pacman', repos='http://cran.us.r-project.org')

pacman::p_load("parallel","cluster")

system ("R CMD SHLIB tmscorelg.f")
