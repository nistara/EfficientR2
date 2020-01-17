# run as
#  e = new.env()
#  source("rm.R", e)
# to avoid the rm() clearing the global environment.
# Also, can save the resulting environment to RDS file.
rm(list = ls(all = TRUE))
