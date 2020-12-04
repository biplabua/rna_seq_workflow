## ----eval = FALSE----------------------------------------------------------
#  library(ensembldb)
#  ## Load the EnsDb package that should be installed on the MySQL server
#  library(EnsDb.Hsapiens.v86)
#  
#  ## Call the useMySQL method providing the required credentials to create
#  ## databases and inserting data on the MySQL server
#  edb_mysql <- useMySQL(EnsDb.Hsapiens.v86, host = "localhost",
#                        user = "userwrite", pass = "userpass")
#  
#  ## Use this EnsDb object
#  genes(edb_mysql)
#  

## ----eval = FALSE----------------------------------------------------------
#  library(ensembldb)
#  library(RMariaDB)
#  
#  ## Connect to the MySQL database to list the databases.
#  dbcon <- dbConnect(MariaDB(), host = "localhost", user = "readonly",
#                     pass = "readonly")
#  
#  ## List the available databases
#  listEnsDbs(dbcon)
#  
#  ## Connect to one of the databases and use that one.
#  dbcon <- dbConnect(MariaDB(), host = "localhost", user = "readonly",
#                     pass = "readonly", dbname = "ensdb_hsapiens_v75")
#  edb <- EnsDb(dbcon)
#  edb
#  

## ----sessionInfo-----------------------------------------------------------
sessionInfo() 

