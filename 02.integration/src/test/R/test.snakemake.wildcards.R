logger::log_threshold(logger::TRACE)
log_file <- snakemake@log[[1]]
logger::log_appender(logger::appender_file(log_file))


afnm <- snakemake@input[["afnm"]]
str(afnm)
logger::log_info(afnm)

bfnm <- snakemake@input[["bfnm"]]
str(bfnm)
print(bfnm)

params <- snakemake@params[["a"]]
str(params)
print(params)

print(snakemake@wildcards)



