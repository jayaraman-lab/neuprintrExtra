
cacheEnv <- new.env()
assign("storedTypes",neuprint_get_meta(1) %>% mutate(databaseType=character()),envir=cacheEnv)
assign("storedMeta",neuprint_get_meta(1) %>% mutate(databaseType=character()),envir=cacheEnv)
assign("storedRoiInfo",data.frame(bodyid=double(),roi=character(),pre=integer(),post=integer(),downstream=integer(),stringsAsFactors = FALSE),envir=cacheEnv)