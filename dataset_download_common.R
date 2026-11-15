# Shared R/shell download utilities. Callers hold their dataset's flock lock.
ParseDownloadOptions <- function(defaults,flags=character(),repeated=character()) {
  args <- commandArgs(trailingOnly=TRUE)
  i <- 1L
  while(i<=length(args)) {
    arg <- args[i]
    if(!startsWith(arg,"--")) stop("Expected an option: ",arg)
    parts <- strsplit(sub("^--","",arg),"=",fixed=TRUE)[[1]]
    key <- parts[1]
    if(!key %in% names(defaults)) stop("Unknown option: ",arg)
    if(key %in% flags) {
      if(length(parts)!=1L) stop("Flag takes no value: ",arg)
      defaults[[key]] <- TRUE
    } else {
      if(length(parts)>1L) value <- paste(parts[-1],collapse="=") else {
        i <- i+1L
        if(i>length(args)) stop("Missing option value: ",arg)
        value <- args[i]
      }
      if(key %in% repeated) defaults[[key]] <- c(defaults[[key]],value) else
        defaults[[key]] <- value
    }
    i <- i+1L
  }
  defaults
}
DownloadCommand <- function(command,args) {
  status <- system2(command,shQuote(as.character(args)))
  if(status!=0L) stop(command," failed with exit status ",status)
}
VerifyDownload <- function(path,bytes,md5) {
  if(!file.exists(path) || file.info(path)$size!=as.numeric(bytes) ||
     unname(tools::md5sum(path))!=md5) stop("Source size/MD5 mismatch: ",path)
  invisible(path)
}
FetchDownload <- function(path,url,bytes,md5,download=TRUE) {
  if(!file.exists(path)) {
    if(!download) stop("Missing source: ",path,"; pass --download to retrieve it")
    partial <- paste0(path,".part")
    DownloadCommand("curl",c("--fail","--location","--silent","--show-error",
      "--retry","8","--retry-delay","5","--connect-timeout","30",
      "--continue-at","-","--output",partial,url))
    VerifyDownload(partial,bytes,md5)
    if(!file.link(partial,path)) stop("Cannot publish ",path)
    unlink(partial)
  }
  VerifyDownload(path,bytes,md5)
}
PublishDownloadJSON <- function(record,path) {
  temp <- tempfile("receipt-",tmpdir=dirname(path))
  on.exit(unlink(temp))
  jsonlite::write_json(record,temp,auto_unbox=TRUE,pretty=TRUE)
  if(!file.link(temp,path)) stop("Cannot publish receipt: ",path)
}
CopyDownloadMetadata <- function(source,target) {
  if(normalizePath(source,mustWork=TRUE)==normalizePath(target,mustWork=FALSE)) return(invisible())
  if(file.exists(target)) {
    if(unname(tools::md5sum(source))!=unname(tools::md5sum(target))) {
      stop("Refusing to replace different metadata: ",target)
    }
  } else if(!file.copy(source,target)) stop("Cannot copy metadata: ",target)
}
RunDownloadTasks <- function(tasks,workers,fun) {
  workers <- as.integer(workers)
  stopifnot(length(workers)==1L,!is.na(workers),workers>0L,length(tasks)>0L)
  result <- parallel::mclapply(tasks,function(task) tryCatch({
    fun(task)
    list(ok=TRUE)
  },error=function(e) list(ok=FALSE,message=conditionMessage(e))),
    mc.cores=min(workers,length(tasks)),mc.preschedule=FALSE,mc.set.seed=FALSE)
  ok <- vapply(result,function(x)is.list(x)&&isTRUE(x$ok),logical(1))
  if(any(!ok)) stop(paste(vapply(result[!ok],function(x) {
    if(is.list(x)) x$message else as.character(x)
  },character(1)),collapse="\n"))
}
