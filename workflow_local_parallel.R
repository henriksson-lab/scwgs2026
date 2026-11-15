# Bounded local array execution using the same Zorn JobScripts as LocalRunner.
# Loaded after job_local.R; opt-in per workflow. Resources are PER array task.
setClass("DatasetArrayRunner",contains="LocalRunner",slots=c(workers="integer"))
setMethod("RunJob","DatasetArrayRunner",function(runner,jobname,bascetInstance,cmd,arraysize) {
  stopifnot(is.jobscript(cmd),runner@direct,runner@workers>=1L)
  tmproot <- GetBascetTempDir(bascetInstance)
  task <- function(i) {
    taskId <- i-1L
    tempdir <- file.path(tmproot,as.character(taskId))
    dir.create(tempdir,recursive=TRUE,showWarnings=FALSE)
    env <- list(BASCET_TEMPDIR=tempdir)
    if(bascetInstance@logToFile) {
      dir.create("logs",showWarnings=FALSE)
      env$BASCET_LOGFILE <- paste0("logs/",jobname,".",taskId,".log")
    }
    tryCatch({
      runJobScriptLocal(cmd,task_id=taskId,env=env)
      list(ok=TRUE)
    },error=function(e) list(ok=FALSE,message=paste0("Local array task ",taskId,
                                                     " (",jobname,"): ",conditionMessage(e))))
  }
  message("Running ",jobname," with at most ",runner@workers," concurrent array tasks")
  result <- parallel::mclapply(seq_len(arraysize),task,
                               mc.cores=min(runner@workers,arraysize),
                               mc.preschedule=FALSE,mc.set.seed=FALSE)
  ok <- vapply(result,function(x)is.list(x)&&isTRUE(x$ok),logical(1))
  if(any(!ok)) {
    failures <- vapply(result[!ok],function(x) {
      if(is.list(x)&&!is.null(x$message)) x$message else paste(x,collapse=" ")
    },character(1))
    stop(paste(failures,collapse="\n"),call.=FALSE)
  }
  new_no_job()
})
DatasetArrayRunner <- function(runner,workers) {
  stopifnot(is(runner,"LocalRunner"),length(workers)==1L,!is.na(workers),workers>=1L)
  if(workers==1L) return(runner)
  if(.Platform$OS.type!="unix") stop("Parallel local arrays require Unix fork support")
  new("DatasetArrayRunner",ncpu=runner@ncpu,mem=runner@mem,direct=runner@direct,
      showScript=runner@showScript,workers=as.integer(workers))
}
