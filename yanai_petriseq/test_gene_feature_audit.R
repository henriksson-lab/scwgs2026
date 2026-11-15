#!/usr/bin/env Rscript
# Exercise the audit on the complete pilot against independently read native H5.
# Only the isolated /tmp fixture is mutated by the negative controls.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table)})
arg <- grep("^--file=",commandArgs(),value=TRUE)
helperDir <- dirname(normalizePath(sub("^--file=","",arg[1])))
repoDir <- dirname(helperDir)
source(file.path(helperDir,"verify_variable_gene_coverage.R"))
source(file.path(helperDir,"gene_profile_fingerprint.R"))
source(file.path(helperDir,"verify_cache_native_counts.R"))
main <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
pilot <- file.path(main,"validation_full_D1_lib1")
fixture <- tempfile("yanai-feature-audit-",tmpdir="/tmp")
dir.create(file.path(fixture,"umap"),recursive=TRUE)
object <- readRDS(file.path(pilot,"cache_adata_gene.RDS"))
stopifnot(ncol(object)==10018L,nrow(object)==4550L)
features <- VariableFeatures(object)
native <- file.path(pilot,"genecount.1.h5")
rawGenes <- as.character(rhdf5::h5read(native,"/var/_index"))
genes <- gsub("_","-",rawGenes,fixed=TRUE)
cells <- as.character(rhdf5::h5read(native,"/obs/_index"))
counts <- methods::new("dgCMatrix",Dim=as.integer(c(length(genes),length(cells))),
  Dimnames=list(genes,cells),p=as.integer(rhdf5::h5read(native,"/X/indptr")),
  i=as.integer(rhdf5::h5read(native,"/X/indices")),
  x=as.numeric(rhdf5::h5read(native,"/X/data")))
methods::validObject(counts)
stopifnot(all(features %in% genes),all(colnames(object) %in% cells))
# Reconstruct the entire pilot VST selection independently, including Zorn's
# feature sorting before Seurat normalizes gene names.
sortedGenes <- as.character(sort(rawGenes))
independentCounts <- counts[match(sortedGenes,rawGenes),colnames(object),drop=FALSE]
rownames(independentCounts) <- sortedGenes
set.seed(42L)
independent <- CreateSeuratObject(counts=independentCounts,project="Yanai",min.cells=0,min.features=0)
independent <- independent[Matrix::rowSums(GetAssayData(independent,layer="counts"))>0,]
independent <- NormalizeData(independent)
independent <- FindVariableFeatures(independent,selection.method="vst",nfeatures=2000L)
stopifnot(identical(unname(VariableFeatures(independent)),unname(features)))
rm(independent,independentCounts);gc()
counts <- counts[features,colnames(object),drop=FALSE]
coverage <- data.table(cell_id=colnames(counts),variable_genes=as.numeric(Matrix::colSums(counts>0)),
                       variable_gene_reads=as.numeric(Matrix::colSums(counts)))
rm(counts);gc()
stopifnot(file.symlink(native,file.path(fixture,"genecount.1.h5")),
          file.copy(file.path(pilot,"yanai_gene_count_libraries.csv"),fixture),
          file.copy(file.path(pilot,"umap","gene_cell_selection.csv"),file.path(fixture,"umap")))
out <- file.path(fixture,"umap")
writeLines(features,file.path(out,"gene_variable_features.reconstructed.txt"))
genome <- ifelse(grepl("^b[0-9]+$",features),"Escherichia coli",
  ifelse(grepl("^SAUSA300-RS",features),"Staphylococcus aureus","Bacillus subtilis"))
fwrite(data.table(feature=features,reference_species=genome),
       file.path(out,"gene_variable_features.tsv"),sep="\t")
coveragePath <- file.path(out,"gene_variable_gene_coverage.tsv")
fwrite(coverage,coveragePath,sep="\t")
profilePath <- file.path(out,"gene_variable_gene_coverage.json")
profile <- list(status="fixture_pending",expressed_genes=nrow(object),
  input_fingerprint=YanaiGeneProfileFingerprint(fixture),feature_genomes=as.list(table(genome)))
jsonlite::write_json(profile,profilePath,auto_unbox=TRUE,pretty=TRUE)
# This fixture tests feature/count auditing; actual production warning parsing
# is separately checked against the production log.
writeLines("Audit fixture",file.path(fixture,"all_workflow.log"))
Sys.setenv(YANAI_ROOT=fixture)
stopifnot(system2(file.path(R.home("bin"),"Rscript"),
  shQuote(file.path(helperDir,"annotate_variable_genes.R")))==0L)
positive <- VerifyYanaiVariableGeneCoverage(fixture,object,repoDir)
nativeCheck <- VerifyYanaiCacheNativeCounts(fixture,object)
stopifnot(nativeCheck$exact_every_cell_gene_count)
changedObject <- object
changedCounts <- GetAssayData(changedObject,layer="counts")
changedCounts[1,1] <- changedCounts[1,1]+1
SeuratObject::LayerData(changedObject,layer="counts") <- changedCounts
cacheCountRejected <- tryCatch({VerifyYanaiCacheNativeCounts(fixture,changedObject);FALSE},
  error=function(e) grepl("Cached gene counts differ",conditionMessage(e),fixed=TRUE))
stopifnot(isTRUE(cacheCountRejected))
rm(changedObject,changedCounts);gc()
stopifnot(positive$gene_feature_coverage$exact_actual_feature_list,
          positive$gene_feature_coverage$exact_actual_variable_gene_counts)
original <- coverage$variable_gene_reads[1]
coverage[1,variable_gene_reads:=original+1]
fwrite(coverage,coveragePath,sep="\t")
countRejected <- tryCatch({VerifyYanaiVariableGeneCoverage(fixture,object,repoDir);FALSE},
  error=function(e) grepl("variable_gene_reads",conditionMessage(e),fixed=TRUE))
stopifnot(isTRUE(countRejected))
coverage[1,variable_gene_reads:=original]
fwrite(coverage,coveragePath,sep="\t")
wrong <- features
wrong[1:2] <- wrong[2:1]
writeLines(wrong,file.path(out,"gene_variable_features.reconstructed.txt"))
featureRejected <- tryCatch({VerifyYanaiVariableGeneCoverage(fixture,object,repoDir);FALSE},
  error=function(e) grepl("VariableFeatures",conditionMessage(e),fixed=TRUE))
stopifnot(isTRUE(featureRejected))
writeLines(features,file.path(out,"gene_variable_features.reconstructed.txt"))
profile$input_fingerprint <- "stale_fixture"
jsonlite::write_json(profile,profilePath,auto_unbox=TRUE,pretty=TRUE)
stale <- VerifyYanaiVariableGeneCoverage(fixture,object,repoDir)
stopifnot(stale$gene_feature_coverage$status=="stale_optional_diagnostic")
profile$input_fingerprint <- YanaiGeneProfileFingerprint(fixture)
jsonlite::write_json(profile,profilePath,auto_unbox=TRUE,pretty=TRUE)
VerifyYanaiVariableGeneCoverage(fixture,object,repoDir)
report <- list(status="passed",scope="every selected variable gene in the complete native-H5 pilot",
  cells=ncol(object),features=length(features),native_h5_positive_control=TRUE,
  independently_sorted_native_vst_matches=TRUE,
  all_cached_genes_match_native=TRUE,changed_cached_count_rejected=cacheCountRejected,
  changed_count_rejected=countRejected,changed_feature_order_rejected=featureRejected,
  stale_optional_diagnostic_handled=TRUE,production_files_modified=FALSE,fixture=fixture)
jsonlite::write_json(report,file.path(fixture,"test_result.json"),auto_unbox=TRUE,pretty=TRUE)
print(report)
