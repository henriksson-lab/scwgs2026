"""Integration checks for Yanai read counting through the actual local Zorn wrapper.

Requires the local Bascet release build, R/Zorn, samtools and h5py.
Run: python3 -m unittest discover -s yanai_petriseq -p test_gene_counts.py -v
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import h5py


def matrix(path):
    with h5py.File(path) as handle:
        cells = handle["obs/_index"].asstr()[:].tolist()
        genes = handle["var/_index"].asstr()[:].tolist()
        indices = handle["X/indices"][:]
        offsets = handle["X/indptr"][:]
        data = handle["X/data"][:]
        values = {(cell,gene):0 for cell in cells for gene in genes}
        for i,cell in enumerate(cells):
            for j in range(offsets[i], offsets[i+1]):
                values[cell,genes[indices[j]]] = int(data[j])
        return cells,genes,values


class GeneCountIntegration(unittest.TestCase):
    def test_88_library_group_and_shard_order(self):
        """Exercise real natural shard discovery and helper metadata beyond shard 9."""
        zorn = Path(os.getenv("YANAI_ZORN_DIR", "~/github/zorn")).expanduser()
        repo = Path(__file__).resolve().parent.parent
        if not shutil.which("Rscript"):
            self.skipTest("Rscript is required")
        with tempfile.TemporaryDirectory(prefix="yanai-88-gene-count-") as directory:
            root = Path(directory)
            for i in range(1,89):
                (root/f"aligned_pos.{i}.bam").touch()
                (root/f"aligned_pos.{i}.bam.bai").touch()
            script = root/"ordering.R"
            script.write_text("\n".join([
                f"source({json.dumps(str(zorn/'R/zorn.R'))})",
                f"source({json.dumps(str(zorn/'R/bascet_file.R'))})",
                f"source({json.dumps(str(repo/'yanai_petriseq/count_genes.R'))})",
                f"setwd({json.dumps(str(root))})",
                "prefixes <- paste0('library',seq_len(88))",
                "debstat <- list(group_prefix=prefixes,list_picked_cells=lapply(prefixes,function(x)paste0(x,'_a1_b1_c1')))",
                "BascetCountFeature <- function(bascetRoot,inputName,outputName,gffFile,useFeature,attrGeneId,attrGeneName,countMode,cellFiles,runner,bascetInstance) {",
                " stopifnot(countMode=='reads',attrGeneId=='locus_tag',attrGeneName=='Name',useFeature=='gene')",
                " for(i in seq_len(88)) { stopifnot(identical(readLines(cellFiles[i]),debstat$list_picked_cells[[i]])); file.create(file.path(bascetRoot,paste0(outputName,'.',i,'.h5'))) }; NULL",
                "}",
                "WaitForJob <- function(job) invisible(NULL)",
                "result <- RunYanaiGeneCounts(getwd(),debstat,'annotation.gff',NULL,NULL)",
                "stopifnot(identical(result$library_prefix,prefixes),identical(result$alignment,paste0('aligned_pos.',seq_len(88),'.bam')),identical(result$counts,paste0('genecount.',seq_len(88),'.h5')))",
                "saved <- read.csv('yanai_gene_count_libraries.csv')",
                "stopifnot(nrow(saved)==88,identical(saved$library_prefix,prefixes))",
            ])+"\n")
            subprocess.run(["Rscript",str(script)],check=True)

    def test_reads_overlap_boundaries_flags_and_zero_cells(self):
        bascet_dir = Path(os.getenv("YANAI_BASCET_DIR", "~/github/bascet")).expanduser()
        binary = bascet_dir / "target/release/bascet"
        zorn = Path(os.getenv("YANAI_ZORN_DIR", "~/github/zorn")).expanduser()
        if not binary.is_file() or not shutil.which("samtools") or not shutil.which("Rscript"):
            self.skipTest("Local Bascet, samtools and Rscript are required")
        repo = Path(__file__).resolve().parent.parent
        with tempfile.TemporaryDirectory(prefix="yanai-gene-count-") as directory:
            root = Path(directory)
            gff = root / "genes.gff"
            gff.write_text(
                "##gff-version 3\n"
                "chr\tNCBI\tgene\t10\t30\t.\t+\t.\tID=gene-A;locus_tag=A;Name=GeneA\n"
                "chr\tNCBI\tgene\t25\t45\t.\t-\t.\tID=gene-B;locus_tag=B;Name=GeneB\n"
                "chr\tNCBI\tgene\t60\t60\t.\t+\t.\tID=gene-C;locus_tag=C;Name=GeneC\n"
            )
            first = "D1_lib1_a1_b1_c1"
            second = "D1_lib1_a2_b1_c1"
            zero = "D1_lib1_a3_b1_c1"
            records = []
            def add(cell, flag, pos, cigar="1M", seq="A", mapq=60):
                records.append("\t".join(map(str, (
                    cell+":AAAAAAA:"+str(len(records)), flag, "chr", pos,
                    mapq, cigar, "*", 0, 0, seq, "I"*len(seq)
                ))))
            add(first,0,10)           # first gene base must count
            add(first,0,15)
            add(first,0,15)           # duplicate UMI still counts as a read
            add(first,1024,15)        # duplicate flag retained in raw read mode
            add(first,0,20,mapq=0)    # MAPQ=0 primary alignment retained
            add(first,0,25,"10M","A"*10) # full count to both overlapping genes
            add(second,16,40)         # opposite strand also counted
            add(second,0,45)          # last gene base counts
            add(second,0,60)          # single-base gene
            add(first,256,15)         # secondary excluded
            add(first,2048,25)        # supplementary excluded
            records.append(f"{zero}:AAAAAAA:unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI")
            sam = root / "fixture.sam"
            sam.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:100\n"+"\n".join(records)+"\n")
            bam = root / "aligned_pos.1.bam"
            subprocess.run(["samtools","sort","-o",str(bam),str(sam)],check=True)
            subprocess.run(["samtools","index",str(bam)],check=True)
            script = root / "count.R"
            script.write_text("\n".join([
                f"zorn <- {json.dumps(str(zorn))}",
                "for(f in c('job_general.R','job_local.R','job_slurm.R','jobscript.R','bascet_file.R','zorn.R','shell.R','container.R','refgenome.R')) source(file.path(zorn,'R',f))",
                f"source({json.dumps(str(repo/'yanai_petriseq/count_genes.R'))})",
                f"setwd({json.dumps(str(root))})",
                f"instance <- getBascetDevDir({json.dumps(str(bascet_dir))})",
                "runner <- LocalRunner(direct=TRUE, showScript=FALSE, ncpu='2', mem='2g')",
                "debstat <- list(group_prefix='D1_lib1',list_picked_cells=list(c("+",".join(json.dumps(c) for c in (first,second,zero))+")))",
                "RunYanaiGeneCounts(getwd(),debstat,'genes.gff',runner,instance)",
            ])+"\n")
            subprocess.run(["Rscript",str(script)],check=True)
            cells,genes,values = matrix(root/"genecount.1.h5")
            self.assertEqual(cells,[first,second,zero])
            self.assertEqual(genes,["A","B","C"])
            self.assertEqual([values[first,g] for g in genes],[6,1,0])
            self.assertEqual([values[second,g] for g in genes],[0,2,1])
            self.assertEqual([values[zero,g] for g in genes],[0,0,0])
            # Historical default is still UMI counting; identical UMIs collapse.
            subprocess.run([str(binary),"countfeature","-i",str(bam),"-g",str(gff),
                            "--attr-id","locus_tag","--attr-name","Name",
                            "--cells",str(root/"gene_cells.1.txt"),"-@","2",
                            "-o",str(root/"umi.h5")],check=True)
            _,_,umi = matrix(root/"umi.h5")
            self.assertEqual(umi[first,"A"],1)
            self.assertTrue((root/"yanai_gene_count_libraries.csv").is_file())


if __name__ == "__main__":
    unittest.main()
