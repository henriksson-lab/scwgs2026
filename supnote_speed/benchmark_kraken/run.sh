
time kraken2 --db /data/henlab/kraken/standard-8 --thread 10 --output out --gzip-compressed --paired empty.R1.fq.gz empty.R2.fq.gz
#real    0m3,987s
#user    0m0,287s
#sys     0m3,960s

3.987 s  database loading time


time kraken2 --db /data/henlab/kraken/standard-8 --thread 10 --output out --gzip-compressed --paired fastp.1.R1.fq.gz fastp.1.R2.fq.gz
#39844353 sequences (8433.53 Mbp) processed in 154.733s (15450.3 Kseq/m, 3270.23 Mbp/m).
#  32309725 sequences classified (81.09%)
#  7534628 sequences unclassified (18.91%)

#real    2m38,801s
#user    8m30,562s
#sys     0m18,078s

38.801 s  ===> 34.814 s to process 40M reads

39 844 353  #40M reads


2742 cells

=> 10932.35  loading time if individual
(10932.35 + 34.814) / 38.801

282x faster


