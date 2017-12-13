nohup awk '{if ($4==36 && $3-$2 > 100) print $0;}' mm9_36_chr_fareads.bdg > mm9_36_chr_fareads_ge100mapregion.bed &
nohup awk '{if ($4==36 && $3-$2 > 100) print $0;}' hg19_36_chr_fareads.bdg > hg19_36_chr_fareads_ge100mapregion.bed &

nohup awk '{if ($4==50 && $3-$2 > 100) print $0;}' mm9_50_chr_fareads.bdg > mm9_50_chr_fareads_ge100mapregion.bed &
nohup awk '{if ($4==50 && $3-$2 > 100) print $0;}' hg19_50_chr_fareads.bdg > hg19_50_chr_fareads_ge100mapregion.bed &

