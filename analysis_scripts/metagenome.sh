# fastp
cat ../sample_name | parallel -j 10 /share/data1/mjx/bin/flow_fastp_rmhost.sh /share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/00.cleandata/{}_1_clean.fq.gz,/share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/00.cleandata/{}_2_clean.fq.gz chicken {}

# megahit
cat ../sample_name | parallel -j 3 /share/data1/mjx/bin/flow_megahit.sh /share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/01.fastp_rmhost/{}_dehost_1.fq.gz,/share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/01.fastp_rmhost/{}_dehost_2.fq.gz 21,41,61,81,101,121,141 {}

# binning
cat ../sample_name | parallel -j 5 flow_metabat2.sh /share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/01.fastp_rmhost/{}_dehost_1.fq.gz,/share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/01.fastp_rmhost/{}_dehost_2.fq.gz /share/data1/mjx/proj/01.Broiler_Metagenome_2023_3/02.assembly/{}/{}.contigs.fa {} bin &

# checkm2
checkm2 predict --input ../03.metabat2/bin_fna/ --output-directory ckm2_res --tmpdir /tmp -x .fna --threads 36 --force
awk -F '\t' '$2>50 && $3<=10' ckm2_res/quality_report.tsv > ckm2_res_mq
awk -F '\t' '$2>=70 && $3<=5 && ($2-$3*5)>55' ckm2_res/quality_report.tsv > ckm2_res_qs
awk -F '\t' '$2>=90 && $3<=5' ckm2_res/quality_report.tsv > ckm2_res_hq
cut -f1 ckm2_res_qs | parallel -j 30 -q perl -e 'open I, "<../03.metabat2/bin_fna/$ARGV[0].fa";open O, ">bin_fna_983/$ARGV[0].fna";$i=1;while(<I>){chomp;if(/>/){print O ">$ARGV[0]_$i\n";$i++}else{print O "$_\n"}}' {}

# dRep
dRep compare genome_fna_cls_99 -g genome_fna_list -pa 0.9 -sa 0.99 -nc 0.30 -cm larger -p 40 --S_algorithm fastANI
dRep compare strains_fna_cls_95 -g strains_fna_list -pa 0.9 -sa 0.95 -nc 0.30 -cm larger -p 40 --S_algorithm fastANI

# gene predict
le ../05.dRep/strains_fna_list | grep -v HC | parallel -j 56 prodigal -f gff -p meta -a faa/{/.}.faa -d ffn/{/.}.ffn -o gff/{/.}.gff -i {} >/dev/null 2>/dev/null &
cat ../../02.broiler_bacteria_wgs_20230310/02.checkm2/isolate_fna_82_list | parallel -j 56 prodigal -f gff -p single -a faa/{/.}.faa -d ffn/{/.}.ffn -o gff/{/.}.gff -i {} >/dev/null 2>/dev/null &

# gtdb and tree
gtdbtk classify_wf --genome_dir ../05.dRep/strains_fna/ --out_dir gtdbtk_classify -x fna --cpus 100 --pplacer_cpus 24 --skip_ani_screen
pigz -kdc gtdbtk_classify/align/gtdbtk.bac120.user_msa.fasta.gz > gtdbtk.bac120.user_msa.fasta
gtdbtk infer --msa_file gtdbtk.bac120.user_msa.fasta --out_dir gtdbtk_infer --cpus 56 &

# profile
cat ../05.dRep/species_fna/* > species.fna
bowtie2-build species.fna species --threads 112
cat ../sample_name | parallel -j 8 ./flow_align.sh ../01.fastp_rmhost/{}_rmhost_1.fq.gz,../01.fastp_rmhost/{}_rmhost_2.fq.gz species align/{}
cat ../sample_name | parallel -j 8 ./flow_cvm.sh align/{}_sort.bam cvm/{}
grep 'overall' *log | perl -ne 'chomp;$_=~/(\S+).log:(\S+%)? /;print "$1\t$2\n"' > ../map_rate.tsv

# geneset
# metagenome gene predict
cat ../sample_name | parallel -j 10 -q perl -e 'open I, "seqkit seq -g -m 500 $ARGV[0] |";open O, ">prodigal/$ARGV[1]_m500.fna";$i=1;while(<I>){chomp;if(/>/){print O ">$ARGV[1]_$i\n"; $i+=1;}else{print O "$_\n"}} close O' ../02.assembly/{}.contigs.fa {}
cat ../sample_name | parallel -j 5 flow_prodigal.sh prodigal/{}_m500.fna prodigal/{}_m500 10
# 选择完整的基因，partial=00，长度选择大于100bp
cat ../../sample_name | parallel -j 10 -q perl -e 'open I, "seqkit seq -g -m 100 $ARGV[0] |";open O, ">$ARGV[1]";while(<I>){chomp;if(/partial=00/){$a=1}elsif(/>/ and !/partial=00/){$a=0};if($a==1){print O "$_\n"}}' {}_m500.ffn {}_genes.ffn
cat ../../sample_name | while read i;do perl -ne 'chomp;print "$1\n" if $_=~/>(\S+)/' ${i}_genes.ffn | seqkit grep -f - ${i}_m500.faa -o ${i}_genes.faa;done
cat *_genes.ffn | sed 's/.#.*//g' > gene.ffn
cat *_genes.faa | sed 's/.#.*//g' > gene.faa
# isolate gene partial=00，长度选择大于100bp
cat ../../06.prodigal/ffn/*ffn | seqkit seq -g -m 100 | perl -e 'while(<>){chomp;if(/partial=00/){$a=1}elsif(/>/ and !/partial=00/){$a=0};if($a==1){print "$_\n"}}' | sed 's/.#.*//g' > gene.ffn
cat ../../06.prodigal/faa/*faa | sed 's/.#.*//g' > gene.faa
# 所有基因
cat prodigal/gene.ffn isolates/gene.ffn > gene.ffn
cat prodigal/gene.faa isolates/gene.faa > gene.faa
# 基因聚类
mmseqs easy-cluster prodigal/gene.ffn mmseqs/gene /tmp --cluster-mode 2 --cov-mode 1 --min-seq-id 0.9 -c 0.9 --kmer-per-seq-scale 0.8 --threads 112
cp mmseqs/gene_rep_seq.fasta geneset.ffn
grep ">" geneset.ffn | sed 's/[>| ]//g' | seqkit grep -f - gene.faa > geneset.faa 
seqkit fx2tab -n -l -o geneset_ffn_len geneset.ffn
# abundance
bowtie2-build geneset.ffn geneset --threads 112 
# cat ../sample_name | parallel -j 7 ./flow_bowtie2_calu_rc.sh ../01.fastp_rmhost/{}_rmhost_1.fq.gz geneset rcs_all/{}
cat ../sample_name | parallel -j 8 ./flow_bowtie2_calu_rc_20m.sh ../01.fastp_rmhost/{}_rmhost_1.fq.gz geneset rcs/{}
# taxa
blastn -db /share/data1/database/ncbi/nt -query geneset.ffn -outfmt '6 std staxids' -num_threads 32 -num_alignments 20 -out nt/geneset_btn

# KEGG
diamond blastp -d /share/data1/database/KEGG/KEGG20230401_Prokaryotes.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 10 -p 12 -q ../10.geneset/geneset.faa -o geneset_kegg_btp >/dev/null 2>/dev/null
perl -ne '@s=split /\s+/;print $_ unless $a eq $s[0];$a=$s[0];' geneset_kegg_btp | perl -ne 'chomp;if($_=~/(.*)\s+.*\|(K\d+)\s+/){print "$1\t$2\n"}' > geneset_kegg_res
# KO abundance
transform_aggregate_kegg.pl geneset_kegg_res ../10.geneset/geneset_tpm KO_tpm

# CAZy
diamond blastp -d /share/data1/database/CAZy/CAZy_20220806.dmnd --outfmt 6 --min-score 60 --query-cover 50 --max-target-seqs 5 -p 12 -q ../10.geneset/geneset.faa -o geneset_CAZy_btp >/dev/null 2>/dev/null &
perl -ne '@s=split /\s+/;print $_ unless $a eq $s[0];$a=$s[0];' geneset_CAZy_btp | perl -ne 'chomp;if($_=~/(.*?)\s+.*?\|(.*?)\s/){print "$1\t$2\n"}' > geneset_CAZy_res
# abundance
transform_aggregate_CAZyme.pl geneset_CAZy_res ../10.geneset/geneset_tpm CAZy_tpm
transform_aggregate_CAZyme_filter.pl geneset_CAZy_res ../10.geneset/geneset_tpm CAZy_tpm_filter_m10 10
# transform_aggregate_CAZyme_filter.pl geneset_CAZy_res ../10.geneset/geneset_tpm CAZy_tpm_filter_m100 100
# annotation
perl -e '%h;open I, "geneset_CAZy_res"; while(<I>){chomp;@s=split/\s+/;$h{$s[0]}=1} while(<>){chomp;@s=split/\s+/;print "$_\n" if exists $h{$s[0]}}' ../11.KEGG/geneset_kegg_res > CAZy_genes_in_KEGG_annotation

# CARD
diamond blastp -d /share/data1/database/CARD/CARD_v3.2.7.dmnd --outfmt 6 --min-score 60 --query-cover 75 --id 80 --max-target-seqs 10 -p 56 -q ../10.geneset/geneset.faa -o geneset_CARD_btp >/dev/null 2>/dev/null
 perl -ne '@s=split /\s+/;print $_ unless $a eq $s[0];$a=$s[0];' geneset_CARD_btp | perl -ne 'chomp;if($_=~/(\S+)\s+ARO:(\d+)/){print "$1\t$2\n"}' > geneset_CARD_res
transform_aggregate_kegg.pl geneset_CARD_res ../10.geneset/geneset_tpm ARGs_tpm
transform_aggregate_kegg_filter.pl geneset_CARD_res ../10.geneset/geneset_tpm ARGs_tpm_m10 10
# transform_aggregate_kegg_filter.pl geneset_CARD_res ../10.geneset/geneset_tpm ARGs_tpm_m100 100
csvtk join -t --left-join --na NA -f1 geneset_CARD_res ../10.geneset/nt/geneset_latest_taxonomy.tsv > geneset_CARD_res_taxonomy

# MGEs
blastn -query ../10.geneset/geneset.ffn -db /share/data1/database/MGEdb/MGEs_FINAL_99perc_trim -evalue 1e-5 -max_target_seqs 5 -num_threads 32 -outfmt 6 -perc_identity 80 -qcov_hsp_perc 70 -out geneset_MGE_btn
perl -ne '@s=split /\s+/;print $_ unless $a eq $s[0];$a=$s[0];' geneset_MGE_btn | perl -ne 'chomp;if($_=~/(.*?)\s+(.*?)\s+/){print "$1\t$2\n"}' > geneset_MGE_res
transform_aggregate_kegg.pl geneset_MGE_res ../10.geneset/geneset_tpm MGEs_tpm
csvtk join -t --left-join --na NA -f1 geneset_MGE_res ../10.geneset/nt/geneset_latest_taxonomy.tsv > geneset_MGE_res_taxonomy

# find ARG-MGE combination
perl -e 'open I, "<$ARGV[0]";%h;while(<I>){chomp;$h{$_}=1} open I, "<$ARGV[1]";while(<I>){chomp;@s=split/\s+/;print "$_\n" if exists $h{$s[0]}}' args_gene.tsv gene_cluster.tsv > args_gene_all.tsv &
csvtk join -t -H --left-join -f 1 args_gene_all.tsv geneset_CARD_res > xx
find_gene_interval.py -bg ../10.geneset/geneset_gff.tsv -g1 args_gene_all_input.tsv -g2 mge_gene_all_input.tsv -o arg_interval &
