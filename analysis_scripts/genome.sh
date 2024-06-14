# spades
le ../sample_name | parallel -j 4 spades.py --careful --cov-cutoff auto -t 20 -1 ../00.cleandata/{}_1_clean.fq.gz -2 ../00.cleandata/{}_2_clean.fq.gz -o {}

# checkm2
checkm2 predict --input isolate_98 --output-directory ckm2_res --tmpdir /tmp -x .fna --threads 36 --force
awk -F '\t' '$2>50 && $3<=10' ckm2_res/quality_report.tsv > ckm2_res_mq
awk -F '\t' '$2>=70 && $3<=5 && ($2-$3*5)>55' ckm2_res/quality_report.tsv > ckm2_res_qs
awk -F '\t' '$2>=90 && $3<=5' ckm2_res/quality_report.tsv > ckm2_res_hq
cut -f1 ckm2_res_qs | parallel -j 10 -q perl -e 'open I, "seqkit seq -g -m 500 ../01.assembly/$ARGV[0].scaffolds.fasta |";open O, ">isolate_fna_82/$ARGV[0].fna";while(<I>){chomp;if(/>NODE_(\d+)_length/){print O ">$ARGV[0]_$1\n"}else{print O "$_\n"}}' {}
