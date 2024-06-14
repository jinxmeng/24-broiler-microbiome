#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-11-17, 12:43:43
# modified date: 2023-12-03, 22:47:25

import argparse
import numpy as np
import re

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', type=str, metavar='profile', help='input profile')
    parser.add_argument('-t', type=str, metavar='taxon', help='input taxonomy')
    parser.add_argument('-k', type=str, metavar='KO', help='input gene_KO_res')
    parser.add_argument('-l', type=str, metavar='taxa level', default='family', help='input taxa_level, contain: \'phylum, class, order, family, genus, species\'')
    parser.add_argument('-o', type=str, metavar='output', help='output result')
    args = parser.parse_args()
    return args

def parse_KO(KO_f):
    KO_dict = {}
    with open(KO_f, 'r') as f:
        for line in f:
            line_list = line.strip().split('\t')
            KO_dict[line_list[0]] = line_list[1]
    return KO_dict

def parse_taxon(taxon_f, taxa_level):
    taxa_dict = {}
    if taxa_level == 'species':
        pattern = re.compile(r's__(.*?)(?:,|$)')
    else:
        pattern = re.compile(taxa_level + r'__(.*?);')
    with open(taxon_f, 'r') as f:
        for line in f:
            line_list = line.strip().split('\t')
            taxa_dict[line_list[0]] = re.findall(pattern, line_list[1])
    return taxa_dict

def main(in_f, KO_f, taxon_f, taxa_level, out_f):
    taxa_info = {'phylum':'p', 'class':'c', 'order':'o', 'family':'f', 'genus':'g', 'species':'s'}
    taxa_level = taxa_info[taxa_level]
    KO_dict = parse_KO(KO_f)
    taxa_dict = parse_taxon(taxon_f, taxa_level)
    sum_dict = {}
    header = ''
    with open(in_f, 'r') as f:
        header = next(f).strip()
        for line in f:
            line_list = line.strip().split('\t')
            gene = line_list[0]
            if gene in KO_dict:
                KO_num = KO_dict[gene]
            else:
                KO_num = "Unknown"
            if KO_num not in sum_dict:
                sum_dict.setdefault(KO_num, {})
            if gene in taxa_dict:
                taxon = taxa_dict[gene]
            else:
                taxon = ['Unknown']
            if len(taxon) == 1:
                if taxon[0] not in sum_dict[KO_num].keys():
                    sum_dict[KO_num][taxon[0]] = np.array(line_list[1:], dtype=float)
                else:
                    sum_dict[KO_num][taxon[0]] += np.array(line_list[1:], dtype=float)
            else:
                for i in taxon:
                    if i not in sum_dict[KO_num].keys():
                        sum_dict[KO_num][i] = 1/(len(taxon)) * np.array(line_list[1:], dtype=float)
                    else:
                        sum_dict[KO_num][i] += 1/(len(taxon)) * np.array(line_list[1:], dtype=float)

    with open(out_f, 'w') as f:
        header = header.split('\t')
        header = 'KO\ttaxa\t' + '\t'.join(header[1:])
        f.writelines(header + '\n')
        for i in sum_dict.keys():
            for k, v in sum_dict[i].items():
                f.writelines(i + '\t' + k + '\t' + '\t'.join([str(x) for x in v]) + '\n')

if __name__ == '__main__':
    args = get_args()
    main(args.i, args.k, args.t, args.l, args.o)

