#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-12-21, 22:04:06
# modified date: 2023-12-25, 19:53:56

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-bg', type=str, metavar='background', help='gff-like table with field contig|gene|start|end|direction.')
    parser.add_argument('-g1', type=str, metavar='marker gene table 1', help='a tab-delimited table with column contig|gene|label')
    parser.add_argument('-g2', type=str, metavar='marker gene table 2', help='a tab-delimited table with column contig|gene|label')
    parser.add_argument('-l', type=int, metavar='interval_len', default=5000, help='interval lenth for up- and down-stream of marker gene 1')
    parser.add_argument('-o', type=str, metavar='out_prefix', help='output graph file')
    args = parser.parse_args()
    return args

def get_dict(in_f, k_col = 0, v_col = 1):
    dct = {}
    with open (in_f, 'r') as f:
        for l in f:
            s = l.split()
            k, v = s[k_col], s[v_col]
            if k not in dct:
                dct[k] = v
            else:
                continue
    return dct
def get_list(in_f, col = 0):
    lst = []
    with open (in_f, 'r') as f:
        for l in f:
            x = l.split()[col]
            if x not in lst:
                lst.append(x)
            else:
                continue
    return lst

def get_bg(in_f, ref_list):
    gene_dict = {}
    interval_dict = {}
    with open (in_f, 'r') as f:
        for l in f:
            contig, gene, start, end, direct = l.split()
            if contig in ref_list:
                if contig not in gene_dict:
                    gene_dict[contig] = []
                    interval_dict[contig] = []
                gene_dict[contig].append(gene)
                interval_dict[contig].append([int(start), int(end), str(direct)])
            else:
                continue
    return gene_dict, interval_dict

def out_g(out_dict, out_g, ref_list):
    out_g = open(out_g, 'w')
    for i in out_dict.keys():
        contig = '_'.join(i.split('_')[0:-1])
        seqs = sorted([x.split('_')[-1] for x in out_dict[i]], key = int)
        seqs = [contig + '_' + str(x) for x in seqs] 
        out_g.write(i + '\t')
        for j in seqs:
            x = j
            if x == i and x in ref_list:
                x += ' ARGs|MGEs'
            elif x == i and x not in ref_list:
                x += ' ARGs'
            elif x != i and x in ref_list:
                x += ' MGEs'
            out_g.write('---(' + x + ')---')
        out_g.write('\n')
    out_g.close()
    

def out_f(out_dict, out_f, ref_list):
    out_f = open(out_f, 'w')
    for i in out_dict.keys():
        contig = '_'.join(i.split('_')[0:-1])
        seqs = sorted([x.split('_')[-1] for x in out_dict[i]], key = int)
        seqs = [contig + '_' + str(x) for x in seqs] 
        for j in seqs:
            x = j
            if x == i and x in ref_list:
                x += '\tARGs|MGEs'
            elif x == i and x not in ref_list:
                x += '\tARGs'
            elif x != i and x in ref_list:
                x += '\tMGEs'
            else:
                x += '\tOther'
            out_f.write('\t'.join([i, x]) + '\t' + '\t'.join(str(s) for s in out_dict[i][j]) + '\n')
    out_f.close()

def main(bg, g1, g2, l, o):
    g1_contig = get_list(g1, col = 0)
    g2_contig = get_list(g2, col = 0)
    tmp_contig = list(set(g1_contig).intersection(set(g2_contig))) # 只保留有交集的contig
    
    g1_dict = get_dict(g1, k_col = 1, v_col = 0) # 只有哪些基因在交集的contig才会被分析
    tmp_k = []
    for k, v in g1_dict.items():
        if v not in tmp_contig:
            tmp_k.append(k)
    for i in tmp_k:
        del g1_dict[i]
    
    g2_gene = get_list(g2, col = 1) # MGE基因列表
    # 背景基因信息 bg_gene_dict[contig] = [gene_name] bg_interval_dict[contig] = [gene_loc]
    bg_gene_dict, bg_interval_dict = get_bg(bg, tmp_contig)
    
    # 滑动基因进行ARGs上下游某个区间判断基因
    out_dict = {}
    for k, v in g1_dict.items(): 
        gene_index = bg_gene_dict[v].index(k)
        gene_interval = bg_interval_dict[v][gene_index]
        n_gene = len(bg_gene_dict[v])
        
        '''
        1. 如果目的基因单独存在某个contig上，不进行区间判断
        2. 如果目的基因位于某个contig的最左端，只判断另一侧 
        3. 如果目的基因位于某个contig的最右端，只判断另一侧
        4. 如果目的基因位于contig的中间，两侧均判断
        '''
        out_dict[k] = {}
        out_dict[k][k] = gene_interval
        i = 1
        if gene_index == 0 and n_gene == 1:
            continue
        elif gene_index == 0 and n_gene > 1:
            lim_right = gene_interval[1] + l
            while bg_interval_dict[v][gene_index + i][0] <= lim_right:
                out_dict[k][bg_gene_dict[v][gene_index + i]] = bg_interval_dict[v][gene_index + i]
                i+=1
                if gene_index + i > (n_gene - 1):
                    break
        elif gene_index == (n_gene - 1) and n_gene > 1:
            lim_left = gene_interval[0] - l
            while bg_interval_dict[v][gene_index - i][1] >= lim_left:
                out_dict[k][bg_gene_dict[v][gene_index - i]] = bg_interval_dict[v][gene_index - i]
                i+=1
                if gene_index - i < 0:
                    break
        else:
            lim_left = gene_interval[0] - l
            lim_right = gene_interval[1] + l
            while bg_interval_dict[v][gene_index - i][1] >= lim_left:
                out_dict[k][bg_gene_dict[v][gene_index - i]] = bg_interval_dict[v][gene_index - i]
                i+=1
                if gene_index - i < 0:
                    break
            i = 1
            while bg_interval_dict[v][gene_index + i][0] <= lim_right:
                out_dict[k][bg_gene_dict[v][gene_index + i]] = bg_interval_dict[v][gene_index + i]
                i+=1
                if gene_index + i > (n_gene - 1):
                    break

    # 判断MGE基因是否在这些范围内，在的话，就保留下来这个out_dict的子字典
    # 做集合运算
    # 循环修改字典时，迭代对象不能时字典
    y = set(g2_gene)
    for k in list(out_dict.keys()): 
        x = set(out_dict[k].keys())
        if len(x & y) == 0:
            del out_dict[k]
    
    # 输出结果
    out_f(out_dict, o + '_find.tsv', g2_gene)
    out_g(out_dict, o + '_graph.tsv', g2_gene)

if __name__ == '__main__':
    args = get_args()
    main(args.bg, args.g1, args.g2, args.l, args.o)

