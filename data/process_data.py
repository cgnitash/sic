# data from https://hivdb.stanford.edu/pages/genopheno.dataset.html

import pandas as pd
import math
import csv
from collections import defaultdict

data = pd.read_csv('PI_DataSet.txt', sep='\t', header=0)
#data = pd.read_csv('dummy.tsv', sep='\t', header=0)

sequence = 'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'

#print(data.head())
#print(type(sequence))
#print(data.iloc[[0]])
#line = data.iloc[[0]]
#print(line['FPV'][0])

treatments = ['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV']
treat_data = defaultdict(list)

for _, line in data.iterrows():
    chars = list(sequence)
    for i in range(1, 100):
        pos = line['P' + str(i)][0]
        if pos != '-':
            chars[i - 1] = pos
    cur_seq = ''.join(chars)

    for treat in treatments:
        treat_val = line[treat]
        if not math.isnan(treat_val):
            treat_data[treat].append((cur_seq, treat_val, line['SeqID']))

for treat, isolate_list in treat_data.items():
    with open('PI_' + treat + '.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        for isolate in isolate_list:
            writer.writerow(isolate)
print(type(treat_data))


'''
print(diff_posn)

all_seqs = []
all_seqs.append(sequence)
for pos, chars in diff_posn:
    temp_seqs = []
    for c in chars:
        for seq in all_seqs:
            temp_seqs.append(seq)
            temp_seqs.append(str(seq[:pos-1] + c + seq[pos:]))
    #print(temp_seqs)
    all_seqs = temp_seqs.copy()

print(len(all_seqs))
#for seq in all_seqs:
    #print(seq)
'''
