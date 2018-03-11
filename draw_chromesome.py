#!/usr/bin/python3
import sys
import random
import pickle
import numpy as np
import matplotlib.pyplot as plt
import cv2

cmap = plt.get_cmap('jet')
with open('data/gene_id.pkl', 'rb') as f:
    id_to_gene = pickle.load(f)

with open('data/chrome.pkl', 'rb') as f:
    chromes, gene2chr = pickle.load(f)

DIV = 100000
'''
a = []
for k, v in chromes.items():
    if 'M' in k:
        continue
    print(k, v[1]-v[0])
    a.append(v[1] - v[0])

print(min(a))
print(max(a))
sys.exit(0)

 35146039   / 100000
249219373   

2500
'''

'''
tests = [
         (lambda x: x>0.99, [], "AUC > 0.99"),
         (lambda x: x>0.98, [], "AUC > 0.98"),
         (lambda x: x>0.95, [], "AUC > 0.95"),
         (lambda x: x<0.08, [], "AUC < 0.08"),
         (lambda x: x<0.05, [], "AUC < 0.05")]
         '''


#all_genes = []
count = {}
found = 0
miss = 0

images = {}
for k, v in chromes.items():
    images[k] = np.zeros((2500,), dtype=np.float32)
    pass

with open('missing.xxx', 'w') as mis:
    for l in open('auc/bl', 'r'):
        gene, score = l.strip().split('\t')
        score = float(score)
        #if score > 0.2: #< 0.95:
        #    continue
        loci = []
        names = id_to_gene[gene]
        for x in names.split('/'):
            x = x.strip()
            if len(x) == 0:
                continue
            if x in gene2chr:
                loci = gene2chr[x]
                break
            pass
        for chrome, start, end, strand in loci:
            chrome = chrome + strand
            b, e = chromes[chrome]
            start -= b
            end -= b
            start = start //DIV
            end = end//DIV
            images[chrome][start:end] = score
            pass
        pass

names = ['chr%d' % x for x in range(1,23)]
names.append('chrX')
names.append('chrY')

def gen_image (pos, pos_l, neg, neg_l):
    pos1 = cmap(pos, bytes=True)[:, :3][:,::-1]
    neg1 = cmap(neg, bytes=True)[:, :3][:,::-1]
    pos1[:,0][pos == 0] = 255
    pos1[:,1][pos == 0] = 255
    pos1[:,2][pos == 0] = 255
    pos = pos1
    neg1[:,0][neg == 0] = 255
    neg1[:,1][neg == 0] = 255
    neg1[:,2][neg == 0] = 255
    neg = neg1
    W = max(pos_l, neg_l)
    H = 64
    image = np.zeros((H, W, 3), dtype=np.uint8)
    for i in range(28):
        image[i, :, :] = pos[:W, :]
    for i in range(28):
        image[63-i, :, :] = neg[:W, :]
        pass
    return image

with open('chromes/index.html', 'w') as f:
    f.write('<HTML><BODY>\n')
    for name in names:
        pos = images[name + '+']
        b, e = chromes[name + '+']
        pos_l = (e-b)//DIV + 1
        neg = images[name + '-']
        b, e = chromes[name + '-']
        neg_l = (e-b)//DIV + 1
        f.write('%s<br/><img src="%s.png"></img><hr/>\n' % (name, name))
        cv2.imwrite('chromes/%s.png' % name, gen_image(pos, pos_l, neg, neg_l))
    f.write('</BODY></HTML>\n')


