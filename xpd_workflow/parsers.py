import os
import numpy as np


def parse_xrd_standard(filename):
    fn = os.path.abspath(filename)
    with open(fn, 'r') as f:
        data = f.read().split('n^2')[-1]
    tthl = []
    il = []
    hkll = []
    for line in data.split('\n'):
        if line.split():
            tth, _, i, _, h, k, l, _, _, _ = line.split()
            l = l.strip(')')
            tthl.append(float(tth))
            il.append(float(i))
            hkll.append('{},{},{}'.format(h, k, l))
    return np.asarray(tthl), np.asarray(il), hkll

if __name__ == '__main__':
    a = '/mnt/bulk-data/research_data/Pr2NiO4orthorhombicPDF#97-008-1577.txt'
    print(parse_xrd_standard(a))
