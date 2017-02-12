#!/usr/bin/env python

from math import sqrt
from os import remove
from sys import argv

fmean = open(argv[1])
B = fmean.readlines()
fmean.close()

fvar = open(argv[2])
C = fvar.readlines()
fvar.close()

def numerify(s):
  D = len(s) - 2
  l = [int(s[0])]
  for d in range(D):
    l.append(float(s[d + 1]))
  l.append(float(s[-1]))
  return l

B = [numerify(b.strip().split()) for b in B]
C = [numerify(c.strip().split()) for c in C]

ks = list(set([b[0] for b in B]))

for k in ks:
  T    = [b[1:-1] for b in B if b[0] == k]
  mean = [b[-1]   for b in B if b[0] == k]
  var  = [c[-1]   for c in C if c[0] == k]
  std  = [sqrt(max(0, v)) for v in var]

  Z = zip(T, mean, std)
  dat = [z[0] + [z[1], z[1]-z[2], z[1]+z[2]] for z in Z]
  dstr = '\n'.join([' '.join(map(str, d)) for d in dat])

  fgp = open(argv[3] + '.k{}.dat'.format(k), 'w')
  fgp.write(dstr)
  fgp.close()

