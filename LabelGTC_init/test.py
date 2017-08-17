#!/usr/bin/env python
from minSGT import getMinSGT
from ete3 import Tree


tlist = "((A1__A, C1__C),(B1__B, (C2__C,D2__D)));((B2__B,D3__D),(B1__B,C3__C))"
gtreelist = [Tree(x+";") for x in tlist.split(';')]
sptree = Tree("((A,B),(C,D));")
clades_to_preserve = "((A1__A, C1__C);(C2__C,D2__D);(B1__B,C3__C))"
treated_trees = "((C1__C, A1__A),(B1__B, (C2__C,D2__D)));((B1__B, (C2__C,D2__D));(B1__B,C3__C))"

gcontent = "".join([gt.write(format=9) for gt in gtreelist])
scontent = sptree.write(format=9).strip(';')

print getMinSGT(gcontent, scontent, True, clades_to_preserve, treated_trees, "")
