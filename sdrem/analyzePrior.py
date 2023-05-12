#!/usr/bin/env python
import pdb,sys,os
from junlib.File import *

priors=TabFile("nodePriors.txt").read("\t")
for i in priors:
	if float(i[1])>0.5:
		print(f"{i[0]}\tyes")
	else:
		print(f"{i[0]}\tno")
