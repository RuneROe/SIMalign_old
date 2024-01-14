from pymol import cmd, stored
import sys
import numpy as np

# input = input("Input file: ")
input = sys.argv[1]
cmd.reinitialize()
stored.pLDDT = []
cmd.load(input)
pLDDT = []
cmd.select("sele",input.split("/")[-1].split(".")[0]+" and name CA")
# cmd.iterate("sele",pLDDT.append("b"))
cmd.iterate("sele","stored.pLDDT.append(b)")
# cmd.iterate("sele","pLDDT.append(b)")
pLDDT = np.array(stored.pLDDT)
print(input.split("/")[-1], np.mean(pLDDT))
