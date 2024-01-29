import sys
import numpy as np
import math 


if len(sys.argv) < 4:
    print("\nUSAGE:\nprogram.py importfile.csv outfile.txt target_line_index")
    print("\nEXAMPLE:\nprogram.py book.csv out.txt 2730\nThis will compare all reads to GPA1 which have index value 2730")
    sys.exit(1)
infile = open(sys.argv[1],"r")
outfile = open(sys.argv[2],"w")
lines = infile.readlines()
target_index = int(sys.argv[3])
out = "Target: "+lines[target_index].split(";")[0]+"\n\n"
print(out)
outfile.write(out)
target = lines[target_index].split(";")[1:5]  #  <------ Takes only WT levels
target = [np.log2(float(x)) for x in target]  #  <------ log2 transform target

score = []
max_av = 0
index_list = []
for j, line in enumerate(lines):
    values = line.split(";")[1:5]   #  <------ Takes only WT levels
    values = [np.log2(float(x)) for x in values]  #  <------ log2 transform target
    av = np.mean(np.array(values))
    diff = []
    for i, v in enumerate(values):
        diff.append(v/target[i])
    sd = np.std(np.array(diff))
    if sd > 1:
        sd = 0
    else:
        sd = 1-sd
    score.append([av,sd])
    if av > max_av:
        max_av = av
    index_list.append(j)

score_norm = []
for val in score:
    number = (1/4)*(val[0]/max_av)+(3/4)*(val[1])  #  <------ 50:50 magnitute and similarity to target
    # number = val[1]
    if math.isnan(number):
        score_norm.append(0)
    else:
        score_norm.append(number)

sort = sorted(index_list, key=lambda x: score_norm[x], reverse=True)


for i in range(len(score_norm)):
    outfile.write(lines[sort[i]].split(";")[0]+"\t"+"\n")
outfile.close()
infile.close()