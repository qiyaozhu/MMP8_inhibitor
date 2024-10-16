import os

with open("proc_0_pareto_analysis.txt", "r") as f:
    lines = f.readlines()

for i in range(5,110):
    fname = lines[i].split("_0001")[0]
    print("cp "+fname+"_0001_0001.pdb pareto_fronts/")
    os.system("cp "+fname+"_0001_0001.pdb pareto_fronts/")


