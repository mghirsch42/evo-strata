
# Converts .rtf files to .csv files that can be used in the analysis

# file = "P_babaulti_LG19_XY_SNPs"
file = "P_curvifrons_LG19_XY_SNPs"
header_len = 7

data = []

with open(file+".rtf", "r") as f:
    lines = f.readlines()
    lines = lines[header_len:]
    lines = [l.replace("\n", "").replace("\t", ",").replace("}", "").replace("\\", "") for l in lines]
    lines[0] = lines[0][11:]

with open(file+".csv", "w") as f:
    f.write("id,start,stop,count\n")
    f.write("\n".join(lines))