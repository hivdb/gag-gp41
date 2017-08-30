#!/usr/bin/python
import csv
import sys
from io import StringIO

if len(sys.argv) != 3:
    print(
        "\nUsage: MEDSproc.py <resultsFile.csv> <alphaLevel> \n"
        "(alphaLevel is the p-value threshold before Bonferroni for 20 tests "
        "- recommended 0.01 or 0.05)"
    )
    exit()

file = open(sys.argv[1])
filetext = file.read()
csv_input_file = StringIO(filetext)

input_thresh = float(sys.argv[2])
thresh = input_thresh / 20  # bonferroni
aa_list = list("FLIMVSPTAYXHQNKDECWRG")
p_indices = [11 + 6 * i - 1 for i in range(21)]


def csvprint(strlist):
    print(",".join(map(str, strlist)))


csvReader = csv.reader(csv_input_file, delimiter=',', quotechar='|')
discard = next(csvReader)

MEDSheader = [
    "Site", "AA", "altL", "alt_p", "alt_omega", "AltFgNonSyn",
    "AltBgNonSyn", "AltSyn", "DivL", "DivFgNonSyn", "DivBgNonSyn", "DivSyn"
]

print("MEDS")
csvprint(MEDSheader)

for row in csvReader:
    p_row = [float(row[i]) for i in p_indices]
    w_row = [float(row[i + 1]) for i in p_indices]
    for i in range(20):
        if((p_row[i] < thresh) & (w_row[i] > 0.5)):
            # only want accelerated changes
            temp_list = [row[0]]
            temp_list.append(aa_list[i])
            # alt details
            temp_list = temp_list + \
                [row[p_indices[i] + j - 1] for j in range(2)]
            # transforming AA multiplier
            if (float(row[p_indices[i] + 1]) >= 1):
                temp_list.append("Infinity")
            else:
                temp_list.append((1 / (1 - float(row[p_indices[i] + 1]))) - 1)
            temp_list = temp_list + [row[p_indices[i] + j - 1]
                                     for j in range(3, 6)]
            # null(div) details
            temp_list.append(row[4])
            for k in row[6:9]:
                temp_list.append(k)
            csvprint(temp_list)

print()
print("FEEDS")

csv_input_file = StringIO(filetext)
csvReader2 = csv.reader(csv_input_file, delimiter=',', quotechar='|')
discard = next(csvReader2)
FEEDSheader = ["Site", "DivL", "p", "DivFgNonSyn",
               "DivBgNonSyn", "DivSyn", "NullL", "NullBgNonSyn", "NullSyn"]
csvprint(FEEDSheader)

for row in csvReader2:
    # only want positive selection
    if((float(row[5]) < input_thresh) & (float(row[6]) > float(row[8]))):
        temp_list = [row[0]]
        temp_list = temp_list + row[4:9] + row[1:4]
        csvprint(temp_list)
