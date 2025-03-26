import glob
import csv
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("raw")
parser.add_argument("sample_path")
args = parser.parse_args()

raw = args.raw
sample_path = args.sample_path

# pipeline = {"BWA_MEM.csv":None,
#             "FASTQC.csv":None,
#             "MKDUP.csv":"BWA_MEM.csv",
#             "RECAL.csv":"MKDUP.csv",
#             "SPLIT_BAM.csv":"RECAL.csv",
#             "STRELKA.csv":"RECAL.csv",
#             "SAGE.csv":"RECAL.csv",
#             "MuSE2.csv":"RECAL.csv",
#             "CONPAIR.csv":"RECAL.csv",
#             "MOSDEPTH.csv":"RECAL.csv",
#             "Mutect2.csv":"SPLIT_BAM.csv",
#             }

pipeline = {"BWA_MEM.csv":None,
            "FASTQC.csv":None,
            "MKDUP.csv":"BWA_MEM.csv",
            "RECAL.csv":"MKDUP.csv",
            "STRELKA.csv":"RECAL.csv",
            "SAGE.csv":"RECAL.csv",
            "MuSE2.csv":"RECAL.csv",
            "CONPAIR.csv":"RECAL.csv",
            "MOSDEPTH.csv":"RECAL.csv",
            }

results = {}

with open(sample_path, newline='') as all_samples:
  sample_dt = pd.read_csv(all_samples)
  sample_dt['ID'] = sample_dt['patient'] + "_" + sample_dt['status']
  sample_num = sample_dt.shape[0]
  normal_tumor_pair_num = sample_dt[sample_dt['status'] == 'tumor'].shape[0]
print("Total sample number:" + str(sample_dt.shape[0]))
print("Total patient number:" + str(sample_dt.groupby('patient').ngroups))
print("Total normal tumor pair number:" + str(sample_dt.groupby('patient').ngroups))

expected = {"BWA_MEM.csv":sample_num,
            "FASTQC.csv":sample_num,
            "MKDUP.csv":sample_num,
            "RECAL.csv":sample_num,
            "STRELKA.csv":normal_tumor_pair_num,
            "SAGE.csv":normal_tumor_pair_num,
            "MuSE2.csv":normal_tumor_pair_num,
            "CONPAIR.csv":normal_tumor_pair_num,
            "MOSDEPTH.csv":normal_tumor_pair_num,
            "Mutect2.csv":normal_tumor_pair_num,
            }
difference = {}

files = raw.strip("[]").split(',')

for file in files:
  with open(file.lstrip(), newline='') as csvfile:
    file_name = file.split('/')[-1]
    dt = pd.read_csv(csvfile)
    if file_name in ["BWA_MEM.csv","FASTQC.csv","MKDUP.csv","BWA_MEM.csv","RECAL.csv","MKDUP.csv","RECAL.csv"]:
        dt['ID'] = dt['patient'] + "_" + dt['status']
    results[file_name] = dt

for p in pipeline:
  print("Process: "+p[:-4])
  print("Expected: "+str(expected[p]))
  print("Got: "+str(results[p].shape[0]))
  if expected[p] != results[p].shape[0]:
    if p in results and expected[p] != results[p].shape[0]:
      difference[p] = sample_dt[~sample_dt.ID.isin(results[p].ID)]
    else:
      difference[p] = sample_dt[~sample_dt.patient.isin(results[p].patient)]


print("Summary")
print("Missing Files:")
if not difference:
  print("None")
else:
  for p in difference:
    print("Process: "+p[:-4])
    print(difference[p])
  print("Recommended command:")
  print("nextflow run main.nf --sample rerun.csv -with-conda")

  pd.concat([dt for dt in difference.values()]).drop_duplicates().to_csv("rerun.csv",index=False)