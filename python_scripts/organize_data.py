import pandas as pd
import os
import numpy as np

project_path = "/n/data1/hms/dbmi/farhat/jiazheng"
eqtl_directory = f"{project_path}/eqtl"
metadata_directory = f"{eqtl_directory}/metadata"
dna_results_path = f"{eqtl_directory}/dna_analysis"
rna_results_path = f"{eqtl_directory}/rna_analysis/align_results"

# metadata = pd.read_csv(f"{eqtl_directory}/metadata/eqtl_metadata.csv")
# metadata["Lineage"] = metadata["Lineage"].apply(lambda x: x.split('.')[0])

### Expression
rna_samples = os.listdir(rna_results_path)

i = 0
counts = pd.read_csv(f"{rna_results_path}/{rna_samples[i]}/{rna_samples[i]}_counts.txt", 
                     skiprows=2, sep='\t', 
                     names=["Geneid","Chr","Start","End","Strand","Length","Counts"])
counts[rna_samples[i]] = counts["Counts"] / counts["Length"]
counts[rna_samples[i]] = 1e9 * counts[rna_samples[i]] / counts["Counts"].sum()
counts[rna_samples[i]] = counts[rna_samples[i]].apply(np.log)
counts = counts[["Geneid", rna_samples[i]]]

for i in range(1, len(rna_samples)):
    tmp = pd.read_csv(f"{rna_results_path}/{rna_samples[i]}/{rna_samples[i]}_counts.txt", 
                     skiprows=2, sep='\t', 
                     names=["Geneid","Chr","Start","End","Strand","Length","Counts"])
    tmp[rna_samples[i]] = tmp["Counts"] / tmp["Length"]
    tmp[rna_samples[i]] = 1e9 * tmp[rna_samples[i]] / tmp["Counts"].sum()
    tmp[rna_samples[i]] = tmp[rna_samples[i]].apply(np.log)
    tmp = tmp[["Geneid", rna_samples[i]]]
    counts = counts.merge(tmp, how='left', on=["Geneid"])

counts = counts.set_index("Geneid")
counts = counts.transpose()
counts = counts.reset_index(names=["RNA_Accession"])

rna_metadata = pd.concat([pd.read_csv(f"{metadata_directory}/jody_rna.csv"), 
                          pd.read_csv(f"{metadata_directory}/chiner_oms_rna.csv")], axis=0)

counts = rna_metadata.merge(counts, how="left", on="RNA_Accession")
counts.to_csv(f"{eqtl_directory}/results/GeneExpressionLogFPKM.csv", index=False)
del counts

# Variants
sample_names = np.array(os.listdir(dna_results_path))

variants = []
failed = []
for sample_name in sample_names:
    input_filename = f"{dna_results_path}/{sample_name}/{sample_name}_filtered.vcf"
    if not os.path.isfile(input_filename):
        failed.append(sample_name)
        continue
    vcf = pd.read_csv(input_filename, skiprows=34, sep='\t', names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER',
                                                                    'INFO','FORMAT','SAMPLE'])
    VarID = vcf['REF'] + "_" + vcf['POS'].astype(str) + "_" + vcf['ALT']
    VarID = VarID.to_frame(name='VarID')
    VarID['POS'] = vcf['POS']
    VarID['FILTER'] = vcf['FILTER']
    VarID['SAMPLE'] = sample_name
    variants.append(VarID)

print(failed)
variants = pd.concat(variants)
variants['COUNT'] = [1 if Filter=='PASS' else 0 for Filter in variants['FILTER'].values.flatten()]
variants['VarID'].unique().shape


variants = variants.pivot_table(index="VarID", columns='SAMPLE', values='COUNT', aggfunc="max")
variants.shape

variants = pd.DataFrame(sample_names, columns=['SAMPLE']).merge(variants.transpose().reset_index(), how='left', on=['SAMPLE'])
variants = variants.set_index("SAMPLE").transpose()

variants = variants.replace(np.nan, 0).astype('int')
variants.shape
# variants = variants.iloc[np.where((variants.mean(axis=1)>0.05) & (variants.mean(axis=1)<0.95))[0], :].transpose()
indexers = np.where((variants.min(axis=1)==0) & (variants.max(axis=1)==1))
variants = variants.iloc[indexers[0], :].transpose()
variants.shape
variants = variants.reset_index(names=["DNA_Accession"])

dna_metadata = pd.concat([pd.read_csv(f"{metadata_directory}/jody_dna.csv"), 
                          pd.read_csv(f"{metadata_directory}/chiner_oms_dna.csv")], axis=0)

variants = dna_metadata.merge(variants, how='left', on="DNA_Accession")
variants.to_csv(f"{eqtl_directory}/results/Variants.csv", index=False)
del variants

### Region of Difference
def depth_extractor(coordinates, data_type="dna", prefix="ERR"):
    [lower, upper] = coordinates
    result_path = f"{eqtl_directory}/{data_type}_analysis"
    norm_coef = 4411532 if data_type=='dna' else 1e6
    sample_names = np.array(os.listdir(result_path))
    sample_names = sample_names[[prefix in sample_name for sample_name in sample_names]]
    depths = []
    for sample_name in sample_names:
        filename = f"{result_path}/{sample_name}/{sample_name}_coverage.txt"
        coverage = pd.read_csv(filename, header=None, sep='\t', quotechar='"')
        coverage = coverage.iloc[:, 0:3]
        coverage.columns = ["ref", "pos", "depth"]
        coverage["depth"] = coverage["depth"].div(coverage["depth"].sum()) * norm_coef
        depth = coverage.query("pos>=@lower&pos<=@upper")[["depth"]]
        depth.columns=[sample_name]
        depths.append(depth)
    depths.append(coverage.query("pos>=@lower&pos<=@upper")["pos"])
    depths = pd.concat(depths, axis=1)
    depths = depths.reset_index(drop=True)
    return depths

RD236a_coordinates = [4056945, 4058396]
RD236a_depths = depth_extractor(RD236a_coordinates, data_type="dna")
RD236a = RD236a_depths.drop(columns=["pos"], inplace=False).mean().to_frame(name='RD236aMeanDepth').reset_index(names="DNA_Accession")
RD236a = dna_metadata.merge(RD236a, how='left', on='DNA_Accession')
RD236a['HasDeletion'] = ["No Deletion" if coverage > 0.1 else "With Deletion" for coverage in RD236a["RD236aMeanDepth"] ]
RD236a.head()
RD236a.to_csv(f"{eqtl_directory}/results/RD236a_MeanDepth.csv", index=False)