# pyGenicCT

Clumping via plink and then p value thresholding through python

## Instructions

This will give show how i have created clumped PRS via plink and python. It is assuming you are working with systems similar to the IEU in bristol and have access to BlueCrystal 4, but should be adjustable to other systems and locations with some minor adjustment. 

### Getting and formatting summary statistics

First you will need some summary statistics in a tsv file format for plink. If you have data from the GWAS catalog it should work as is, where as if you have data from the IEU OpenGWAS it will be in a .vcf file. Plink won't like using that as a summary statistic file, so we need to convert it. 

We can do this within the pyGenicParser package, from loading the VCFObject object. Provide VCF object the path to your vcf file, it can be gzipped VCFObject will read both zipped and non zipped files. If using the OpenGWAS data you will have some headings you don't need. Whilst you do not need to filter them out you can do so if you choose. Then, covert it to a tsv via convert_to_summary, by providing the output path and a file name. OpenGWAS files also store the p values in Log P, where as plink is expecting p values. We can convert them back by giving the file header to log_p_convert. 


```python
from pyGenicParser import VCFObject

path_to_vcf = r"Path_Here"

vcf = VCFObject(path_to_vcf)

# If your using the OpenGWAS data you will have these headers that we don't need
for hv in ["QUAL", "FILTER", "SS_format_0", "EZ_format_0", "SI_format_0", "NC_format_0", "ID_format_0"]:
    vcf.write_headers[hv] = False

# This will give Where ID is the variatn ID, LP is the log p value, and ES is the estimated coefficient.
print([hv for hv, ht in vcf.write_headers.items() if ht])
>> ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AF_info', 'ES_format_0', 'SE_format_0', 'LP_format_0', 'AF_format_0']

vcf.covert_to_summary(r"Output_Path", "File_name", log_p_convert="LP_format_0")

```
### Clumping via plink

The basic commands you need to run the plink file are as follows, this runs a conservative clump based on the parameters for Clumping of twoSampleMR. Here we are using the same LD reference as used by the [IEU OpenGWAS by default](https://github.com/mrcieu/gwasvcf/) which can be downloaded at this github page as the ***1000 Genomes reference panels for LD for each super population***.  If using the OpenGWAS, the only variable you will **need** to change is the path to this LD preference (--bfile) and the name of the summary statistic file you created from the VCFObject (--clump). 

This is also assuming you are using the Slurm method of job submission. Should you be using Bluecrystal phase 3 which uses the PBS job system see [this conversion](https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html#slurm-for-speakers-of-other-scheduler-languages), if using something else you may be able to convert it via [this conversion resource](https://slurm.schedmd.com/rosetta.html). 

```bash 

#!/bin/bash

# Generated by pyGenicCT Version: 0.01.0
#SBATCH --job-name=Clump
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem=10000M

module load apps/plink/1.90

plink \
	--bfile LDReference/EUR \
	--clump-p1 1 \
	--clump-p2 1 \
	--clump-r2 0.001 \
	--clump-kb 10000 \
	--clump File_Name.tsv.gz\
	--clump-snp-field ID \
	--clump-field LP_format_0 \
	--out EUR_CLUMP_Edu \
	--threads 1
```

Once this is completed you will get a few files from plink, which we then want to extract the snps that we have clumped on out off. To do this we can use another script with awk. This will unzipped and then re-zip the data as required. 

```bash
#!/bin/bash

# Generated by pyGenicCT Version: 0.01.0
#SBATCH --job-name=ExtractResults
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --mem=10000M

awk '{print $1,$3}' EUR_CLUMP_Edu.clumped > SNP.valid
gzip -d file_name.tsv.gz

awk '{print $3,$7,$9}' file_name.tsv > SNP.values
gzip Okbay.tsv

```

You can also generate these scripts in python should you wish via the following


```python
from pyGenicCT import plink_clump, extract_results

plink_clump(r"Local script save location",
            "Clump",
            "apps/plink/1.90",
            "/LDReference/EUR",
            "SummaryFile.tsv",
            "EUR_CLUMP",
            snp_field="ID",
            p_field="LP_format_0"
            )

extract_results(r"Local script save location",
                "ExtractResults",
                "SummaryFile.tsv",
                "EUR_CLUMP.clumped",
                3, 7, 9)

```

### Creating the PRS

Now we want to link the snps we found to the summary data so that we can use it to create the PRS. All of these scripts are in python, and call the CTScores object which takes a yaml file as its only parameter. Since we are going to need to set it up later, lets fully define it now.

```yaml
chromosome_index: 0
snp_index: 1
coefficient_index: 2
p_value_index: 3

Valid: /SNP.valid
Values: /SNP.values
valid_snp_name: SNP
values_snp_name: variant_id

meta_path: pysnptools Meta Path
gen_path: Genetic path
write_path: working directory
write_name: PRSScore
base_name: data.chr

threshold: [0.00000005, 0.00005, 0.005, 0.05, 1]
```

The first four variables should not need to be changed unless a version adds additional variables to link_resources(). Valid and Values are the paths to the files you created via extract_results, with the the column headers of the variant name for each file being given after that.

This package uses pysnptools, which creates metadata. You need to provide a writable path for this metadata if you do not have write permissions to the location of the genetic files. **This pipeline is currently only configured to work with Bgen data.** You then need to provide the directory that holds the chromosome data, set the write name of the final output, and the base name of the genetic files; for example **data.chr**02.bgen. Finally you need to set the P value thresholds you want to iterate though.

Now using this yaml file as a soul argument, submit a python script with the following below via sbatch and this will link the resources for you ready for use. 

```python
from CTScores import CTScores
import sys

print(f"submitting {sys.argv[1]}")

CTScores(sys.argv[1]).link_resources()

```

This will make a file called Snps.csv in your working directory. All you need to do now is run the script with
create_score_levels and this will now create the p_value thresholds PRS of the levels you defined in threshold within the .yaml file as a csv called whatever the write_name was set as. Keep in mind, this assumes that the IID's are embedded into the Bgen files, as the IIDs will be extract from within it rather than via .sample. This may be extended to allow for this in future. 

```python
from CTScores import CTScores
import sys

print(f"submitting {sys.argv[1]}")

CTScores(sys.argv[1]).create_score_levels()

```
