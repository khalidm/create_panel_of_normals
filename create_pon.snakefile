
VERSION="0.1"

configfile: "cfg/config_pon.yaml"
cluster = json.load(open("cfg/cluster.json"))

# NOTE: no mitochondria MT because they aren't in our exome
GATK_CHROMOSOMES=('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')

def germline_samples():
  samples = set(config['samples'])
  tumours = set(config['tumours'])
  return list(samples.difference(tumours))

def germline_sample(wildcards):
  return config["tumours"][wildcards.tumour]

#ruleorder:  mutect2_sample_pon > mutect2_GenomicsDBImport > mutect2_CreateSomaticPanelOfNormals

### final outputs ###
rule all:
  input:
    expand("pon_db/{sample}.mutect2.pon.vcf.gz", sample=config['samples']),
    #expand("pon_db/{chromosome}", chromosome=GATK_CHROMOSOMES),
    expand("pon_db/{chromosome}.completed", chromosome=GATK_CHROMOSOMES),
    expand("pon_db/mutect2.pon.{chromosome}.vcf.gz", chromosome=GATK_CHROMOSOMES),

rule mutect2_sample_pon:
  input:
    reference=config["genome"],
    bam="out/{germline}.sorted.dups.bam",
    regions=config["regions"],
  output:
    "pon_db/{germline}.mutect2.pon.vcf.gz",
  log:
    stderr="log/{germline}.mutect2.pon.stderr"
  shell:
    "{config[module_java]} && "
    "mkdir -p pon_db && "
    "tools/gatk-4.1.2.0/gatk Mutect2 -R {input.reference} -I {input.bam} -O {output} --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter --max-mnp-distance 0 2>{log.stderr}"

rule mutect2_GenomicsDBImport:
  input:
    reference=config["genome"],
    regions_chr=config["regions_name"] + "_{chromosome}.bed",
    #vcfs=rules.mutect2_sample_pon.output
    vcfs=expand("pon_db/{germline}.mutect2.pon.vcf.gz", germline=germline_samples()),
  output:
    comp="pon_db/{chromosome}.completed",
  params:
    vcfs_p=' '.join(['-V {}'.format(vcf) for vcf in expand("pon_db/{germline}.mutect2.pon.vcf.gz", germline=germline_samples())]),
    dir="pon_db/{chromosome}",
    cores=cluster["mutect2_GenomicsDBImport"]["n"]
  log:
    stderr="log/mutect2.gdb.{chromosome}.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.1.2.0/gatk --java-options '-Xmx12g -Xms12g' GenomicsDBImport -R {input.reference} -L {input.regions_chr} --reader-threads {params.cores} --genomicsdb-workspace-path {params.dir} {params.vcfs_p} 2>{log.stderr} && "
    "touch {output[0]}"

## Create a pon VCF per chromosome
rule mutect2_CreateSomaticPanelOfNormals:
  input:
    reference=config["genome"],
    pon_db="pon_db/{chromosome}.completed",
  output:
    "pon_db/mutect2.pon.{chromosome}.vcf.gz",
  log:
    stderr="log/mutect2.gdb.vcf.{chromosome}.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.1.2.0/gatk --java-options '-Xmx12g -Xms12g' CreateSomaticPanelOfNormals -R {input.reference} -V gendb://pon_db/{wildcards.chromosome} -O {output} 2>{log.stderr}"
