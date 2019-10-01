
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

### final outputs ###
rule all:
  input:
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

# Create a GenomicsDB from the normal Mutect2 calls
rule mutect2_GenomicsDBImport:
  input:
    reference=config["genome"],
    regions=config["regions"],
    vcfs=expand("pon_db/{germline}.mutect2.pon.vcf.gz", germline=germline_samples()),
  output:
    directory("pon_db/{chromosome}")
  params:
    vcfs=' '.join(['-V {}'.format(vcf) for vcf in expand("pon_db/{germline}.mutect2.pon.vcf.gz", germline=germline_samples())]),
    cores=cluster["mutect2_GenomicsDBImport"]["n"]
  log:
    stderr="log/mutect2.gdb.{chromosome}.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.1.2.0/gatk --java-options '-Xmx12g -Xms12g' GenomicsDBImport -R {input.reference} -L {wildcards.chromosome} --reader-threads {params.cores} --genomicsdb-workspace-path {output} {params.vcfs} 2>{log.stderr}"

# Create a pon VCF per chromosome
rule mutect2_GenomicsDBImport_vcf:
  input:
    reference=config["genome"],
    pon_db="pon_db/{chromosome}",
  output:
    "pon_db/mutect2.pon.{chromosome}.vcf.gz"
  params:
    vcfs=' '.join(['-V {}'.format(vcf) for vcf in expand("pon_db/{germline}.mutect2.pon.vcf.gz", germline=germline_samples())]),
    cores=cluster["mutect2_GenomicsDBImport"]["n"]
  log:
    stderr="log/mutect2.gdb.vcf.{chromosome}.stderr"
  shell:
    "{config[module_java]} && "
    "tools/gatk-4.1.2.0/gatk --java-options '-Xmx12g -Xms12g' CreateSomaticPanelOfNormals -R {input.reference} -V gendb://out/pon_db/{wildcards.chromosome} -O {output} 2>{log.stderr}"
