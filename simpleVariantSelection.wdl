workflow SimpleVariantSelection {
  File gatk
  File refFasta
  File refIndex
  File refDict
  String name
  
  call haplotypeCaller {
    input:
      sampleName=name,
      RefFasta=refFasta,
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict
  }
  call select as selectSNPs {
    input:
      sampleName=name,
      RefFasta=refFasta,
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict,
      type="SNP",
      rawVCF=haplotypeCaller.rawVCF
  }
  call select as selectIndels {
    input:
      sampleName=name,
      RefFasta=refFasta,
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict,
      type="INDEL",
      rawVCF=haplotypeCaller.rawVCF
  }
}

task haplotypeCaller {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File inputBAM
  File bamIndex
  
  command {
    java -jar ${GATK} \
    HaplotypeCaller \
    -R ${RefFasta} \
    -I ${inputBAM} \
    -O ${sampleName}.raw.indels.snps.vcf
  }
  output {
    File rawVCF = "${sampleName}.raw.indels.snps.vcf"
  }
}

task select {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  String type
  File rawVCF

  command {
    java -jar ${GATK} \
      SelectVariants \
      -R ${RefFasta} \
      -V ${rawVCF} \
      -select-type ${type} \
      -O ${sampleName}_raw.${type}.vcf
  }
  output {
    File rawSubset = "${sampleName}_raw.${type}.vcf"
  }
}
