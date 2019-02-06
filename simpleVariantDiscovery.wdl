workflow SimpleVariantDiscovery {
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
  call hardFilterSNP {
    input:
      sampleName=name,
      RefFasta=refFasta,
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict,
      rawSNPs=selectSNPs.rawSubset
  }
  call hardFilterIndel {
    input:
      sampleName=name,
      RefFasta=refFasta,
      GATK=gatk,
      RefIndex=refIndex,
      RefDict=refDict,
      rawIndels=selectIndels.rawSubset
  }
  call combine {
    input:
      sampleName=name,
      GATK=gatk,
      filteredSNPs=hardFilterSNP.filteredSNPs,
      filteredIndels=hardFilterIndel.filteredIndels
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

task hardFilterSNP {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File rawSNPs

  command {
    java -jar ${GATK} \
      VariantFiltration \
      -R ${RefFasta} \
      -V ${rawSNPs} \
      --filter-expression "RS > 60.0" \
      --filter-name "snp_filter" \
      -O ${sampleName}.filtered.snps.vcf
  }
  output {
    File filteredSNPs = "${sampleName}.filtered.snps.vcf"
  }
}

task hardFilterIndel {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File rawIndels


  command {
    java -jar ${GATK} \
      VariantFiltration \
      -R ${RefFasta} \
      -V ${rawIndels} \
      --filter-expression "FS > 200.0" \
      --filter-name "indel_filter" \
      -O ${sampleName}.filtered.indels.vcf
  }
  output {
    File filteredIndels = "${sampleName}.filtered.indels.vcf"
  }
}

task combine {
  File GATK
  String sampleName
  File filteredSNPs
  File filteredIndels

  command {
    java -jar ${GATK} \
      MergeVcfs \
      -I ${filteredSNPs} \
      -I ${filteredIndels} \
      -O ${sampleName}.filtered.snps.indels.vcf
  }
  output {
    File filteredVCF = "${sampleName}.filtered.snps.indels.vcf"
  }
}


