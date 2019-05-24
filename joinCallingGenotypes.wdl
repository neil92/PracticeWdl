workflow jointCallingGenotypes {
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File gatk
  File refFasta
  File refIndex
  File refDict

  scatter (sample in inputSamples) {
    call HaplotypeCallerERC {
      input: GATK=gatk,
        RefFasta=refFasta,
        RefIndex=refIndex,
        RefDict=refDict,
        sampleName=sample[0],
        bamFile=sample[1],
        bamIndex=sample[2]
    }
  }
  call CombineGVCFs {
    input: GATK=gatk,
    RefFasta=refFasta,
    RefIndex=refIndex,
    RefDict=refDict,
    sampleName="CEUtrio",
    GVCFs=HaplotypeCallerERC.GVCF
  }
  call GenotypeGVCFs {
    input: GATK=gatk,
    RefFasta=refFasta,
    RefIndex=refIndex,
    RefDict=refDict,
    sampleName=CombineGVCFs.outputSampleName,
    GVCF=CombineGVCFs.CombinedGVCF
  }
}

task HaplotypeCallerERC {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File bamFile
  File bamIndex

  command {
    java -jar ${GATK} \
        HaplotypeCaller \
        -ERC GVCF \
        -R ${RefFasta} \
        -I ${bamFile} \
        -O ${sampleName}_rawLikelihoods.g.vcf
  }
  output {
    File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
  }
}
task CombineGVCFs {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Array[File] GVCFs
  
  command {
    java -jar ${GATK} \
        CombineGVCFs \
        -R ${RefFasta} \
        --variant ${sep=" --variant " GVCFs} \
        -O "${sampleName}_rawLikelihoods.g.vcf"
  }
  output {
    File CombinedGVCF = "${sampleName}_rawLikelihoods.g.vcf"
    String outputSampleName = "${sampleName}"
  }
}
task GenotypeGVCFs {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  File GVCF

  command {
    java -jar ${GATK} \
        GenotypeGVCFs \
        -R ${RefFasta} \
        -V ${GVCF} \
        -O ${sampleName}_rawVariants.vcf
  }
  output {
    File rawVCF = "${sampleName}_rawVariants.vcf"
  }
}
