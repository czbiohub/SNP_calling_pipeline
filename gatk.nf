#!/usr/bin/env netflow

params.sample_id = "sample_id"
params.input_bam = "$baseDir/data/bam/${params.sample_id}.bam"
params.output_vcf = "$baseDir/data/variants/${params.sample_id}.vcf.gz"
params.threads = 16
params.genome = "hg38"
params.genome_suffix = "-plus"

genome_tar = file("s3://czi-hca/ref-genome/${params.genome}${params.genome_suffix}.tgz")

process decompress {
    input:
	file genome from genome_tar

    output:
	file "${params.genome}${params.genome_suffix}/${params.genome}${params.genome_suffix}.fa" into genome_file

    """
    tar xvf ${genome}
    """
}

process samtools_fasta_index {
    input:
	file genome from genome_file

    output:
	file "${genome}.fai" into genome_index

    """
    samtools faidx ${genome}
    """
}

process make_genome_dictionary {
    input:
	file genome from genome_file

    output:
	file "${genome}.dict" into genome_dict

    """
    gatk CreateSequenceDictionary --REFERENCE ${genome} --OUTPUT ${genome}.dict
    """
}

// process add_or_replace_groups {
//     input:
//	sorted bam file

//     output:
//	another sorted bam file?

//     """
//     gatk AddOrReplaceReadGroups -I {{sorted_bam}} -O {{output}} -RGID 4 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM 20
//     """
// }

// process samtools_index {
//     input:
//	group bam?

//     output:
//	samtools index files

//     """
//     samtools index {{group_bam}} {{output}}
//     """
// }

// process generate_dbSNP_index {
//     """
//     gatk IndexFeatureFile -F {{dbSNP}} -O {{output}}
//     """
// }

// process haplotype_caller {
//     """
//     gatk HaplotypeCaller -R {{d}}/{{genome}}.fa -O {{output}} \
//			 -I {{d}}/{{sample_id}}.bam \
//			 --disable-read-filter MappingQualityReadFilter \
//			 --disable-read-filter GoodCigarReadFilter \
//			 --disable-read-filter NotSecondaryAlignmentReadFilter \
//			 --disable-read-filter MappedReadFilter \
//			 --disable-read-filter MappingQualityAvailableReadFilter \
//			 --disable-read-filter NonZeroReferenceLengthAlignmentReadFilter \
//			 --disable-read-filter NotDuplicateReadFilter \
//			 --disable-read-filter PassesVendorQualityCheckReadFilter \
//			 --disable-read-filter WellformedReadFilter \
//			 --native-pair-hmm-threads {{threads}} \
//			 --dbsnp {{d}}/{{genome}}.vcf
//     """
// }
