# Makefile
# authors: Olga Botvinnik & Lincoln Harris
# date: 10.3.18

test:
	reflow run gatk_pipeline.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/170510/A14_1001000339/A14_1001000339_S64.homo.Aligned.out.sorted.bam \
		-sample_id=A14_1001000339 \
		-output_vcf=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/vcf/

test_gvcf:
	reflow run gatk_pipeline_gvcf.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/170510/J21_1001000339/J21_1001000339_S133.homo.Aligned.out.sorted.bam \
		-sample_id=J21_1001000339 \
		-output_vcf=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/gvcf/J21_1001000339.g.vcf

test_gvcf1:
	reflow run gatk_pipeline_gvcf.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/180423/D12_B003528/D12_B003528_S24.homo.Aligned.out.sorted.bam \
		-sample_id=D12_B003528 \
		-output_vcf=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/gvcf/D12_B003528.g.vcf

test_gvcf2:
	reflow run gatk_pipeline_gvcf.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/181120/N6_B001016/N6_B001016_S270.homo.Aligned.out.sorted.bam \
		-sample_id=N6_B001016 \
		-output_vcf=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/gvcf/