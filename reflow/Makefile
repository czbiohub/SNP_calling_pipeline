# Makefile
# authors: Olga Botvinnik & Lincoln Harris
# date: 10.3.18

# test_gatk -- run GATK reflow pipeline for a single cell. 
#		need to specify full s3 paths here
test_gatk:
	reflow run gatk_pipeline.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/nonImmune_bams_9.27/180226/C3_1001000403/C3_1001000403_S123.homo.Aligned.out.sorted.bam \
		-sample_id=C3_1001000403 \
		-output_vcf=s3://darmanis-group/singlecell_lungadeno/nonImmune_bams_9.27/180226/C3_1001000403/C3_1001000403.vcf

# test_batch -- trying out a reflow batch job
#		correct command for this? 
test_batch: 
	reflow runbatch config.json

