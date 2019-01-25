# Makefile
# authors: Lincoln Harris
# date: 1.23.18

# test_gatk -- run picard tools reflow pipeline for a single cell. 
#		need to specify full s3 paths here
test:
	reflow run picard_markDuplicates.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/170125/G12_1001000292/G12_1001000292_S72.homo.Aligned.out.sorted.bam \
		-sample_id=G12_1001000292 \
		-output_file=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/removeDups_out/G12_1001000292/G12_1001000292_marked_dups.bam

test1:
	reflow run picard_markDuplicates.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/180226/C3_1001000403/C3_1001000403_S123.homo.Aligned.out.sorted.bam \
		-sample_id=C3_1001000403 \
		-output_file=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/removeDups_out/C3_1001000403/C3_1001000403_marked_dups.bam

test2:
	reflow run picard_markDuplicates.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/170508/B5_1001000327/B5_1001000327_S259.homo.Aligned.out.sorted.bam \
		-sample_id=B5_1001000327 \
		-output_file=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/removeDups_out/B5_1001000327/B5_1001000327_marked_dups.bam

test3:
	reflow run picard_markDuplicates.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/170508/C3_1001000327/C3_1001000327_S295.homo.Aligned.out.sorted.bam \
		-sample_id=C3_1001000327 \
		-output_file=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/removeDups_out/C3_1001000327_marked_dups.bam

test4:
	reflow run picard_markDuplicates.rf \
		-input_bam=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/180604/O17_B002090/O17_B002090_S353.homo.Aligned.out.sorted.bam \
		-sample_id=O17_B002090 \
		-output_file=s3://darmanis-group/singlecell_lungadeno/non_immune/nonImmune_bams_9.27/removeDups_out/O17_B002090_marked_dups.bam
