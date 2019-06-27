#!/usr/bin/env python

import argparse
import sys, os, re
import logging
import subprocess
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)


##########################
# Add options to inputs
##########################

#------------
## Required
#------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description = "Run VADiR.")
parser.add_argument('--ID', type = str , required = True, help="ID")
parser.add_argument('--RNA', type = str , required = True, help="RNA directory")
parser.add_argument('--out_dir', type = str , required = True, help="output directry")
parser.add_argument('--vadir', type = str , required = True, help="vadir directry")
parser.add_argument('--config', type = str , required = True, help="config file")

#-----------
## Optional 
#-----------
parser.add_argument('--skip_refine', required = False, default = False, action = "store_true", help="Run refine")
parser.add_argument('--skip_star', required = False, default = False, action = "store_true", help="Run STAR")
parser.add_argument('--skip_snpir', required = False, default = False, action = "store_true", help="Run SNPiR if True")
parser.add_argument('--skip_var', required = False, default = False, action = "store_true", help="Run variant calling in SMPiR")
parser.add_argument('--skip_mutect', required = False, default = False, action = "store_true", help="Run MuTect2 step.")
parser.add_argument('--skip_VADiR', required = False, default = False, action = "store_true", help="Run MuTect2 step.")

parser.add_argument('--threads', type = str , required = False, default = 4, help="Threads")
parser.add_argument('--vcf', type = str , required = False, default = None, help="vcf input file. will not use UG output.")
parser.add_argument('--bam', type = str , required = False, default = None, help="refined bam input file. will not use vadir refined bam.")
args = parser.parse_args()

#-------------------------
# parse arguments
#-------------------------
ID = args.ID
out_dir = args.out_dir
vadir = args.vadir
CONFIG = args.config
RNA = args.RNA
skip_star = args.skip_star
skip_refine = args.skip_refine
skip_var = args.skip_var
skip_mutect = args.skip_mutect
vcf = args.vcf
bam = args.bam
threads = args.threads
skip_VADiR = args.skip_VADiR
skip_snpir = args.skip_snpir
print(args)

###############
# Constants 
###############
ref_dir = "{}/refs".format(vadir)

# Refine
ref_genome_fasta = ref_dir + "/ucsc.hg19.fasta"

# Unified Genotyper
coding_bed = ref_dir + "/coding.bed"
dbsnp = ref_dir + "/dbsnp_138.hg19.vcf"

# Tools 
bedtools = "{}/tools/bedtools-2.25.0".format(vadir)
# pblat = "/broad/software/free/Linux/redhat_6_x86_64/pkgs/blat_35/bin/x86_64/blat"
pblat = "/seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/VaDiR/tools/pblat/pblat"

if vcf == None:
    ug_out = "{}/unifiedgenotyper/snpir/{}.ug.raw.vcf".format(out_dir, ID)
else:
    ug_out = vcf

if bam == None:
    refined_bam =  "{}/refined-mapping/{}.refined.bam".format(out_dir, ID)
else:
    refined_bam = bam

# TEST sanity check 
print("reference dir: ", ref_dir)
print("reference genome: ", ref_genome_fasta)
print("Configuration: ", CONFIG)
print("Raw VCF: ", ug_out)
print("Refined bam: ", refined_bam)
print("skip star, skip refine, skip variation : ", skip_star, skip_refine, skip_var)

###############
# cd to vadir
###############
cd = "cd {}".format(vadir)
logging.info(cd)
os.chdir(vadir)


logging.info("\n STEP: SETUP and CONGIF \n")

cmd_config = "/ahg/regevdata/projects/trinity/mbrown/VaDiR/SETUP.sh {}".format(CONFIG)
logging.info(cmd_config)
subprocess.call(cmd_config, shell = True)

# Make output directrory 
os.system("mkdir -p {}".format(out_dir))
os.system("mkdir -p {}/unifiedgenotyper/snpir".format(out_dir))
os.system("mkdir -p {}/snpir".format(out_dir))

##################################################################################################################################################################
##################################################################################################################################################################
# STEP 1
# Alignment of RNA 
##################################################################################################################################################################
##################################################################################################################################################################

if skip_star == False:
    logging.info("\n STEP 1 \n Running STAR alignment!!!!!!! \n")

    cmd_align = "{}/src/1_RNA_star_basic.sh {} {} {} {} 4".format(vadir, RNA, ID, out_dir, CONFIG)
    logging.info(cmd_align)
    subprocess.call(cmd_align, shell = True)



##################################################################################################################################################################
##################################################################################################################################################################
# STEP 2
# refine_mapping
##################################################################################################################################################################
##################################################################################################################################################################

if skip_refine == False:
    logging.info("\n STEP 2: \n Refine Mapping \n")

    cmd_1 = "{}/src/2_RNA_refine_mapping.sh {} {} {} 4".format(vadir,ID,out_dir,CONFIG)
    logging.info(cmd_1)
    subprocess.call(cmd_1, shell = True)

##########################################################################################################################################################
##########################################################################################################################################################
################################################################# STEP 3.A ###############################################################################
################################################################## SNPiR #################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
if skip_snpir == False:
    #---------------------------------------------------------------------------------------------------------
    # STEP 3.1
    # Call raw variants:
    #       GATK UnifiedGenotyper
    #---------------------------------------------------------------------------------------------------------
    logging.info("\n STEP 3.A: SNPiR VARIANT CALLING")

    cd = "cd {}/tools/SNPiR".format(vadir)
    logging.info(cd)
    os.chdir("{}/tools/SNPiR".format(vadir))
    print(os.getcwd())

    # glm: genotype_likelihoods_model
    if skip_var == False:
        logging.info("\n GATK UnifiedGenotyper - Calling raw variants")
        cmd_2 = "java -jar /seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/VaDiR/tools/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper \
            -R {} \
            -I {} \
            -L {} \
            -glm SNP \
            -stand_call_conf 0 \
            -stand_emit_conf 0 \
            --dbsnp {} \
            -o {} \
            -nt 4".format(ref_genome_fasta, refined_bam, coding_bed, dbsnp, ug_out)
        logging.info(cmd_2)
        subprocess.call(cmd_2, shell = True)

        # check number of variants 
        phase_1 = open(ug_out,"r").readlines()
        print("\t SNPiR: phase_1: ", len(phase_1))

    #---------------------------------------------------------------------------------------------------------
    # # Filter raw variants according to SNPiR workflow
    #---------------------------------------------------------------------------------------------------------
    #---------------------------------
    # STEP 3.2
    # revised_convertVCF.sh
    #   Filter out variants that have a Variant Call Quality > 20
    #---------------------------------

    logging.info("\n STEP 2 \n SNPiR - Converting VCF to custom SNPiR format & filtering out variants with quality < 20")

    cmd_3 = "{}/tools/SNPiR/revised/revised_convertVCF.sh \
        {} \
        {}/snpir/1_filtered_phase_1 \
        20".format(vadir, ug_out, out_dir)
    logging.info(cmd_3)
    subprocess.call(cmd_3, shell = True)#.communicate()

    # check number of variants 
    file_phase_2 = "{}/snpir/1_filtered_phase_1".format(out_dir)
    phase_2 = open(file_phase_2, "r").readlines()
    print("\t SNPiR: phase_2: ", len(phase_2))

    #---------------------------------
    # STEP 3.3
    # filter_mismatch_first6bp_ln.pl
    #   Removes variations that occur within the first 6bp of teh read 
    #---------------------------------
    
    logging.info("\n STEP 3: \nFiltering out mismatches in first 6 bp of reads")

    cmd_4 = "perl {}/tools/SNPiR/revised/filter_mismatch_first6bp_ln.pl \
        -infile {}/snpir/1_filtered_phase_1 \
        -outfile {}/snpir/2_filtered_first6bp.txt \
        -bamfile {}".format(vadir,out_dir,out_dir,refined_bam)
    logging.info(cmd_4)
    subprocess.call(cmd_4, shell = True)#.communicate()

    # check number of variants 
    file_phase_3 = "{}/snpir/2_filtered_first6bp.txt".format(out_dir)
    phase_3 = open(file_phase_3, "r").readlines()
    print("\t SNPiR: phase_3: ", len(phase_3))

    #---------------------------------
    # STEP 3.4
    # bedtools subtract
    #   Removing Variations in repetitive regions using Repeat Mask Annotation 
    #---------------------------------
    # Using BEDtools subtract to remove sites in repetitive regions based on RepeatMasker annotation
    logging.info("\n STEP 4: \n Using BEDtools subtract to remove sites in repetitive regions based on RepeatMasker annotation")

    cmd_5 = "awk \'{{OFS=\"\\t\";$2=$2-1\"\\t\"$2;print $0}}\' {}/snpir/2_filtered_first6bp.txt \
        | {}/bedtools subtract \
            -a stdin \
            -b {}/tools/SNPiR/genome_ref/hg19_rmsk.bed \
        | cut -f1,3-7 \
        > {}/snpir/3_filtered_rmsk.txt".format(out_dir, bedtools, vadir, out_dir)
    logging.info(cmd_5)
    subprocess.call(cmd_5, shell = True)#.communicate()

    # check number of variants 
    file_phase_4 = "{}/snpir/3_filtered_rmsk.txt".format(out_dir)
    phase_4 = open(file_phase_4, "r").readlines()
    print("\t SNPiR: phase_4: ", len(phase_4))

    #---------------------------------
    # STEP 3.5
    # revised_filter_intron_near_splicejuncts.pl
    # Remove variants that are due to intronic alignment of RNA-seq reads 
    #---------------------------------
    # log "SNPiR - Filtering intronic candidates within 4 bp of splicing junctions"

    logging.info("\n STEP 5: \n Filtering intronic candidates within 4 bp of splicing junctions")

    cmd_6 = "perl {}/tools/SNPiR/revised/revised_filter_intron_near_splicejuncts.pl \
        -infile /{}/snpir/3_filtered_rmsk.txt \
        -outfile {}/snpir/4_filtered_sj.txt \
        -genefile {}/tools/SNPiR/revised/gene_annotation_table".format(vadir, out_dir, out_dir, vadir)
    logging.info(cmd_6)
    subprocess.call(cmd_6, shell = True)#.communicate()

    # check number of variants 
    file_phase_5 = "{}/tools/SNPiR/revised/gene_annotation_table".format(vadir)
    phase_5 = open(file_phase_5, "r").readlines()
    print("\t SNPiR: phase_5: ", len(phase_5))

    # # In ACF's installation of SNPiR, the gene annotation table is space-separated
    # # instead of tab-separated, so the above command uses a fixed one instead of:
    # #-genefile "$snpir_dir/genome_ref/gene_annotation_table"

    # #---------------------------------
    # # STEP 3.6
    # # filter_homopolymer_nucleotides.pl
    # #---------------------------------
    # log "SNPiR - Filtering candidates in homopolymer runs"

    logging.info("\n STEP 6: \n Filtering candidates in homopolymer runs")

    cmd_7 = "perl {}/tools/SNPiR/filter_homopolymer_nucleotides.pl \
        -infile {}/snpir/4_filtered_sj.txt \
        -outfile {}/snpir/5_filtered_homopolymer.txt \
        -refgenome {}".format(vadir,out_dir,out_dir,ref_genome_fasta)
    logging.info(cmd_7)
    subprocess.call(cmd_7, shell = True)#.communicate()

    # check number of variants 
    file_phase_6 = "{}/snpir/5_filtered_homopolymer.txt".format(out_dir)
    phase_6 = open(file_phase_6, "r").readlines()
    print("\t SNPiR: phase_6: ", len(phase_6))



    # ---------------------------------
    # STEP 3.7
    # pblat_candidates_ln.pl
    # Remapping variant reads using pblat 
    # ---------------------------------
    # log "SNPiR - Using PBLAT to ensure unique mapping"

    logging.info("\n STEP 7: \n Using PBLAT to ensure unique mapping")

    cmd_8 = "perl {}/tools/SNPiR/revised/pblat_candidates_ln.pl \
        -infile {}/snpir/5_filtered_homopolymer.txt \
        -outfile {}/snpir/6_filtered_pblat.txt \
        -bamfile {} \
        -refgenome {} \
        -pblatpath {} \
        -threads {} ".format(vadir,out_dir,out_dir,refined_bam, ref_genome_fasta, pblat, threads)
    logging.info(cmd_8)
    # os.system(cmd_8)
    subprocess.call(cmd_8, shell = True)#.communicate()


    # check number of variants 
    file_phase_7 = "{}/snpir/6_filtered_pblat.txt".format(out_dir)
    phase_7 = open(file_phase_7, "r").readlines()
    print("\t SNPiR: phase_7: ", len(phase_7))
    

    # #---------------------------------
    # # STEP 3.8
    # # Bedtools Subtract
    # #---------------------------------
    # log "SNPiR - Using BEDtools subtract to filter out known RNA editing sites"

    logging.info("\n STEP 8: \n Using BEDtools subtract to filter out known RNA editing sites")

    cmd_9 = "awk \'{{OFS=\"\\t\";$2=$2-1\"\\t\"$2;print $0}}\' {}/snpir/6_filtered_pblat.txt \
        | {}/bedtools subtract \
            -a stdin \
            -b {}/tools/SNPiR/genome_ref/Human_AG_all_hg19.bed \
        > {}/snpir/7_filtered_knownedits.bed".format(out_dir, bedtools, vadir, out_dir)
    logging.info(cmd_9)
    # os.system(cmd_9)
    subprocess.call(cmd_9, shell = True)#.communicate()


    # #---------------------------------
    # # STEP 3.9
    # # Bedtools Inspector 
    # #---------------------------------
    # log "Using SNPiR output & BEDtools intersect to extract final variants from raw VCF"

    logging.info("\n STEP 9: \n Using SNPiR output & BEDtools intersect to extract final variants from raw VCF")

    cmd_10 = "{}/bedtools intersect \
        -a {} \
        -b {}/snpir/7_filtered_knownedits.bed \
        -wa \
        -header \
        > {}/snpir/{}.ug.snpir.filtered.vcf".format(bedtools, ug_out, out_dir, out_dir, ID)
    logging.info(cmd_10)
    # os.system(cmd_10)
    subprocess.call(cmd_10, shell = True)#.communicate()


    # check number of variants 
    file_phase_9 = "{}/snpir/{}.ug.snpir.filtered.vcf".format(out_dir, ID)
    phase_9 = open(file_phase_9, "r").readlines()
    print("\t SNPiR: phase_9: ", len(phase_9))
    


##########################################################################################################################################################
##########################################################################################################################################################
################################################################## STEP 3.C ##############################################################################
################################################################### MuTect2 ##############################################################################
##########################################################################################################################################################
##########################################################################################################################################################

if skip_mutect == False:
    logging.info("\n STEP 3.C \n Running MuTect2. \n")
    # cmd_mutect = "{}/src/4c_filtered_calls_mutect.sh {} {} {} {}".format(vadir, ID, out_dir, CONFIG, threads)
    # logging.info(cmd_mutect)
    # subprocess.call(cmd_mutect, shell = True)#.communicate()

    cmd_source = "source {}".format(CONFIG)
    logging.info(cmd_source)
    subprocess.call(cmd_source, shell = True)#.communicate()

    log_ref = "Reference genome FASTA:  {}".format(ref_genome_fasta)
    log_dbsnp = "Known SNPs:  {}".format(dbsnp)
    log_bam = "Input RNA BAM:  {}".format(refined_bam)

    mutect_output = out_dir + "/mutect"
    log_output = "Output Directory:  {}".format(mutect_output)

    cmd_output = "mkdir -p {}".format(mutect_output)
    raw_vcf = mutect_output + "/{}.mt.raw.vcf".format(ID)

    logging.info(log_ref)
    logging.info(log_dbsnp)
    logging.info(log_bam)
    logging.info(log_output)
    subprocess.call(cmd_output, shell = True)#.communicate()


    #####################
    # Mutect2
    #####################
    if skip_var == False:
        logging.info("\nVariant calling using MuTect2\n")
        # cmd_mutect = "java -jar /seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/VaDiR/tools/GATK/GenomeAnalysisTK.jar -T MuTect2 \
        cmd_mutect = "java -jar /broad/software/free/Linux/redhat_6_x86_64/pkgs/GATK-3.6/GenomeAnalysisTK.jar -T MuTect2 \
            -R {} \
            -I:tumor {} \
            -dontUseSoftClippedBases \
            -stand_call_conf 20.0 \
            -stand_emit_conf 20.0 \
            --max_alt_alleles_in_normal_count 50 \
            --dbsnp {} \
            -o {} \
            -nct {}".format(ref_genome_fasta, refined_bam, dbsnp, raw_vcf, threads)
        logging.info(cmd_mutect)
        subprocess.call(cmd_mutect, shell = True)#.communicate()


    #################################################################################
    # Additionally, filter variants using the last 3 steps of the SNPiR workflow     
    #################################################################################

    logging.info("Preparing for SNPiR by removing variants that failed previous filter conditions")
    
    filtered_phase_1="{}/1_rmfailed_for_snpir.vcf".format(mutect_output)

    #"$PERL" -lane 'print $_ if $_ ~~ /^#/ or $F[6] eq "PASS"' "$raw_vcf"
    cmd_mutect_1 = "perl -lane \'print $_ if $_ ~~ /^#/ or $F[6] eq \"PASS\"\' {} > {}".format(raw_vcf, filtered_phase_1)
    logging.info(cmd_mutect_1)
    subprocess.call(cmd_mutect_1, shell = True)#.communicate()

    # change the directory to the sniper tool 
    snpir_dir = "{}/tools/SNPiR".format(vadir)
    logging.info(snpir_dir)
    os.chdir(snpir_dir)

    phase_1 = open(filtered_phase_1,"r").readlines()
    print("\t MUTECT2: phase_1: ", len(phase_1))

    # #------------------------------------------
    # # Converting VCF to custom SNPiR format
    # #------------------------------------------
    logging.info("SNPiR - Converting VCF to custom SNPiR format")
    
    filtered_phase_2 = "{}/2_snpir_format.txt".format(mutect_output)

    cmd_mutect_2 = "{}/revised/revised_convertVCF.sh \
                    {} \
                    {}".format(snpir_dir, filtered_phase_1, filtered_phase_2)
    logging.info(cmd_mutect_2)
    subprocess.call(cmd_mutect_2, shell = True)#.communicate()

    phase_2 = open(filtered_phase_2,"r").readlines()
    print("\t MUTECT2: phase_2: ", len(phase_2))

    # #------------------------------------------
    # # Filtering intronic candidates within 3 bp of splicing junctions
    # #------------------------------------------
    logging.info("SNPiR - Filtering intronic candidates within 3 bp of splicing junctions")

    filtered_phase_3 = "{}/3_filtered_sj.txt".format(mutect_output)

    cmd_mutect_3 = "perl {}/revised/revised_filter_intron_near_splicejuncts.pl \
        -infile {} \
        -outfile {} \
        -genefile {}/revised/gene_annotation_table".format(snpir_dir, filtered_phase_2, filtered_phase_3, snpir_dir)
    logging.info(cmd_mutect_3)
    subprocess.call(cmd_mutect_3, shell = True)#.communicate()

    # check number of variants 
    phase_3 = open(filtered_phase_3,"r").readlines()
    print("\t MUTECT2: phase 3: ", len(phase_3))
    



    # #------------------------------------------
    # # Phase 4
    # #------------------------------------------
    # logging.info("SNPiR - Using PBLAT to ensure unique mapping")
    filtered_phase_4="{}/4_filtered_pblat.txt".format(mutect_output)

    cmd_mutect_4 = "perl {}/revised/pblat_candidates_ln.pl \
        -infile {} \
        -outfile {} \
        -bamfile {} \
        -refgenome {} \
        -pblatpath {} \
        -threads {} ".format(snpir_dir,filtered_phase_3, filtered_phase_4, refined_bam, ref_genome_fasta, pblat, threads)
    logging.info(cmd_mutect_4)
    subprocess.call(cmd_mutect_4, shell = True)

    # check number of variants 
    phase_4 = open(filtered_phase_4,"r").readlines()
    print("\t MUTECT2: phase_4: ", len(phase_4))


    #------------------------------------------
    # bedtools subtract
    #------------------------------------------
    logging.info("SNPiR - Using BEDtools to separate out candidates that are known RNA editing sites")
    
    filtered_phase_5="{}/5_filtered_knownedits.bed".format(mutect_output)

    cmd_mutect_5 = "awk \'{{OFS=\"\\t\";$2=$2-1\"\\t\"$2;print $0}}\' {} | {}/bedtools subtract \
        -a stdin \
        -b {}/genome_ref/Human_AG_all_hg19.bed > {}".format(filtered_phase_4, bedtools, snpir_dir, filtered_phase_5)
    logging.info(cmd_mutect_5)
    subprocess.call(cmd_mutect_5, shell = True)

    # check number of variants 
    phase_5 = open(filtered_phase_5,"r").readlines()
    print("\t MUTECT2: phase_5: ", len(phase_5))


    #------------------------------------------
    # bedtools intersect
    #------------------------------------------
    logging.info("Using SNPiR output & BEDtools intersect to extract final variants from raw VCF")
    final_vcf="{}/sample_id.mt.snpir.filtered.vcf".format(mutect_output)

    cmd_mutect_6 = "{}/bedtools intersect \
        -a {} \
        -b {} \
        -wa \
        -header \
        > {}".format(bedtools, raw_vcf, filtered_phase_5, final_vcf)
    logging.info(cmd_mutect_6)
    subprocess.call(cmd_mutect_6, shell = True)#.communicate()

    # check number of variants 
    phase_6 = open(final_vcf,"r").readlines()
    print("\t MUTECT2: phase_6: ", len(phase_6))


    #------------------------------------------
    #
    #------------------------------------------
    final_mod_vcf="{}/{}.mt.snpir.filtered.mod.vcf".format(mutect_output, ID)

    cmd_mutect_7 = "grep \"#\" {} > {}".format(final_vcf, final_mod_vcf)

    logging.info(cmd_mutect_7)
    subprocess.call(cmd_mutect_7, shell = True)#.communicate()

    # check number of variants 
    phase_7 = open(final_mod_vcf,"r").readlines()
    print("\t MUTECT2: phase_7: ", len(phase_7))

    
    #------------------------------------------
    # 
    #------------------------------------------
    cmd_mutect_8 = "grep -v \"#\" {} | \
        awk -v OFS=\'\\t\' \'{{ \
            for(i=9;i<=11;++i){{ \
                        split($(i),a,\":\"); \
                            b[i]=a[1]\":\"a[2]}} \
                    print $1,$2,$3,$4,$5,$6,$7,$8,b[9],b[10],b[11]; \
                    }}\' \
        >> {}".format(final_vcf, final_mod_vcf)
    logging.info(cmd_mutect_8)
    subprocess.call(cmd_mutect_8, shell = True)#.communicate()


    # check number of variants 
    phase_8 = open(final_mod_vcf,"r").readlines()
    print("\t MUTECT2: phase_8: ", len(phase_8))


                                ##################################################################################
                #####################################################################################################################
##########################################################################################################################################################
################################################################### STEP VaDiR ###########################################################################
##########################################################################################################################################################
                #####################################################################################################################
                                ##################################################################################

if skip_VADiR == False:
    logging.info("\n Running VaDiR \n")
    #####################################################################################
    #####################################################################################
    # Step 5: MERGE VCFs
    #####################################################################################
    #####################################################################################
    logging.info("\nStep 5: MERGE VCFs\n")
    cmd_vadir_5 = "{}/src/5_merge_vcfs.sh {} {} {}".format(vadir, ID, out_dir, CONFIG)
    logging.info(cmd_vadir_5)
    subprocess.call(cmd_vadir_5, shell = True)


    # #####################################################################################
    # # Step 6: Samtools mpileup
    # #####################################################################################
    # cmd_vadir_6 = "src/6_mpileup.sh {} {} {}".format(ID, out_dir, CONFIG)
    logging.info("\nStep 6: Samtools mpileup\n")
    mkdir = "mkdir -p {}/sample_vcf/mpileup".format(out_dir)
    logging.info(mkdir)
    subprocess.call(mkdir, shell = True)

    input_bam="{}/refined-mapping/{}.refined.bam".format(out_dir, ID )
    inFile="{}/sample_vcf/{}.loc.txt".format(out_dir, ID )
    outDir="{}/sample_vcf/mpileup".format(out_dir, ID )
    # samtools="/seq/RNASEQ/mbrown/Anaconda/miniconda3/bin/samtools"
    samtools="/seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/VaDiR/tools/samtools-1.3/samtools"
    cmd_vadir_6 = "{} mpileup \
              -Buv -t DP,AD -q 40 \
              -f {} \
              -l {} \
              {} \
              -o {}/{}.q40.mpileup.vcf".format(samtools, ref_genome_fasta, inFile, input_bam, outDir, ID )
    
    logging.info(cmd_vadir_6)
    subprocess.call(cmd_vadir_6, shell = True)
    vadir_6 = "{}/{}.q40.mpileup.vcf".format(outDir, ID)
    step_6 = open(vadir_6,"r").readlines()
    step_6_vcf = [i for i in step_6 if i[0][0] != "#"]
    print("\t VaDiR: STEP 6: ", len(step_6_vcf), "\n")
    
    # #####################################################################################
    # # Step 7: ANNOTATING THE VARIANTS
    #      Annotate variants from the mpileup.vcf
    #           uses *annovar* to annotate the vcf file.
    #           Input:  mpileup.vcf
    #           output: annotated.txt
    # #####################################################################################
    logging.info("\nStep 7: ANNOTATING THE VARIANTS\n")
    cmd_vadir_7 = "{}/src/7_endo.annotate.sh {} {} {}".format(vadir, ID, out_dir, CONFIG)
    logging.info(cmd_vadir_7)
    subprocess.call(cmd_vadir_7, shell = True)

    vadir_7 = "{}/{}.q40.mpileup.sort.vcf.annotated.hg19_multianno.txt".format(out_dir, ID)
    step_7 = open(vadir_7,"r").readlines()
    step_7_vcf = [i for i in x if i[0][0] != "#"]
    print("\t VaDiR: STEP 7: ", len(step_7_vcf), "\n")

    # print("ANNOTATING THE VARIANTS")
    # # name1="$output_dir/sample_vcf/mpileup/$sample_id.q40.mpileup"
    # name1="{}/sample_vcf/mpileup/{}.q40.mpileup".format( out_dir, ID )
    # logging.info(name1)

    # #------------------
    # # 7.1: sort
    # # outout: ID.sort.vcf
    # #------------------
    # cmd_vadir_7_1 = "sort -k1,1d -k2,2n {}.vcf > {}.sort.vcf".format(name1, name1) #"$name1.vcf" > "$name1.sort.vcf"
    # logging.info(cmd_vadir_7_1)
    # subprocess.call(cmd_vadir_7_1, shell = True)

    # #------------------
    # # 7.2: annovar
    # # outout: ID.sort.avinput
    # #------------------
    # annovar_dir = "{}/tools/annovar".format(vadir)
    # # "$annovar_dir/convert2annovar.pl" \
    # #     --includeinfo \
    # #     --snpqual 0 \
    # #     -format vcf4old \
    # #     --outfile "$name1.sort.avinput" \
    # #     "$name1.sort.vcf" \
    # # # >> $output_dir/sample_vcf/mpileup/log.txt 2>&1
    # cmd_vadir_7_2 = "{}/convert2annovar.pl \
    #     --includeinfo \
    #     --snpqual 0 \
    #     -format vcf4old \
    #     --outfile {}.sort.avinput \
    #     {}.sort.vcf".format(annovar_dir, name1, name1)
    # logging.info(cmd_vadir_7_2)
    # subprocess.call(cmd_vadir_7_2, shell = True)

    # #------------------
    # # 7.3: 
    # # outout: ID.sort.mod.avinput
    # #------------------
    # # sed -e s/\,\<\\*\>//g "$name1.sort.avinput" \
    # # | sed -e s/\<\\*\>/-/g \
    # # > "$name1.sort.mod.avinput"
    # cmd_vadir_7_3 = "sed -e s/\\,\\<\\*\\>//g {}.sort.avinput | sed -e s/\\<\\*\\>/-/g > {}.sort.mod.avinput".format(name1, name1)
    # logging.info(cmd_vadir_7_3)
    # subprocess.call(cmd_vadir_7_3, shell = True)

    # #------------------
    # # 7.4: 
    # # outout: ID.sort.vcf.annotated
    # #------------------
    
    # # "$annovar_dir/table_annovar.pl" \
    # # "$name1.sort.mod.avinput" \
    # # "$annovar_dir/humandb/" \
    # # -buildver hg19 \
    # # -out $name1.sort.vcf.annotated \
    # # -remove \
    # # -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic70 \
    # # -operation g,r,r,f,f,f,f,f,f,f,f \
    # # --nastring . \
    # # --otherinfo \
    # # >> $output_dir/sample_vcf/mpileup/log.txt 2>&1
    
    # cmd_vadir_7_4 = "{}/table_annovar.pl \
    # {}.sort.mod.avinput \
    # {}/humandb/ \
    # -buildver hg19 \
    # -out {}.sort.vcf.annotated \
    # -remove \
    # -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic70 \
    # -operation g,r,r,f,f,f,f,f,f,f,f \
    # --nastring . \
    # --otherinfo ".format( annovar_dir, name1, annovar_dir, name1)

    # logging.info(cmd_vadir_7_4)
    # subprocess.call(cmd_vadir_7_4, shell = True)
    
    # #####################################################################################
    # #####################################################################################
    # # Step 8: 
    #           Final Hard Filtering Step
    # #####################################################################################
    # #####################################################################################
    logging.info("\n Step 8: Filter\n")
    cmd_vadir_8 = "$PERL {}/src/8_filter.pl {} {}".format(vadir, ID, out_dir)
    logging.info(cmd_vadir_8)
    subprocess.call(cmd_vadir_8, shell = True)

    vadir_8 = "{}/{}.annotated.filtered.txt".format(out_dir, ID)
    step_8 = open(vadir_8,"r").readlines()
    step_8_vcf = [i for i in x if i[0][0] != "#"]
    print("\t VaDiR: STEP 7: ", len(step_8_vcf), "\n")


    vcfFile = "{}/sample_vcf/{}.namesloc.txt.sort".format(out_dir, ID)
    tumorFile = "{}/sample_vcf/mpileup/{}.q40.mpileup.sort.vcf.annotated.hg19_multianno.txt".format(out_dir, ID)
    output = "{}/sample_vcf/{}.annotated.filtered.txt".format(out_dir, ID)


    #------------
    # VCF File
    #
    # Foramt:
    #   chr1    897325  G   C   Tier1,Tier2
    #   chr1    981931  A   G   Tier1,Tier2
    #------------
    # Open File
    vcfFile = open(vcfFile,"r").readlines()
    # VCF information
    vcf_info = []
    # Counter
    count2 = 0
    for i in vcfFile:
        if i[0] != "#":
            temp = i.strip().split("\t")
            tmpname = temp[0]+":"+temp[1]+";"+temp[2]+";"+temp[3]
            vcf_info.append([tmpname,temp[4]])
            count2 +=1
    print("done: read variant names: ",len(vcf_info))

    #-------------------
    # mpileup File
    #-------------------
    # Open File
    tumorFile = open(tumorFile,"r").readlines()
    # save the first line 
    firstline = tumorFile[0]
    tumorFile = tumorFile[1:]
    tumor_info = []
    count = 0
    for i in tumorFile:
        temp = i.strip().split("\t")
        tumor_info.append(temp)
        lp = len(temp)
        count += 1

    print("done: read tumor.mpileup: ",len(tumor_info))
    print("lp ",lp)


    #--------------------------------
    # compare variants with mpileup
    #--------------------------------

    numTumor = []
    num = 0
    for i in range(count2):
        name = tumor_info[i][lp-10] + ":" + tumor_info[i][lp-9] + ";" + tumor_info[i][lp-7] + ";" + tumor_info[i][lp-6]
        if(len(tumor_info[i][lp-6]) > 1):
            temp2 = tumor_info[i][lp-6].split(",")
            name = tumor_info[i][lp-10] + ":" + tumor_info[i][lp-9] + ";" + tumor_info[i][lp-7] + ";" + temp2[0]

        for j in range(num, count): 
            if vcf_info[i][0] == name:
                vcf_info[i].insert(lp+1, vcf_info[i][1])
                numTumor.append(j)
                num = j + 1

    print("done: compare variants with mpileup:\t", len(vcf_info), ":", len(numTumor))

    ##################################################################################################################
    # Hard Filtering
    # 1.) Germline read depth (DP) Greater than 5
    #           DP = filtered depth, number of filtered reads that support each of the reported alleles
    #           AD = unfiltered allele depth. number of reads that support each of the reported alleles
    # 2.) Germline Variant Fraction (VFA) less than 0.03
    ##################################################################################################################
    numTrue=[]
    ctTrue = 0
    pos = 0
    check = 0
    for p in range(0,len(numTumor)):
        # print(tumor_info[numTumor[p]])
        normInfo = tumor_info[numTumor[p]][-1].split(":")
        DPnorm = int(normInfo[1])
        # DP 
        if DPnorm >= 5:
            # CHECK: count the number of variants w/ DP >= 5
            check += 1
            normrefalt = normInfo[2].split(",") # AD
            # Append which variants pass this test 
            numTrue.append(p)
            
            # VFA
            frequency = int(normrefalt[1])/DPnorm
            if frequency <= 0.03:
                # numTrue.append(p)
                ctTrue+=1

    print(ctTrue)
    print(check)
    print("done: DP5 variants: ", len(numTrue))

    #------------
    # Write File 
    #------------

    output = open(output, "w")
    temp = firstline.rstrip()
    print(temp)
    output.write(temp)

    output.write("\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\tcombination\n")

    for p in numTrue:
        line = tumor_info[numTumor[p]]
        line = "\t".join(line)+"\n"
        # print(line)
        output.write(line)

    ##########################
    # Write VCF Formated file 
    ##########################
    test = "{}/sample_vcf/{}.annotated.filtered.txt_test.vcf".format(out_dir, ID)
    output2 = open(test, "w")
    output2.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    for p in numTrue:
        line = tumor_info[numTumor[p]][-10:]
        line = "\t".join(line)+"\n"
        output2.write(line)


# vcfFile.close()
# tumorFile.close()
# output.close()

