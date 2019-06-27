#!/usr/bin/env python

import argparse
import sys, os, re
import logging
import subprocess
import csv
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)




##########################
# Read the vcf inputs
##########################
def readVCF(infile, STR_VCF_DELIMITER, CHR_COMMENT):
    vcf_header = []
    vcf_body = []
    vcf_file = open(infile, "r").readlines()

    for line in vcf_file:
        # Keep comments
        ## Header 
        if line[0][0] == CHR_COMMENT:
            vcf_header.append(line)
            continue
        ## Variants
        else:
            vcf_body.append(line)
            
    return vcf_header, vcf_body

      ##############################################
############################################################
########################## Step 2 ##########################
############################################################
      ##############################################
def step_2(outdir, infile, quality_filter):
    
    logging.info("SNPiR - Step 2: Converting VCF to custom SNPiR format & filtering out variants with quality < 20")

    ##############
    # Constatnts #
    ##############
    VCF_DELIMITER = "\t"
    CHR_COMMENT = "#"

    # keep track of which lines are removed 
    removed=[]
    # output file
    outputFile = "{}/step2.txt".format(outdir)
    ################
    # Read in files 
    ################
    head, body = readVCF(infile, VCF_DELIMITER, CHR_COMMENT)

    # counters to see how many variants pass and fail 
    passed, failed = 0, 0

    outputFile = open(outputFile, "w")
    for i in range(len(body)):
        line = body[i].split("\t")

        # filter quality 
        if float(line[5]) >= quality_filter:
            # if GT and AD depths
            if re.search(r'[GT:AD]',line[8]):
                tmp = line[9].split(":")
                alleleDepth = tmp[1].split(",")

                referenceDepth = int(alleleDepth[0])
                altDepth = int(alleleDepth[1])

                if referenceDepth != 0 or altDepth != 0:
                    totalDepth = (referenceDepth + altDepth)

                    depths = str(totalDepth) + "," + str(altDepth)
                    altFraction = str(altDepth / totalDepth)
                    
                    newLine_list=[line[0], line[1], depths, line[3], line[4], altFraction]
                    newLine = "\t".join(newLine_list)

                    outputFile.write(newLine+"\n")

                    passed += 1

                else:
                    removed.append(i)
            else:
                removed.append(i)
        else:
            removed.append(i)
        failed = len(removed)
    outputFile.close()

    print("Variants kept:", passed)
    print("Variants filtered:", failed)

      ##############################################
############################################################
########################## Step 3 ##########################
############################################################
      ##############################################
def step_3(outdir, vadir, bamFile):
    
    logging.info("SNPiR - STEP 3: Filtering out mismatches in first 6 bp of reads")

    ##############
    # Constants #
    ##############
    quality_offset = 33
    minimal_base_quality = 25
    minmismatch = 1
    
    step3_file = "{}/step2.txt".format(outdir)
    step3_outfile = "{}/step3.txt".format(outdir)
    outfile_failed = "{}/step3_failed.txt".format(outdir)

    TEMP = step3_outfile + '_tmp'

    infile = open(step3_file, "r")
    outfile = open(step3_outfile , "w")
    outfile_failed = open (outfile_failed, "w")


    # counters 
    not_filtered, removed= 0, 0

    for i in infile.readlines():
        line = i.split("\t")
        chrom, position = line[0], line[1]
        bamposition = chrom + ':' + position+'-'+position

        cmd = "samtools view {} {} > {}".format(bamFile, bamposition, TEMP)
        subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()

        editnuc = line[4]
        newmismatch = 0
        mismatchreadcount = 0
        newcov, newmismatch = 0, 0 
        f = open(TEMP,"r")

        for j in f.readlines():
            bamfields = j.strip().split("\t")
            alignment, readstart, cigar, sequence, qualities = bamfields[1], bamfields[3], bamfields[5], bamfields[9], bamfields[10]
            sequencebases = sequence#.split("//")
            qualscores = qualities#.split(":")[-1]

            currentpos, readpos = int(readstart), 1
            base_readpos = []
            
            # leters 
            cigarletters = re.findall(r'[MIDNSHP=X]',cigar)
            # numbers
            cigarnums = re.split('[MIDNSHP]', cigar)


            for k in range(len(cigarletters)):
                position = int(position)
                if currentpos > position:
                    break
                if cigarletters[k] == "I" or cigarletters[k] == "S":
                # if cigarletters[i] =~ m/[SI]/) {
                    readpos = readpos + int(cigarnums[k])

                elif cigarletters[k] == "D" or cigarletters[k] == "N": 
                    currentpos = currentpos + int(cigarnums[k])
                
                elif cigarletters[k] == "M":
                    for j in range(int(cigarnums[k])):
                        if currentpos == position:
                            base_readpos = readpos
                        currentpos += 1
                        readpos += 1
            
            if base_readpos:
                # set the reverse indicator to 0
                revers_strand = 0 
                # revers_strand = 1 if (alignment & 16)
                # check the SAM lag, (1000 means this is a revers strand)
                a = bin(int(alignment))
                rev_test = a[-5]
                if len(a) >= 5 and int(rev_test) == 1:
                    revers_strand = 1
                if (revers_strand == 0 and int(base_readpos) > 6) or (revers_strand == 1 and base_readpos < readpos - 5):
                    if ord(str(qualscores)[base_readpos-1]) >= minimal_base_quality + quality_offset:
                        newcov += 1
                        if (sequencebases[base_readpos-1] == editnuc):
                            newmismatch+=1
                        else:
                            basequalFail=1
                    else:
                        readPosFail=1

        # remove the temp file 
        # cmd = "rm {}".format(TEMP)
        # subprocess.Popen(cmd, shell=True)

        #-----------------------
        # output lines 
        #-----------------------
        if newmismatch >= minmismatch: 
            not_filtered += 1
            varfreq = (newmismatch/newcov)
            new_output_line = (line[0] + "\t" + line[1]+ "\t" + str(newcov) +","+ str(newmismatch) + "\t" + line[3] + "\t" + line[4] + "\t" + str(varfreq)+"\n")
            outfile.write(new_output_line)
        #-----------------------
        # output removed/filtered lines 
        #-----------------------
        if newmismatch < minmismatch:
            removed += 1
            if basequalFail:
                reason = ' low basequal;'
            if readPosFail:
                reason = ' mismatch at readend;'
            
            new_output_line = (line[0] + "\t" + line[1]+ "\t" + str(newcov) +","+ str(newmismatch) + "\t" + line[3] + "\t" + line[4] + "\t" + "NAN\treason:"+reason + "\n")
            outfile_failed.write(new_output_line)

    # report how many removed and remained 
    info_message = "Variants removed: {} \n Varirants Remaining: {}".format(removed, not_filtered)
    logging.info(info_message)

    infile.close()
    outfile.close()
    outfile_failed.close()

      ##############################################
############################################################
########################## Step 4 ##########################
############################################################
      ##############################################

def step_4(outdir, vadir):
    
    logging.info("\n STEP 4: \n Using BEDtools subtract to remove sites in repetitive regions based on RepeatMasker annotation")

    ##############
    # Constatnts #
    ##############
    ## Store the vcf in list 
    new_vcf = []
    
    # String to write lines to 
    temp=""

    ################
    # Read in files 
    ################
    step4_output_path = "{}/step4.txt".format(outdir)
    step4_infile = "{}/step3.txt".format(outdir)
    # repeat masker bed file 
    bed_path = '{}/tools/SNPiR/genome_ref/hg19_rmsk.bed'.format(vadir)

    vcf_file = open(step4_infile, "r")
    bed_file = open(bed_path, "r")
    
    # loop over the vcf file variants 
    for i in vcf_file.readlines():
        line = i.split("\t")
        start_location = int(line[1])-1
        line.insert(1, str(start_location))
        new_vcf.append("\t".join(line))
    temp = "".join(new_vcf)
    # print(temp)
    
    # temp is the stdin
    # Use bedtools subtract to remove variations found in rmsk.bed file 
    cmd = "bedtools subtract -a stdin -b {}/tools/SNPiR/genome_ref/hg19_rmsk.bed | cut -f1,3-7 > {}".format(vadir, step4_output_path)
    print(cmd)
    subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf8', shell=True).communicate(temp)


    # Close the files 
    vcf_file.close()
    bed_file.close()

      ##############################################
############################################################
########################## Step 5 ##########################
############################################################
      ##############################################

def intron_splice_region(variant_pos, variant_genes, splice_dist):
    
    ###################################################################
    # For each variant, 
    #    1) loop through its chromosome's gene annotations
    #    2) & check if variant is in intronic region near a splice junction
    ###################################################################
    for i in variant_genes:
        gene_start = int(i[4])
        gene_end = int(i[5])
        if variant_pos < gene_start:
            break
        
        # falls in the region 
            ## if if in the gene
        if gene_start <= variant_pos and variant_pos <= gene_end:
            exon_count = int(i[8])
            exon_starts = i[9].split(",")
            exon_ends = i[10].split(",")

            # For each exon, check if falls within 4 bp of the start or end 
            for j in range(exon_count):

                # Does the variant fall in the intronic reagon 
                exon_start = int(exon_starts[j])
                start_intronic_offset = exon_start - splice_dist
                ## if the position of the variant is not in the intronic reagon, then break 
                if variant_pos < start_intronic_offset:
                    break

                exon_end = int(exon_ends[j])
                end_intronic_offset = exon_end + splice_dist

                # Is the variant within the first splice_distance bp's of the START of the exon 
                start_intron_region = start_intronic_offset < variant_pos and variant_pos <= exon_start
                # Is the variant within the first splice_distance bp's of the END of the exon
                end_intron_region = exon_end < variant_pos and variant_pos <= end_intronic_offset

                # If the variant is found to fall within the splice_distance bp's of the start or end of the exon region 
                if start_intron_region or end_intron_region:
                    return(True)

    return(False)

def step_5(outdir, vadir):

    logging.info("\n STEP 5: \n Filtering intronic candidates within 4 bp of splicing junctions")

    ###############
    # Constants 
    ###############

    # set file paths 
    step5_outfile = "{}/step5.txt".format(outdir)
    step5_infile = "{}/step4.txt".format(outdir)

    genefile_path = "{}/tools/SNPiR/revised/gene_annotation_table".format(vadir)

    # Open Files Read and Write
    infile = open(step5_infile, 'r')
    genefile = open(genefile_path, 'r')
    outfile = open(step5_outfile, 'w')
    outfile_failed = open(step5_outfile + "_failed.txt", 'w')

    splice_dist = 4
    infile_list=[]
    variants = []

    # counters to see how many variants pass and fail 
    passed, failed = 0, 0


    ##############################################
    # create dictionary that holds the chromosome 
    # and the line from gene annotations table 
    ##############################################
    dic = {}

    # loop over gene file
    for i in genefile.readlines():
        x = i.rstrip('\n').split("\t")
        gene_chrom = x[2]
        if gene_chrom in dic:
            dic[gene_chrom].append(x)
        else: 
            dic.setdefault(gene_chrom, []).append(x)

    # loop over infile 
    variant_chrom_genes_ref=[]
    for i in infile.readlines():
        x = i.rstrip('\n').split("\t")
        chrom = x[0]
        pos = int(x[1])
        variant_chrom_genes_ref = dic[chrom]

        if intron_splice_region(pos, variant_chrom_genes_ref, splice_dist) == True:
            outfile_failed.write(i)
            failed += 1
        else:
            outfile.write(i)
            passed += 1

    outfile.close()
    outfile_failed.close()
    infile.close()
    genefile.close()

    print("Variants kept:", passed)
    print("Variants filtered:", failed)


      ##############################################
############################################################
########################## Step 6 ##########################
############################################################
      ##############################################

def step_6(outdir, vadir):

    logging.info("\n STEP 6: \n Filtering candidates in homopolymer runs")

    ###############
    # Constants 
    ###############
    # Set File Paths 
    infile_path = "{}/step5.txt".format(outdir)
    outfile_path = "{}/step6.txt".format(outdir)
    refgenome_path = "{}/refs/ucsc.hg19.fasta".format(vadir)

    # distance from splice region
    splice_dist = 4
    infile_list = []

    # BED file constants 
    left_buffer, right_buffer = 4, 4

    # counters to see how many variants pass and fail 
    passed, failed = 0,0
    #----------------------
    # Open Files 
    #----------------------
    infile = open(infile_path, 'r')
    refgenome = open(refgenome_path, 'r')
    outfile = open(outfile_path, 'w')
    outfile_failed = open(outfile_path + "_failed.txt", 'w')

    temp = ""
    temp_bed_path = "{}/tmp.bed".format(outdir)
    temp_bed = open(temp_bed_path, "w")

    fastaFromBed = "{}/tools/bedtools-2.25.0/fastaFromBed".format(vadir)
    cmd = "{} -fi {} -bed stdin -fo stdout".format(fastaFromBed, refgenome_path)

    #-----------------------------------------------------------------------------  
    # 1) Read in the input file (output from step 5) and run through each line 
    # 2) Create the Bed file (Prep for Bedtools fastaFromBed)
    # 3) Create the FASTA file by running Run FastaFromBed 
    # 4) loop over the output file sequences to see if homopolymers
    #-----------------------------------------------------------------------------
    #------------------------
    # 1) Create the bed file
    #------------------------
    for i in infile.readlines():
        #------------------------
        # 2) Create the bed file
        #------------------------
        x = i.strip().split("\t")
        chrom, position, edit_base_nuc = x[0], int(x[1]), x[4]
        # set the constant homopolymer to false for this line 
        homopolymer = False
        # get start and end positions 
        start_position, end_position = position - left_buffer, position + right_buffer + 1
        # create the BED file line 
        new_line = chrom + "\t" + str(start_position) + "\t" + str(end_position) + '\n'

        #------------------------------------
        # 3) Run the command for fastaFromBed
        #------------------------------------
        fasta_output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf8', shell=True).communicate(input=new_line)[0]
        fasta = fasta_output.strip().split("\n")
        sequence = ""
        #----------------------------------------
        # 4) loop over the output file sequences 
        #----------------------------------------
        edit_base_nuc = edit_base_nuc.upper()
        sequence += fasta[1].upper()

            # iterate of the sequence 
            #   check if the 4 consecutive nucleotides are the edited nucleotide 
            #   if true will set th variable homopolymer to TRUE
        for k in range(len(sequence)-4):
            if edit_base_nuc == sequence[k] and edit_base_nuc == sequence[k+1] and edit_base_nuc == sequence[k+2] and edit_base_nuc == sequence[k+3]:
                homopolymer = True
            else:
                homopolymer = False

        #-----------------
        # Write to a file 
        #-----------------
        if homopolymer == False:
            passed+=1
            outfile.write(i)
        else:
            failed+=1
            outfile_failed.write(i)

    # remove temp files 
    cmd = "rm {}.psl*".format(outfile_path)
    subprocess.Popen(cmd, shell=True)

    # close files 
    outfile.close()
    outfile_failed.close()
    infile.close()
    refgenome.close()

    print("Variants kept: ", passed)
    print("Variants filtered: ", failed)


      ##############################################
############################################################
########################## Step 7 ##########################
############################################################
      ##############################################


def step_7(outdir, vadir, threads, bamFile):

    logging.info("\n SNPiR - STEP 7: Running BLAT for realignment to remap reads containing variants.")


    ##############
    # Constants #
    ##############
    
    Phredscore_endcoding_offset = 33
    minimum_base_quality = 25
    minimum_mismatch = 1
    score_limit = 0.95
    
    # File Paths input and output 
    step7_infile_path = "{}/step6.txt".format(outdir)
    step7_outfile_path = "{}/step7.txt".format(outdir)
    outfile_failed_path = "{}/step7_failed.txt".format(outdir)
    refgenome_path = "{}/refs/ucsc.hg19.fasta".format(vadir)

    fa_file_path = "{}.fa".format(step7_outfile_path)
    psl_file_path = "{}.psl".format(step7_outfile_path)
    TEMP = step7_outfile_path + '_tmp'

    pblat_path = '{}/tools/pblat'.format(vadir)

    step_7_infile = open(step7_infile_path, "r")
    fa_file = open(fa_file_path, "w")


    # counters to see how many variants pass and fail 
    passed, failed = 0,0

    ########
    # 1) Prep file for BLAT
    #       use samtools view to get the alignments and get reads containing variants 
    # 2) 
    ########

    ####################################
    # Prep file for BLAT
    ####################################
    # get reads containing variants 
    # Loop over each variant 
    infile = step_7_infile.readlines()
    for i in infile:
        #-----------------------------------------
        # Get the reads that contain the variants 
        # samtools view: prints alignments that are in the given input alignment file
        #-----------------------------------------
        line = i.split("\t")
        chrom, position = line[0], line[1]
        bamposition = chrom + ':' + position+'-'+position

        cmd = "samtools view {} {} > {}".format(bamFile, bamposition, TEMP)
        subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()

        editnuc = line[4]
        newmismatch = 0
        mismatchreadcount = 0
        newcov, newmismatch = 0, 0 

        f = open(TEMP,"r")

        #----------------------------------------
        # read through the input file 
        #----------------------------------------
        for j in f.readlines():
            bam_fields = j.strip().split("\t")
            alignment, readstart, cigar, sequence_bases, quality_scores = bam_fields[1], bam_fields[3], bam_fields[5], bam_fields[9], bam_fields[10]

            current_pos, readpos = int(readstart), 1
            base_readpos = []
            
            # leters 
            cigarletters = re.findall(r'[MIDNSHP=X]',cigar)
            # numbers
            cigarnums = re.split('[MIDNSHP]', cigar)


            for k in range(len(cigarletters)):

                #### it is now faster ####
                position = int(position)
                if current_pos > position:
                    break
                ##### error corrected ####
                if cigarletters[k] == "S" or cigarletters[k] == "I":
                # if cigarletters[i] =~ m/[SI]/) {
                    readpos = readpos + int(cigarnums[k])

                elif cigarletters[k] == "D" or cigarletters[k] == "N": 
                    current_pos = current_pos + int(cigarnums[k])
                
                elif cigarletters[k] == "M":
                    for j in range(int(cigarnums[k])):
                        if current_pos == position:
                            if sequence_bases[readpos-1] == editnuc:
                                if ord(quality_scores[readpos-1]) >= minimum_base_quality + Phredscore_endcoding_offset:
                                    base_readpos = 1
                        current_pos += 1
                        readpos += 1
            if base_readpos:
                # increment miss matched read count by one 
                mismatchreadcount+=1
                # write the line to the FA file 
                fa_line = ">" + chrom + "-" + str(position) + "-" + str(mismatchreadcount) + "\n" + sequence_bases + "\n"
                # print(fa_line)
                fa_file.write(fa_line)
    fa_file.close()
    step_7_infile.close()

    #----------------------------------------
    # RUN PBLAT: 
    #       use pblat (parallel BLAT) for realignment of the reads that contain the variants 
    #----------------------------------------

    message = "Runnning PBLAT to remap variant reads. Threads: {}".format(threads)
    logging.info(message)

    # Path to blat/pblat
    pblat = "{}/pblat".format(pblat_path)
    cmd_pblat = "{} \
                 -threads={} \
                 -stepSize=5 \
                 -repMatch=2253 \
                 -minScore=20 \
                 -minIdentity=0 \
                 -noHead {} {} {}".format(pblat, threads, refgenome_path, fa_file_path, psl_file_path)
    print(cmd_pblat)
    subprocess.Popen(cmd_pblat, shell=True).communicate()

    #----------------------------------------
    # Process the PBLAT output
    #----------------------------------------
    message = "Process the PBLAT output"
    logging.info(message)

    # separate the file 
    psl_dict = {}
    psl = open(psl_file_path, "r")
    for i in psl.readlines():
        line = i.strip().split("\t")
        name = line[9]

        line_list = [line[0], line[13],line[17],line[18],line[20]]
        # add to dictionary 
        psl_dict.setdefault(name, []).append(line_list)
    psl.close()

    keep_dict = {}
    filter_dict = {}

    # Key format example 
    #   'chr1-14542-1'

    for i in psl_dict:
        tmp = i.split("-")
        current_chromosome = tmp[0]
        current_position = int(tmp[1])
        psl_id = tmp[0] + "_" + tmp[1]

        # set the max line to the first line 
        max_scored_line = psl_dict[i][0]
        max_score = max_scored_line[0]
        list_scores = []
        # loop over the contents in the dictionary
        for j in psl_dict[i]:
            if int(j[0]) > int(max_score):
                max_scored_line = j
                max_score = int(j[0])
            list_scores.append(int(j[0]))

        #-------------------------------------------
        # Check to see if the reads are overlapping 
        # make sure the second best hit (based on score) has a score thats is < 95% of the best hit 
        #-------------------------------------------
        overlaping = 0
        # if more than one read was found
        if len(list_scores) == 1:
            list_scores.append(0)
        list_scores.sort(reverse = True) 
        # check chromosomes to make sure match 
        if max_scored_line[1] == current_chromosome:
            
            # check second best hit score 
            if list_scores[1] < (list_scores[0]*score_limit):

                block_count, block_sizes, block_starts = max_scored_line[2], max_scored_line[3].split(","), max_scored_line[4].split(",")
                for k in range(int(block_count)):
                    start_position = int(block_starts[k])+1
                    end_position = int(block_starts[k]) + int(block_sizes[k])
                    # check if positions are overlapping 
                    if (current_position >= start_position) and (current_position <= end_position):
                        overlaping = 1
                if overlaping:
                    if (psl_id in keep_dict): 
                        keep_dict[psl_id] += 1
                    else:
                        keep_dict[psl_id] = 1
        if overlaping == 0:
            if (psl_id in filter_dict): 
                filter_dict[psl_id] += 1
            else:
                filter_dict[psl_id] = 1

    
    step7_outfile = open(step7_outfile_path, "w")
    step7_failed = open(outfile_failed_path, "w")

    for i in infile:
        line = i.strip().split("\t")
        tmp = line[2].split(",")
        coverage, old_alter = int(tmp[0]), int(tmp[1])
        psl_id = line[0] + "_" + line[1]

        if psl_id in keep_dict:
            new_alter = int(keep_dict[psl_id])

            if psl_id in filter_dict:
                discard_counter = filter_dict[psl_id]
            else:
                discard_counter = 0

            new_coverage = coverage - (old_alter - new_alter)
            new_edit_freq = new_alter/new_coverage

            if (new_alter >= minimum_mismatch) and (new_alter > discard_counter):
                new_output_line =line[0] + "\t" + line[1] + "\t" + str(new_coverage) + "," + str(new_alter) + "\t" + line[3] + "\t" + line[4] + "\t" + str(new_edit_freq) + "\n"
                step7_outfile.write(new_output_line)
                passed += 1
            else:
                failed_output_line = i + "\tfailed_freq #mismatches: " + str(new_alter) + "\tminimumMismatchesNecessary: " + str(minimum_mismatch) + "\tdiscarded mismatch reads: " + str(discard_counter) + "\n"
                step7_failed.write(failed_output_line)
                failed += 1
        else:
            failed_output_line = i + "\ttotalcover_failed\n"
            step7_failed.write(failed_output_line)
            failed += 1

    step7_outfile.close()
    step7_failed.close()

    print("Variants kept:", passed)
    print("Variants filtered:", failed)



def step_8(outdir, vadir):
    #---------------------------------
    # STEP 3.8
    # Bedtools Subtract
    #---------------------------------
    logging.info("\n STEP 8: \n Filtering out currently known RNA-editing sites using BEDtools subtract.")

    #--------------------
    # Comstants 
    #--------------------
    infile = "{}/step7.txt".format(outdir)
    TMP = "{}/step8_tmp.txt".format(outdir)
    outfile = "{}/step8.bed".format(outdir)
    str_infile =""
    #--------------------
    # read in the file 
    #--------------------
    inputFile = open(infile, "r")

    # edit the bed file to have start and stop positions, also round the last value to third decimal 
    # pass the file as a variable to stdin 
    if infile:
        for line in inputFile.readlines(): 
            line = line.strip().split("\t")
            start = int(line[1])-1
            line.insert(1, str(start))
            line[-1] = str(format(round(float(line[-1]),3),'.3f'))
            str_infile+="\t".join(line)+"\n"

    #--------------------
    # Run command 
    #--------------------
    cmd = "bedtools subtract \
        -a stdin \
        -b {}/tools/SNPiR/genome_ref/Human_AG_all_hg19.bed > {}".format(vadir, outfile)
    print(cmd)
    subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate(str_infile)



def step_9(outdir, vadir, infile):
    #---------------------------------
    # STEP 3.9
    # Bedtools intersect
    #---------------------------------
    bedtools = "{}/tools/bedtools-2.25.0".format(vadir)
    logging.info("\n STEP 9: \n Bedtools intersect")

    cmd = "bedtools intersect \
        -a {} \
        -b {}/step8.bed \
        -wa \
        -header \
        > {}/step9.snpir.filtered.vcf".format(infile, outdir, outdir)
    print(cmd)
    subprocess.Popen(cmd, shell=True).communicate()

##############################################################################
##############################################################################
# MAIN: run the steps 
##############################################################################
##############################################################################
def main():
    infile = "/seq/RNASEQ/mbrown/test.vcf"
    outdir = "/seq/RNASEQ/mbrown"
    # infile = "/seq/RNASEQ/mbrown/CTAT/Mutations/SNPiR/variants.vcf"
    # outdir = "/seq/RNASEQ/mbrown/CTAT/Mutations/SNPiR"
    vadir = "/seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/VaDiR"
    quality_filter = 20
    ID = "SRR5665260"
    # refined_bam =  "{}/refined-mapping/{}.refined.bam".format(outdir, ID)
    refined_bam =  "/seq/RNASEQ/mbrown/CTAT/VaDiR/VaDiR/GIAB_DATA/refined-mapping/SRR5665260.refined.bam"

    # step_2(outdir, infile, quality_filter)
    # step_3(outdir, vadir, refined_bam)
    # step_4(outdir, vadir)
    # step_5(outdir, vadir)
    # step_6(outdir, vadir)
    # step_7(outdir, vadir, threads=8, bamFile = refined_bam)
    step_8(outdir, vadir)
    step_9(outdir, vadir, infile)

main()