#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys, re

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)



import datetime
import subprocess

import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger('ctat_mutations')
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"]))
from Pipeliner import Pipeliner, Command, run_cmd, ParallelCommandList




MAYO_SCRIPTS_DIR = "/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/src"
GENEBED = "/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/resources/exons.bed"
GENOME = "/seq/RNASEQ/__ctat_genome_lib_building/August_15_2019/GRCh37_gencode_v19_CTAT_lib_Aug152019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"

SAMTOOLS="/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/bin/samtools/"
BLAT="/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/bin/blatSrc/linux.x86_64/"
BEDTOOLS="/seq/RNASEQ/TOOLS/RVBOOST/RVboost_0.1/bin/bedtools/bin/"

def main(argv):

    usage = "usage: {} target.vcf output.vcf\n\n".format(argv[0])
    if len(argv) < 3:
        print(usage, file=sys.stderr)
        sys.exit(1)

    input_vcf_filename = sys.argv[1]
    output_vcf_filename = sys.argv[2]

    tmpdir = "__" + os.path.basename(output_vcf_filename) + ".tmp_dir"

    chckpts_dir = tmpdir
    pipeliner = Pipeliner(chckpts_dir)
    
    

    
    
    ## DJ annot
    vcf_w_DJ_annot = "{}/wDJ.vcf".format(tmpdir)
    cmd = "perl {}/get_exon_cvg.pl {} {} {} > {}".format(MAYO_SCRIPTS_DIR,
                                                         GENEBED,
                                                         input_vcf_filename,
                                                         BEDTOOLS,
                                                         vcf_w_DJ_annot)

    pipeliner.add_commands([Command(cmd, "dj.ok")])
    ## ED annot
    
    
    #$PERL $MAYO_scripts/vcf_blat_verify.pl -i $output/$sample.raw.vcf -r $REF_GENOME -o $output/$sample.ED.vcf -w 50 -b $BLAT -sam $SAMTOOLS -br $REF_GENOME
    cmd = "perl {}/vcf_blat_verify.pl -i {} -r {} -o {} -w 50 -b {} -sam {} -br {}".format(MAYO_SCRIPTS_DIR,
                                                                                           vcf_w_DJ_annot,
                                                                                        GENOME,
                                                                                           output_vcf_filename,
                                                                                           BLAT,
                                                                                           SAMTOOLS,
                                                                                           GENOME)
    
    pipeliner.add_commands([Command(cmd, "blat_ed.ok")])

    pipeliner.run()

    sys.exit(0)
    


if __name__ == '__main__':
    main(sys.argv)
    
