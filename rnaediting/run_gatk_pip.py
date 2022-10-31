from distutils.cmd import Command
import subprocess
import os
import sys
import logging
from .utils import log_subprocess,init_logging,log_coloredlogs
import datetime
from check_input import check_file
from parse_configfile import Get_software,Configuration


def genome_index(genome, threads = 4):
    """
    This method index genome in different format
    """
    print("%s : Start" %datetime.date)
    init_logging("Indexing genome", "DEBUG")

    cwd = os.getcwd()
    genome_f = check_file(genome).check_file_nonexistent_or_empty()
    myoutdir = os.path.join(cwd,"1_genome_index")
    try:
        log_coloredlogs("Creating genome index folder in %s......" %myoutdir)
        os.mkdir(myoutdir)    
    except FileExistsError as e:
        log_coloredlogs("The output path of '1_genome_index' existing, the old file will be rewrited")

    #hisat2-build
    log_coloredlogs("Building genome index with hisat2-build.")
    hisatx = subprocess.run(["hisat2-build", "–p", str(threads), genome_f],stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    log_subprocess("hisat2-build",hisatx)

    #samtools faidx
    log_coloredlogs("Building genome index with 'samtools faidx'.")
    samtoolsx = subprocess.run(["samtools", "faidx", genome_f], stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    log_subprocess("samtools faidx",samtoolsx)

    #picard CreateSequenceDictionary
    log_coloredlogs("Building genome index with 'picard CreateSequenceDictionary'.")
    picardx = subprocess.run(["picard", "CreateSequenceDictionary", genome_f], stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    log_subprocess("picard CreateSequenceDictionary",picardx)

    print("%s : All done !" %datetime.date)

    return myoutdir

class Rungatk:

    def __init__(self,conf="conf.txt"):

        self.conf = conf
        self.sample_l = Configuration(conf).get_sample_names()
        self.sep = Configuration(conf).get_separator()
        self.strand = Configuration(conf).get_strand()
        self.sufix = Configuration(conf).get_sufix()
        self.single = Configuration(conf).get_single()

    def run_fastp(self,raw_file,outdir="0_cleaned_data",threads = 4, fastp = "fastp",**kwargs):  ############## ""
        """
        filter raw reads data using fastp with its defult parameters
        """
        print("%s : Start" %datetime.date)
        init_logging("Running fastp", "DEBUG")
        #check file and software and output from last run
        ind = check_file(raw_file).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun()
        fastp_ab = Get_software(self.conf,fastp) 

        #
        command = []
        for samples_name in self.sample_l:
            if self.single.lower() == "false":
                command = ["fastp", "-i", ind+os.sep+samples_name+self.sep+self.strand[0]+self.sufix, "-I", ind+os.sep+samples_name+self.sep+self.strand[1]+self.sufix,
                "-o", outdir+os.sep+samples_name+self.sep+self.strand[0]+"_cleaned"+self.sufix, "-O", outdir+os.sep+samples_name+self.sep+self.strand[1]+"_cleaned"+self.sufix,
                "-j", outdir+os.sep+samples_name+".json", "-h", outdir+os.sep+samples_name+".html", "-q", threads]
                if not kwargs:
                    for k,v in kwargs.items():
                        command.append(k)
                        command.append(v)
            else:
                command = [fastp_ab, "-i", ind+os.sep+samples_name+self.sufix, 
                "-o", outdir+os.sep+samples_name+"_cleaned"+self.sufix, 
                "-j", outdir+os.sep+samples_name+".json", "-h", outdir+os.sep+samples_name+".html", "-q", threads]
                if not kwargs:
                    for k,v in kwargs.items():
                        command.append(k)
                        command.append(v)
            fastpx = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("filter with fastp",fastpx)

        return outdir
    
    def run_hisat2(self, genome, single = False, hisat2="hisat2", samtools="samtools", threads=4, ind = "0_cleaned_data", outdir="2_aligned_data", rg="WES",lb="LB:WES",pl="PL:ILLUMINA"):
        
        print("%s : Start" %datetime.date)
        init_logging("Running hisat2", "DEBUG")
        #check file and software and output from last run
        genome = check_file(genome).check_file_nonexistent_or_empty()
        ind = check_file(ind).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun()
        hisat2_ab = Get_software(self.conf,hisat2) 
        samtools_ab = Get_software(self.conf,samtools)

        #
        for samples_name in self.sample_l:
            if self.single.lower() == "false":
                command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", lb, "--rg", pl,
                "-p", threads, "-x", "./1_genome_index"+os.sep+genome, "-1", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[0]+"_cleaned"+self.sufix,
                "-2", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[1]+"_cleaned"+self.sufix, "|", samtools_ab, "sort","-@", threads, "-", "-o",
                outdir+os.sep+samples_name+".sorted.bam"]
            else:
                command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", lb, "--rg", pl,
                "-p", threads, "-x", "./1_genome_index"+os.sep+genome, "-U", "./0_cleaned_data"+os.sep+samples_name+"_cleaned"+self.sufix,
                "|", samtools_ab, "sort","-@", threads, "-", "-o",
                outdir+os.sep+samples_name+".sorted.bam"]
            hisat2x = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("filter with fastp",hisat2x)

        return outdir

    def run_markdup(self,picard = "picard", ind = "2_aligned_data", outdir="3_marked_data"):

        print("%s : Start" %datetime.date)
        init_logging("Running picard MarkDuplicates", "DEBUG")
        #check file and software and output from last run
        ind = check_file(ind).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun()
        picard_ab = Get_software(self.conf,picard) 
        
        for samples_name in self.sample_l:
            command = ["picard_ab", "MarkDuplicates", "-I", ind+os.sep+samples_name+".sorted.bam", "-O", 
            outdir+os.sep+samples_name+".sorted.markdup.bam", "-M", outdir+os.sep+samples_name+".sorted.markdup.txt", "--CREATE_INDEX", "true"]

            markdupx = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("MarkDuplicates with picard",markdupx)
        
        return outdir

class  Callvar(Rungatk):

    def __init__(self, genome, conf="conf.txt", gatk="gatk", ind="3_marked_data", threads=4):
        super().__init__(conf)

        self.genome = check_file(genome).check_file_nonexistent_or_empty()
        self.gatk_ab = Get_software(self.conf,gatk)
        self.ind = check_file(ind).check_dir_nonexistent_or_empty()
        self.threads = threads

    def gatk_call_snp(self, memory="4g",outdir1="4_snp_call",outdir2="5_snp_select",outdir3="6_snp_filter",outdir4="7_snp_filtered"):

        outdir1 = check_file(outdir1).check_file_from_lastrun()
        outdir2 = check_file(outdir2).check_file_from_lastrun()
        outdir3 = check_file(outdir3).check_file_from_lastrun()
        outdir4 = check_file(outdir4).check_file_from_lastrun()


        init_logging("Running snp calling with GATK4", "DEBUG")
        for samples_name in self.sample_l:
            command1 = [self.gatk_ab, "--java-options", "-Xmx"+memory, "HaplotypeCaller", "-R", self.genome, "-I", self.ind+os.sep+samples_name+".sorted.markdup.bam", "-O", 
            outdir1+os.sep+samples_name+".call.vcf", "--native-pair-hmm-threads", self.threads]
            snpcallx1 = subprocess.run(command1, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("HaplotypeCaller",snpcallx1)

            command2 = [self.gatk_ab, "--java-options", "-Xmx"+memory, "SelectVariants", "-R", self.genome, "-V", outdir1+os.sep+samples_name+".call.vcf", "--select-type", "SNP", "-O", 
            outdir2+os.sep+samples_name+".selected.vcf"]
            snpcallx2 = subprocess.run(command2, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("SelectVariants",snpcallx2)

            command3 = [self.gatk_ab, "--java-options", "-Xmx"+memory, "VariantFiltration", "-R", self.genome, "-V", outdir2+os.sep+samples_name+".selected.vcf", "--filter-expression", 
            "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0", "--filter-name", "SNP_filter", "-O", outdir3+os.sep+samples_name+".filter.snp.vcf"]
            snpcallx3 = subprocess.run(command3, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("SelectVariants",snpcallx3)

            command4 = [self.gatk_ab, "--java-options", "-Xmx"+memory, "SelectVariants", "-R", self.genome, "-V", outdir3+os.sep+samples_name+".filter.snp.vcf", "--exclude-filtered", "-O", 
            outdir4+os.sep+samples_name+".filtered.snp.vcf"]
            snpcallx4 = subprocess.run(command4, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("SelectVariants",snpcallx4)

        return outdir4










    