from distutils.cmd import Command
import subprocess
import os
import sys
import logging
from .utils import log_subprocess,init_logging,log_coloredlogs
import datetime
from .check_input import check_file
from .parse_configfile import Get_software,Configuration


def genome_index(genome, config="conf.txt", hisat2_build="hisat2-build", samtools = "samtools", threads = 4):
    """
    This method indexing genome into different format
    """
    print("%s : Start" %datetime.datetime.now())
    init_logging("Indexing genome", "DEBUG")

    cwd = os.getcwd()
    genome_f = check_file(genome).check_file_nonexistent_or_empty()
    myoutdir = os.path.join(cwd,"1_genome_index")
    try:
        log_coloredlogs("Creating genome index folder in %s......" %myoutdir)
        os.mkdir(myoutdir)    
    except FileExistsError as e:
        log_coloredlogs("The output path of '1_genome_index' exist, the old file will be overwrited")

    #hisat2_build_ab = Get_software(config,hisat2_build).get_hisat2_build()
    samtools_ab = Get_software(config,samtools).get_samtools()
    #hisat2-build
    # log_coloredlogs("Building genome index with hisat2-build. this may take some time please wait......")
    # hisatx = subprocess.run([hisat2_build_ab, "-p", str(threads), genome_f, os.path.join(myoutdir,genome)],stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    # log_subprocess("hisat2-build",hisatx)

    #samtools faidx
    # log_coloredlogs("Building genome index with 'samtools faidx'.")
    # samtoolsx = subprocess.run([samtools_ab, "faidx", genome_f, "-o", os.path.join(myoutdir,genome)], stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    # log_subprocess("samtools faidx",samtoolsx)

    #picard CreateSequenceDictionary
    # log_coloredlogs("Building genome index with 'picard CreateSequenceDictionary'.")
    # picardx = subprocess.run(["picard", "CreateSequenceDictionary", genome_f], stdout= subprocess.PIPE,stderr= subprocess.PIPE)
    # log_subprocess("picard CreateSequenceDictionary",picardx)

    print("%s : All done !" %datetime.datetime.now())

    return myoutdir

class Runbcftools:

    def __init__(self,conf="conf.txt"):

        self.conf = conf
        self.sample_l = Configuration(conf).get_sample_names()
        self.sep = Configuration(conf).get_separator()
        self.strand = Configuration(conf).get_strand()
        self.sufix = Configuration(conf).get_sufix()
        self.single = Configuration(conf).get_single()

    def run_fastp(self,*kwargs,raw_file,outdir="0_cleaned_data",threads = 4, fastp = "fastp"):  ############## ""
        """
        filter raw reads data using fastp with its defult parameters
        """
        print("%s : Start" %datetime.datetime.now())
        init_logging("Running fastp", "DEBUG")
        #check file and software and output from last run
        ind = check_file(raw_file).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun(outdir)
        fastp_ab = Get_software(self.conf,fastp).get_fastp()

        #
        command = []
        for samples_name in self.sample_l:
            if self.single.lower() == "false":
                command = ["fastp", "-i", ind+os.sep+samples_name+self.sep+self.strand[0]+self.sufix, "-I", ind+os.sep+samples_name+self.sep+self.strand[1]+self.sufix,
                "-o", outdir+os.sep+samples_name+self.sep+self.strand[0]+"_cleaned"+self.sufix, "-O", outdir+os.sep+samples_name+self.sep+self.strand[1]+"_cleaned"+self.sufix,
                "-j", outdir+os.sep+samples_name+".json", "-h", outdir+os.sep+samples_name+".html", "-q", str(threads)]
                print(" ".join(command))
                if  kwargs:
                    # for k,v in kwargs.items():
                    #     command.append(k)
                    #     command.append(v)
                    for k in kwargs:
                        command.append(k)

            else:
                command = [fastp_ab, "-i", ind+os.sep+samples_name+self.sufix, 
                "-o", outdir+os.sep+samples_name+"_cleaned"+self.sufix, 
                "-j", outdir+os.sep+samples_name+".json", "-h", outdir+os.sep+samples_name+".html", "-q", str(threads)]
                if  kwargs:
                    # for k,v in kwargs.items():
                        # command.append(k)
                        # command.append(v)
                    for k in kwargs:
                        command.append(k)
                        
            print(" ".join(command))
            fastpx = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("filter with fastp",fastpx)

        return outdir
    
    def run_hisat2(self, *args, genome, fastp=True, hisat2="hisat2", samtools="samtools", threads=4, ind = "0_cleaned_data", outdir="2_aligned_data", rg="WES",lb="WES",pl="ILLUMINA"):
        
        print("%s : Start" %datetime.datetime.now())
        init_logging("Running hisat2", "DEBUG")
        #check file and software and output from last run
        genome_ab = check_file(genome).check_file_nonexistent_or_empty()
        genome = os.path.basename(genome_ab)
        ind = check_file(ind).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun(outdir)
        hisat2_ab = Get_software(self.conf,hisat2).get_hisat2()
        samtools_ab = Get_software(self.conf,samtools).get_samtools()

        #
        command=None
        for samples_name in self.sample_l:
            if fastp:
                if self.single.lower() == "false":
                    command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", "LB:"+lb, "--rg", "PL:"+pl,
                    "-p", str(threads), "-x", "./1_genome_index"+os.sep+genome, "-"+"1", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[0]+"_cleaned"+self.sufix,
                    "-"+"2", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[1]+"_cleaned"+self.sufix,"|", samtools_ab, "sort","-@", str(threads), "-", "-o",outdir+os.sep+samples_name+".sorted.bam"]

                    if args:
                        for k in args:
                            command.insert(1,k)
                        print(" ".join(command))
                else:
                    command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", "LB:"+lb, "--rg", "PL:"+pl,
                    "-p", str(threads), "-x", "./1_genome_index"+os.sep+genome, "-U", "./0_cleaned_data"+os.sep+samples_name+"_cleaned"+self.sufix,
                    "|", samtools_ab, "sort","-@", str(threads), "-", "-o",
                    outdir+os.sep+samples_name+".sorted.bam"]
                    if args:
                        for k in args:
                            command.insert(1,k)
                        print(" ".join(command))
            else:
                if self.single.lower() == "false":
                    command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", "LB:"+lb, "--rg", "PL:"+pl,
                    "-p", str(threads), "-x", "./1_genome_index"+os.sep+genome, "-1", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[0]+self.sufix,
                    "-2", "./0_cleaned_data"+os.sep+samples_name+self.sep+self.strand[1]+self.sufix, "|", samtools_ab, "sort","-@", str(threads), "-", "-o",
                    outdir+os.sep+samples_name+".sorted.bam"]
                    if not args:
                        for k in args:
                            command.insert(1,k)
                        print(" ".join(command))
                else:
                    command = [hisat2_ab, "--rg-id", samples_name+rg, "--rg", "SM:"+samples_name, "--rg", "LB:"+lb, "--rg", "PL:"+pl,
                    "-p", str(threads), "-x", "./1_genome_index"+os.sep+genome, "-U", "./0_cleaned_data"+os.sep+samples_name+self.sufix,
                    "|", samtools_ab, "sort","-@", str(threads), "-", "-o",
                    outdir+os.sep+samples_name+".sorted.bam"]
                    if not args:
                        for k in args:
                            command.insert(1,k)
                        print(" ".join(command))
            #print(" ".join(command))
##############################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!, in certain cercumstance, use shell=True may be a better way?
            hisat2x = subprocess.run(" ".join(command), check=True, shell=True, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("Aligning with hisat2",hisat2x)

        return outdir

    def run_indel_realign(self, genome, gatk3="GenomeAnalysisTK",samtools= "samtools", ind="2_aligned_data", tempdir = "realign_temp_intervals", outdir ="2_2_realigned_data"):
        """
        This is optional method, if true it need GATK3 and java
        """
        print("%s : Start" %datetime.datetime.now())
        init_logging("Running Realigner with gatk3 [GenomeAnalysisTK]", "DEBUG")
        #check file and software and output from last run
        genome = check_file(genome).check_file_nonexistent_or_empty()
        ind = check_file(ind).check_dir_nonexistent_or_empty()
        tempdir = check_file(tempdir).check_file_from_lastrun(tempdir)
        outdir = check_file(outdir).check_file_from_lastrun(outdir)
        gatk3_ab = Get_software(self.conf,gatk3).get_gatk3()
        samtools_ab =  Get_software(self.conf,samtools).get_samtools()

        for samples_name in self.sample_l:
            command1 = ["java", "-Xmx4g", "-jar", gatk3_ab, "-T", "RealignerTargetCreator", "-R", genome, "-I", ind+os.sep+samples_name+".sorted.bam",
            "-o", tempdir+os.sep+samples_name+".sorted.intervals.bam"]

            command2 = ["java", "-Xmx4g", "-jar", gatk3_ab, "-T", "IndelRealigner", "-R", genome, "-I", ind+os.sep+samples_name+".sorted.bam",
            "-targetIntervals", tempdir+os.sep+samples_name+".sorted.intervals.bam", "-o", outdir+os.sep+samples_name+".sorted.realigned.bam"]

            command3 = [samtools_ab, "index", outdir+os.sep+samples_name+".sorted.realigned.bam"]

            realignx1 = subprocess.run(command1, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("Realigning with gatk3 [RealignerTargetCreator]",realignx1)

            realignx2 = subprocess.run(command2, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("Realigning with gatk3 [IndelRealigner]",realignx2)

            realignx3 = subprocess.run(command3, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            log_subprocess("Indexing with samtools",realignx3)
        
        return outdir    

    def run_markdup(self, realign = False, picard = "picard", ind = "2_aligned_data", outdir="3_marked_data"):
        
        '''
        This method will marking duplications for alignd files:
        In default, we use not realigned data as input, this will save space and time in terms of SNP calling, the outdir is named as '3_marked_data'.
        not that, For INDEL calling, if you have sufficient time and space, the realigning step will be needed, which will improve the INDEL result, and this option have been set as
        default step, if not you can drop this step by set '--realign False ', the outdir is named as '3_2_marked_realigned_data'
        '''

        print("%s : Start" %datetime.datetime.now())
        init_logging("Running picard MarkDuplicates", "DEBUG")
        #check file and software and output from last run
        ind = check_file(ind).check_dir_nonexistent_or_empty()
        outdir = check_file(outdir).check_file_from_lastrun(outdir)
        picard_ab = Get_software(self.conf,picard).get_picard()
        
        for samples_name in self.sample_l:
            if not realign:
                command = [picard_ab, "MarkDuplicates", "-I", ind+os.sep+samples_name+".sorted.bam", "-O", 
                outdir+os.sep+samples_name+".sorted.markdup.bam", "-M", outdir+os.sep+samples_name+".sorted.markdup.txt", "--CREATE_INDEX", "true"]

                markdupx = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
                log_subprocess("MarkDuplicates with picard",markdupx)
            else:
                command = [picard_ab, "MarkDuplicates", "-I", ind+os.sep+samples_name+".sorted.realigned.bam", "-O", 
                outdir+os.sep+samples_name+".sorted.realigned.markdup.bam", "-M", outdir+os.sep+samples_name+".sorted.realigned.markdup.txt", "--CREATE_INDEX", "true"]

                markdupx = subprocess.run(command, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
                log_subprocess("MarkDuplicates with picard",markdupx)
        
        return outdir

    

class  Callvar(Runbcftools):

    def __init__(self, genome, markdup=True, conf="conf.txt", bcftools="bcftools", ind="3_marked_data", temp_outdir="temp_mpilup", threads=4):
        super().__init__(conf)
        self.genome = check_file(genome).check_file_nonexistent_or_empty()
        self.markdup = markdup
        self.bcftools_ab = Get_software(self.conf,bcftools).get_bcftools()
        self.ind = check_file(ind).check_dir_nonexistent_or_empty()
        self.temp_outdir = check_file(temp_outdir).check_file_from_lastrun(temp_outdir)
        self.threads = threads
        
        
        init_logging("Running mpilup with bcftools", "DEBUG")
        print("%s : Start" %datetime.datetime.now())

        # this method produce different results depending on if  'markdup' was used, with different data as input:
        # markdup=True  ------> 3_marked_data
        # markdup=False  ------>  2_aligned_data
        
        for samples_name in self.sample_l:
            if self.markdup:
                command_snp_mpilup = [self.bcftools_ab, "mpileup", "--threads", str(self.threads), "-d", "10000", "-C", "50", "-E", "-Q", "25", "-q", "20", "-O", "b", "-f", self.genome, 
                self.ind+os.sep+samples_name+".sorted.markdup.bam", "-o", self.temp_outdir+os.sep+samples_name+".bcf"]

                snp_mpilupx = subprocess.run(command_snp_mpilup, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
                log_subprocess("Running mpilup",snp_mpilupx)
            else:
                command_snp_mpilup = [self.bcftools_ab, "mpileup", "--threads", str(self.threads), "-d", "10000", "-C", "50", "-E", "-Q", "25", "-q", "20", "-O", "b", "-f", self.genome, 
                self.ind+os.sep+samples_name+".sorted.bam", "-o", self.temp_outdir+os.sep+samples_name+".bcf"]

                snp_mpilupx = subprocess.run(command_snp_mpilup, stdout= subprocess.PIPE,stderr= subprocess.PIPE)
                log_subprocess("Running mpilup",snp_mpilupx)

    def run_bcftools_snp(self, outdir="4_callsnp_result"):

        # genome = 
        # ind = 
        # temp_outdir = 
        outdir = check_file(outdir).check_file_from_lastrun(outdir)
        # bcftools_ab = 

        init_logging("Running snp calling with bcftools", "DEBUG")
        print("%s : Start" %datetime.datetime.now())
        for samples_name in self.sample_l:

            command_snp_call = [self.bcftools_ab, "call", "--threads", str(self.threads), "-vm", "-V", "indels", "-o", outdir+os.sep+samples_name+".snp.vcf", self.temp_outdir+os.sep+samples_name+".bcf"]
            print(" ".join(command_snp_call))
            snp_callx = subprocess.run(command_snp_call, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
            log_subprocess("Snp calling",snp_callx)
            return outdir

    def run_bcftools_indel(self, realign = True, ind="4_2_realigned_data", outdir="5_call_indel_result"):

        '''
        Rsult will be different regards to True or False of [realign] and [markdup], 
        ind="3_2_marked_realigned_data"  <------>    [realign = True]+[markdup=True]     <------> 5_call_indel_result
        ind="2_2_realigned_data"         <------>    [realign = True]+[markdup=False]    <------> 5_call_indel_result
        temp_outdir="temp_mpilup"        <------>    [realign = False]+[markdup=True]    <------> 5_call_indel_result
        temp_outdir="temp_mpilup"        <------>    [realign = False]+[markdup=False]   <------> 5_call_indel_result
        '''

        init_logging("Running indel calling with bcftools", "DEBUG")
        outdir = check_file(outdir).check_file_from_lastrun(outdir)

        if realign:
            if self.markdup:
                ind = check_file(ind).check_dir_nonexistent_or_empty()

                for samples_name in self.sample_l:
                    command_indel_mpilup_call = [self.bcftools_ab, "mpileup", "--threads", str(self.threads), "-d", "10000", "-C", "50", "-E", "-Q", "25", "-q", "20", "-O", "b", "-f", self.genome, 
                ind+os.sep+samples_name+".sorted.realigned.markdup.bam", "|", self.bcftools_ab, "call", "-vmO", "--threads", str(self.threads), "-V", "snps", "-o", outdir+os.sep+samples_name+".indel.vcf"]

                    indel_mpilupx = subprocess.run(command_indel_mpilup_call, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
                    log_subprocess("Running mpilup",indel_mpilupx)
            else:

                ind = check_file(ind).check_dir_nonexistent_or_empty()

                for samples_name in self.sample_l:
                    command_indel_mpilup_call = [self.bcftools_ab, "mpileup", "--threads", str(self.threads), "-d", "10000", "-C", "50", "-E", "-Q", "25", "-q", "20", "-O", "b", "-f", self.genome, 
                ind+os.sep+samples_name+".sorted.realigned.bam", "|", self.bcftools_ab, "call", "-vmO", "--threads", str(self.threads), "-V", "snps", "-o", outdir+os.sep+samples_name+".indel.vcf"]

                    indel_mpilupx = subprocess.run(command_indel_mpilup_call, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
                    log_subprocess("Running mpilup",indel_mpilupx)
        else:
            for samples_name in self.sample_l:    
                command_indel_call = [self.bcftools_ab, "call", "-vmO", "--threads", str(self.threads), "-V", "snps", "-o", outdir+os.sep+samples_name+".indel.vcf", self.temp_outdir+os.sep+samples_name+".bcf"]
                indel_callx = subprocess.run(command_indel_call, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
                log_subprocess("Running indel calling",indel_callx)

# if __name__ == '__main__':
#     argop = ["-z","2","-A"]
#     Runbcftools(conf="/home/whosy/workstation/RNAedi_testdata/conf.txt").run_fastp(*argop,raw_file="/home/whosy/workstation/RNAedi_testdata/testdata",outdir="0_cleaned_data",threads = 4, fastp = "fastp")