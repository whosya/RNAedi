import configparser
#from itertools import count
import logging
from .check_input import check_file
import os
import re
from .utils import log_coloredlogs,init_logging,can_i_run_software
import sys

#try decorator
init_logging("Will generate configration file", "DEBUG")
#logging.basicConfig(format='[%(levelname)s] [%(asctime)s] [%(filename)s:%(lineno)d] %(message)s', level="DEBUG", datefmt='%H:%M:%S',stream=sys.stdout)
def generate_configfile(rawreads_path, configfile_name="conf.txt",  single = False, custom_file=None, strand = "1,2", separator="_", sufix=None,
software = "fastp,hisat2,hisat2-build,bcftools,gatk,picard,samtools,bwa,GenomeAnalysisTK"): ####################### None
    """
    Get all filenames in rawreads_path, generating configfile from "rawreads_path" provided by users, or
    optionally, 

    """
    configfile_check = check_file(configfile_name)
    if not configfile_check.check_file_nonexistent_or_empty_():
        log_coloredlogs("Configuration file [%s] exist already, also Please check all parameters carefully, and fix them if neccessory." %configfile_name)
        sys.exit(1)

    rawreads_l = os.listdir(rawreads_path)
    rawreads_file = check_file(rawreads_path)
    rawreads_file.check_dir_nonexistent_or_empty()

    if custom_file:
        """
        Need a user-created sample file :with each sample name in each line.
        in this case, three param needed:
        :param custom_file: user-created sample file
        :param separator: (e.g. _ )
        :param sufix: a sufix(e.g. .fq.gz)
        """
        
        custom_file_p = os.path.abspath(custom_file)
        custom_file_check = check_file(custom_file_p)
        custom_file_check.check_file_nonexistent_or_empty()
        
        with open(custom_file_p) as f:
            sample_names_new = ",".join([lines.strip() for lines in f])
            sufix_conf = sufix
        if single:
            strand = "null"
            separator = "null"

    else:
        """
        This method will generate configfile according to the input rawreads files
        #1: check "strand"
        #2: check "sufix"

        """
        sample_names = ",".join(rawreads_l)
        if not single:
            #1
            if int(len(rawreads_l)/2) == int(sample_names.count("_1")):
                # print(sample_names)
                # print(int(len(rawreads_l)/2))
                # print(int(sample_names.count("_1")))
                pass
            elif len(rawreads_l)/2 == sample_names.count("_L"):
                strand = "L,R"
            else:
                log_coloredlogs("your file names were not coincident with our standard format, please check your filename, optionally, you can use --sample_file --separator --sufix\n",err=True)
                sys.exit()
            #2
            posi_sufix = [".fq.gz",".fq"]
            sample_names_new = ""                                                                      ############################### avoid unbound
            sufix_conf = ""
            flag = 0
            for sufixs in posi_sufix:
                if not sample_names.count(sufixs) == 0:
                    
                    if len(rawreads_l) == sample_names.count(sufixs):
                        flag =1
                        sample_names_new = re.sub("_[12RL]"+sufixs,"",sample_names).split(",")
                        sample_names_new = ",".join(list(set(sample_names_new)))
                        log_coloredlogs("All of your sample names were extracted as follows:\n\n%s\n" %sample_names_new)
                        log_coloredlogs("If any of sample names were wrong, please use --sample_file --separator --sufix\n",warn=True)
                        sufix_conf = sufixs
                        break
            if flag == 0:
                
                log_coloredlogs("your file names were not coincident with our standard format, please use --sample_file --separator --sufix\n",err=True)
                sys.exit()
        else:
            #1
            strand = "null"
            separator = "null"
            
            #2
            posi_sufix = [".fq.gz",".fq"]
            sample_names_new = ""                                                                      ############################### avoid unbound
            sufix_conf = ""
            flag = 0
            for sufixs in posi_sufix:
                if not sample_names.count(sufixs) == 0:
                    
                    if len(rawreads_l) == sample_names.count(sufixs):
                        flag =1
                        sample_names_new = sample_names.split(",")
                        sample_names_new = ",".join(list(set(sample_names_new)))
                        log_coloredlogs("All of your sample names were extracted as follows:\n\n%s\n" %sample_names_new)
                        log_coloredlogs("If any of sample names were wrong, please use --sample_file --separator --sufix\n",warn=True)
                        sufix_conf = sufixs
                        break
            if flag == 0:
                
                log_coloredlogs("your file names were not coincident with our standard format, please use --sample_file --separator --sufix\n",err=True)
                sys.exit()
    
    cfgfile = open(configfile_name, "w+")
    Config = configparser.ConfigParser(allow_no_value=True)
    Config.add_section("allsamples_file")

    if not single:
        Config.set("allsamples_file", "single", "false")
        Config.set("allsamples_file", "# Your file type: true/false\n")
    else:
        Config.set("allsamples_file", "single", "true")
        Config.set("allsamples_file", "# Your file type\n")

    Config.set("allsamples_file", "samples", sample_names_new)
    Config.set("allsamples_file", "# Providing all your sample names (with comma as separator)\n")

    Config.set("allsamples_file", "separator", separator)
    Config.set("allsamples_file", "# providing a separator (e.g. '_' in [sample_1.fq]) \n")

    Config.set("allsamples_file", "strand", strand)
    Config.set("allsamples_file", "# providing a pair-end strand(e.g. '1' in [sample_1.fq])\n")

    Config.set("allsamples_file", "sufix", sufix_conf)
    Config.set("allsamples_file", "# providing a common sufix name(e.g. '.fq' in [sample_1.fq])\n")

    Config.add_section("software")
    mysoftwares = can_i_run_software(software)
    print(mysoftwares[1])
    for name,path in mysoftwares[1].items():
        Config.set("software", name, path)
        Config.set("software", "# providing the absolute path of %s\n" %name)


    log_coloredlogs("Configuration file [%s] not found: it will be now generated" %configfile_name)
    log_coloredlogs("Please check all parameters carefully, and fix them if neccessory.")
    Config.write(cfgfile)
    cfgfile.close()

#FOR TEST
if __name__ == '__main__':
    #configfile_name = "cof.txt"
    rawreads_path = "/home/whosy/workstation/rnaediting/test1"
    # generate_configfile(rawreads_path ,custom_file="/home/whosy/workstation/rnaediting/test2/sample.txt")
    # generate_configfile(rawreads_path)
    # generate_configfile(rawreads_path ,custom_file="/home/whosy/workstation/rnaediting/test2/sample.txt",single=True)
    generate_configfile(rawreads_path)

