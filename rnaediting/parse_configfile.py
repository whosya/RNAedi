import configparser
import os
import logging
import re
import sys
from rnaediting import check_input
#import check_input
from .utils import can_i_run_software_, log_coloredlogs,init_logging,can_i_run_software
# from utils import can_i_run_software_, log_coloredlogs,init_logging,can_i_run_software 

class Configuration:
    init_logging("Read in  configration file", "DEBUG")
    def __init__(self, config_path):
        #self.config_path = os.path.join(os.getcwd(),filname)
        self.config_path = config_path
        if check_input.check_file(self.config_path).check_file_nonexistent_or_empty_():
            log_coloredlogs("Config file not exist, or empty, please generate it first.")
            sys.exit()
        else:
            self.config = configparser.ConfigParser()
            self.config.read(self.config_path)
    def get_single(self):
        single = self.config.get("allsamples_file", "single")
        if single == "":
            log_coloredlogs('Field "single" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        return single
    def get_sample_names(self):
        """
        Return string object containing all sample names with comma as its separator
        #1: check existence of samples in configfiles
        #2: filter illegal(with sapace) samplenames in "samples" field in configfiles
        """
        samples = self.config.get("allsamples_file", "samples")
        #1
        if samples == "":
            log_coloredlogs('Field "samples" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        #2
        samples_list = re.split("[\s,]+", samples.strip())
        if "" in samples_list:
            samples_list.remove("")
        wrong_name_l  = []
        for sam_name in samples_list:            
            if len(sam_name.split()) != 1:
                wrong_name_l.append(sam_name)
        if len(wrong_name_l) > 0:
            log_coloredlogs("Informal sample' names, must not contain any spaces or underscores. Please change the following names in the configuration file:%s" %"\n".join(wrong_name_l))
            sys.exit(1)
        #print(samples)
        return samples_list
        
    def get_separator(self):
        """
        Reture string
        """
        separator = self.config.get("allsamples_file", "separator")
        if separator == "":
            log_coloredlogs('Field "separator" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif separator == "null":
            log_coloredlogs('Field "separator" in configuration file is "null", indicating single-paired data, separator option should be droped;'
            'if not, separator symbol should be specified(e.g. "_")' ,err=True)
            sys.exit(1)
        separator = separator.strip()
        #print(separator)
        return separator

    def get_strand(self):
        """
        Reture list
        """
        strands = self.config.get("allsamples_file", "strand")
        if strands == "":
            log_coloredlogs('Field "strand" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif strands == "null":
            log_coloredlogs('Field "strand" in configuration file is "null", indicating single-paired data, strand option should be droped;'
            'if not, strand signal should be specified(e.g. "1,2")' ,err=True)
            sys.exit(1)
        else:
            strands_l = re.split("[\s,]+",strands.strip())
        #print(str(strands_l))
        return strands_l

    def get_sufix(self):
        sufix = self.config.get("allsamples_file", "sufix")
        if sufix == "":
            log_coloredlogs('Field "sufix" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        sufix = sufix.strip()
        #print(sufix)
        return sufix
    
class Get_software(Configuration):
    """
    child class
    """
    def __init__(self,config_path,name):
        super().__init__(config_path)
        self.name = name
    
    def get_fastp(self):
        fastp_path = self.config.get("software",self.name)
        if fastp_path == "":
            log_coloredlogs('Field "fastp" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif fastp_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(fastp_path):
                log_coloredlogs('The path you provided ([fastp_path]) may not exist, please check'.format(fastp_path),err=True)
                sys.exit(1)

        fastp_path = fastp_path.strip()
        #print(sufix)
        return fastp_path

    def get_hisat2(self):
        hisat2_path = self.config.get("software",self.name)
        if hisat2_path == "":
            log_coloredlogs('Field "hisat2" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif hisat2_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(hisat2_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(hisat2_path),err=True)
                sys.exit(1)
                
        hisat2_path = hisat2_path.strip()
        #print(sufix)
        return hisat2_path

    def get_hisat2_build(self):
        hisat2_build_path = self.config.get("software",self.name)
        if hisat2_build_path == "":
            log_coloredlogs('Field "hisat2_build" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif hisat2_build_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun [RNAedi generate-config]",'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(hisat2_build_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(hisat2_build_path),err=True)
                sys.exit(1)
                
        hisat2_build_path = hisat2_build_path.strip()
        #print(sufix)
        return hisat2_build_path

    def get_samtools(self):
        samtools_path = self.config.get("software",self.name)
        if samtools_path == "":
            log_coloredlogs('Field "samtools" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif samtools_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(samtools_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(samtools_path),err=True)
                sys.exit(1)
                
        samtools_path = samtools_path.strip()
        #print(sufix)
        return samtools_path
    
    def get_gatk3(self):
        gatk3_path = self.config.get("software",self.name)
        if gatk3_path == "":
            log_coloredlogs('Field "gatk3" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif gatk3_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(gatk3_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(gatk3_path),err=True)
                sys.exit(1)
                
        gatk3_path = gatk3_path.strip()
        #print(sufix)
        return gatk3_path

    def get_picard(self):
        picard_path = self.config.get("software",self.name)
        if picard_path == "":
            log_coloredlogs('Field "picard" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif picard_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(picard_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(picard_path),err=True)
                sys.exit(1)
                
        picard_path = picard_path.strip()
        #print(sufix)
        return picard_path

    def get_bcftools(self):
        bcftools_path = self.config.get("software",self.name)
        if bcftools_path == "":
            log_coloredlogs('Field "bcftools" in configuration file is empty, please fill in' ,err=True)
            sys.exit(1)
        elif bcftools_path == "null":
            log_coloredlogs('{} executable not found!,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(self.name),err=True)
            sys.exit(1)
        else:
            if not can_i_run_software_(bcftools_path):
                log_coloredlogs('The path you provided ([{}]) may not exist, please check'.format(bcftools_path),err=True)
                sys.exit(1)
                
        bcftools_path = bcftools_path.strip()
        #print(sufix)
        return bcftools_path
# for test
# if __name__ == '__main__' :
#     Configuration().get_sample_names()
#     Configuration().get_separator()
#     Configuration().get_strand()
#     Configuration().get_sufix()
        



