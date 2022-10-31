import logging
import coloredlogs
import math
import datetime
import subprocess
import os
import sys

def init_logging(log_heading, logging_level):
    """
    Initialize the logging environment and print a header to the log.

    :param log_heading: heading for the log
    :return: nothing
    """
    logging.basicConfig(format='[%(levelname)s] [%(asctime)s] [%(filename)s:%(lineno)d] %(message)s', level=logging_level, stream=sys.stdout)
    length = math.ceil((len(log_heading)/2))
    log_coloredlogs('# ' * length)
    log_coloredlogs(log_heading)
    log_coloredlogs(datetime.datetime.today().ctime())
    log_coloredlogs('# ' * length)

def log_subprocess(program, process):
    """
    Log output from a subprocess call to debug log stream

    :param program: program name
    :param process: completed subprocess object
    """
    logging.debug('{} stdout:\n'.format(program) +
                  process.stdout.decode('utf-8'))
    # logging.debug('{} stderr:\n'.format(program) + process.stderr.decode('utf-8'))
    log_coloredlogs('{} stderr:\n'.format(program) + process.stderr.decode('utf-8'),err=True)

def log_coloredlogs(log_content,info =True,debug=False,warn=False,err=False):
    FIELD_STYLES = dict(
    asctime=dict(color='green'),
    hostname=dict(color='magenta'),
    levelname=dict(color='green'),
    filename=dict(color='magenta'),
    name=dict(color='blue'),
    threadName=dict(color='green')
    )

    LEVEL_STYLES = dict(
        debug=dict(color='green'),
        info=dict(color='cyan'),
        warning=dict(color='yellow'),
        error=dict(color='red'),
        critical=dict(color='red')
    )

    logger = logging.getLogger('tos')
    coloredlogs.install(
        level="DEBUG",
        fmt="[%(levelname)s] [%(asctime)s] [%(filename)s:%(lineno)d] %(message)s",
        level_styles=LEVEL_STYLES,
        field_styles=FIELD_STYLES)
    if debug:
        logger.debug(log_content)
    elif warn:
        logger.warn(log_content)
    elif err:
        logger.error(log_content) 
    else:
        logger.info(log_content)                      #only coloring erro and warn info
    #logger.critical('This is critical mode')


def uniq_id():
    """
    Get a unique ID which is not crazy long

    :return: string
    """
    from time import time
    return str(hex(int(time()*10000000))[2:])

def can_i_run_software(software):
    
    #Test if external software is executable

    #:param software: list or string of executable(s)
    #:return: 1 (failure) or 0 (success)
    
    software_d = {}
    if type(software) == str:
        software = software.split(",")
    ex = 0
    for s in software:

        if s in ['picard']:
            command = [s, " MarkDuplicates", "-h"]
        elif s in ['fastp','gatk','samtools','GenomeAnalysisTK','hisat2','hisat2-build','bcftools']:
            command = [s, '--version']
        # elif s in []:
        #     command = [s, '-version']
        # else:
        #     command = [s, '--help']
        
        else:
            command = [s]
        try:
            sp = subprocess.run(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            log_coloredlogs("{} stdout: {}".format(s, sp.stdout.decode("utf-8").strip()))
            log_coloredlogs("{} stderr: {}".format(s, sp.stderr.decode("utf-8").strip()))

            sp = subprocess.run(["which", s],stdout= subprocess.PIPE,stderr= subprocess.PIPE)
            s_path = str(sp.stdout).split("'")[1].strip("\\n")
            software_d[s] = s_path
        except FileNotFoundError:
            logging.error('{} executable not found !,please add its path in your system env and rerun,'
            'optionally, you can add its absolute path in config file mannually'.format(s))
            software_d[s] = "null"
            ex = 1
    return ex,software_d

def can_i_run_software_(software_path):

    s = os.path.basename(software_path)

    x = 1
    if s in ['picard']:
        command = [software_path, "-h"]
    elif s in ['fastp','gatk','samtools','GenomeAnalysisTK','hisat2','hisat2-build','bcftools']:
        command = [software_path, '--version']
    else:
        command = [software_path]
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    except FileNotFoundError:
        x =0
    
    return x











"""
###    my note
1----------------------------------
    # codeml needs input otherwise it prompts the user for input, so a dummy
    # file is created
-------------------------------------
2----------------str(sp.stdout).split("'")[1].strip("\\n")

    # if not os.path.isdir(self.tmp):
            # raise NotADirectoryError(
                    # 'tmp directory {} not found!'.format(self.tmp))

    #click.format_filename(config_file)

    #if trigger_exit:
        logging.error("Please add the missing information to the configuration file and rerun the analysis. Exiting.")
        sys.exit(1)
    #log_coloredlogs(datetime.datetime.today().ctime())

--------------fc_wgd.py
_OUTPUT_BLAST_FILE_PATTERN = '{}.blast.tsv'
_OUTPUT_MCL_FILE_PATTERN = '{}.mcl.tsv'
_OUTPUT_KS_FILE_PATTERN = '{}.ks.tsv'
_PARALOGS_OUTPUT_DIR_PATTERN = 'wgd_{}'
_ORTHOLOGS_OUTPUT_DIR_PATTERN = 'wgd_{}_{}'

_TMP_BLAST = '{}.blast_tmp'
_TMP_KS = '{}.ks_tmp'
---------------------------------------------------------------

-------------------
log_coloredlogs('Removing tmp directory')
        shutil.rmtree(tmp_blast_paralogs)
----------------------

#not len(line) or line.startswith('#'): 
"""