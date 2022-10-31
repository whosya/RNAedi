from email.policy import default
from genericpath import exists
from importlib.resources import path
from rnaediting.utils import log_subprocess,init_logging,log_coloredlogs
import click
import logging
from sys import argv
import re
from rnaediting._version import __version__

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(__version__, prog_name="RNAedi", help="Print version number.")
def cli():
    """
\b    
_____ _   _   ___           _ _ 
| ___ \ \ | | / _ \         | (_)
| |_/ /  \| |/ /_\ \   ___  __| |_ 
|    /| . ` ||  _  |  / _ \/ _` | |
| |\ \| |\  || | | | |  __/ (_| | |
\_| \_\_| \_/\_| |_/  \___|\__,_|_|  

\b
Welcome to RNAedi!
    """
@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Generates configuration file.")

@click.argument('rawreads_path', type=click.Path(exists=True), required=True)
@click.option('--configfile_name', type=click.Path(),default="conf.txt")
@click.option('--single',type=click.BOOL,default=False,  help="single- or paired-end (default:paired)")
@click.option('--custom_file',type=click.Path(exists=True),default=None,  help="User-defined sample names (with comma as separator) (I.E. 'sample' in [sample_1.fq]")
@click.option('--strand',type=str,default="1,2",  help="User-defined strand symbol in rawreads file names (I.E. '1' in [sample_1.fq]")
@click.option('--separator',type=str,default="_",  help="User-defined separator in  (I.E. '.fq' in [sample_1.fq]")
@click.option('--sufix',type=str,default=None,  help="User-defined common sufix (I.E. '.fq' in [sample_1.fq]")

def generate_config(rawreads_path, configfile_name, single, custom_file, strand, separator, sufix):
    # """
    # Generates the configuration file for the rnaediting.
    # The configuration file name is given by argument FILENAME.

    # rawreads_path: the path of your rawreads path
    # configfile_name: configuration file name
    # single: a bool value to specify if your rna-seq datas are single-end or not. default[False] 
    # custom_file:
    # strand:
    # separator:
    # sufix:
    # """
    from rnaediting.generate_config import generate_configfile
    
    generate_configfile(rawreads_path, configfile_name, single, custom_file, strand, separator, sufix)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Run bcftools workflow.")

@click.argument('rawreads_path', type=click.Path(exists=True), required=True)
@click.argument('genome', type=click.Path(exists=True), required=True)
@click.option('-T','--threads', type=int, default=4,  help="Threads used in fastp,hisat2,..... (default:4)")
@click.option('--conf', default='conf.txt', type=click.Path(exists=True),  help="User-defined path to the configuration file (default:'conf.txt')")
@click.option('--other_fastp_options', type=str, default="", help="other options for fastp, in default only neccessary parameters[-i,-I,-o,-O,-j,-h,-q] were used.(comma separating,E.G. '-z=2,-A')")
@click.option('--fastp2', is_flag=True, help="Use this option when your rawdata not filtered (default:True")
@click.option('--markdup', is_flag=True, help="Perform markduplication or not (default:True)")
@click.option('--rg', type=str, default="WES",  help="Set the read group ID (default:WES)")
@click.option('--lb', type=str, default="WES",  help="DNA preparation library identifier (default:WES)")
@click.option('--pl', type=str, default="ILLUMINA",  help="Platform/technology used to produce the read (default:ILLUMINA)")
@click.option('--other_hisat2_options', type=str, default="", help="other options for hisat2, in default only neccessary parameters[--rg-id,--rg,-x,-1,-2,-p,-U] were used.(comma separating)")

def BcfSnp(rawreads_path, genome, threads, conf, other_fastp_options, fastp2, markdup, rg, lb, pl, other_hisat2_options):
    """
    Call snp using bcftools workflow.

    Example:

    RNAedi  bcfsnp testdata/ assembly.fasta --fastp2  --markdup -T 12  --other_hisat2_options -t,--new-summary --other_fastp_options -z=2,-A

     
    """
    fastp_options = re.split(",|=",other_fastp_options.strip())
    hisat2_options = re.split(",|=",other_hisat2_options.strip())

    BcfSnp_(rawreads_path=rawreads_path, genome=genome, other_fastp_options=fastp_options, other_hisat2_options=hisat2_options, threads=threads, conf=conf, fastp_1=fastp2, markdup=markdup, rg=rg, lb=lb, pl=pl)

def BcfSnp_(rawreads_path, genome, other_fastp_options, other_hisat2_options, threads=4, conf='conf.txt', fastp_1=True, markdup=True, rg="WES", lb="WES", pl="ILLUMINA"):

    from rnaediting.run_bcftools_pip import genome_index, Runbcftools, Callvar
    '''
    step1 indexing genome, 
    input: genome, 
    output: 1_genome_index
    '''
    out1 = genome_index(genome,config=conf,threads=threads)
    log_coloredlogs("all index files were saved in %s" %out1)

    '''
    step2 and step3 filtering raw reads(optional) and align, 
    step2 input:rawdata, 
          output:0_cleaned_data 
    step3 input:2_aligned_data or rawdata, 
          output:2_aligned_data
    '''
    if fastp_1:
        ######################################################################################################!!!!!!!!!!!!!!!!!!! if *argv were put in ahead,location para should be specified with shican
        out2 = Runbcftools().run_fastp(*other_fastp_options,raw_file=rawreads_path,threads=threads)
        log_coloredlogs("All filtered files were saved in %s" %(out2))
        out3 = Runbcftools().run_hisat2(*other_hisat2_options, genome=genome, fastp=fastp_1, hisat2="hisat2",  samtools="samtools", threads=threads, ind = out2, outdir="2_aligned_data", rg="WES",lb="WES",pl="ILLUMINA")
        log_coloredlogs("all aligned files were saved in %s" %out3)

    else:
        out2 = rawreads_path
        print(out2)
        out3 = Runbcftools().run_hisat2(genome=genome, fastp=fastp_1, hisat2="hisat2", samtools="samtools", threads=threads, ind = out2, outdir="2_aligned_data", rg="WES",lb="WES",pl="ILLUMINA")
        log_coloredlogs("all aligned files were saved in %s" %out3)

    '''
    step4 and step5 markduplication and call snp,
    step4 input: 2_aligned_data, 
          output: 3_marked_data
    step5 input: 3_marked_data or 2_aligned_data, 
          output: temp_mpilup and 4_callsnp_result
    '''
    if markdup:
        out4 = Runbcftools().run_markdup()
        log_coloredlogs("All marked files were saved in %s\n" %(out4))
        out5 = Callvar(genome, markdup=True, conf="conf.txt", bcftools="bcftools", ind=out4, temp_outdir="temp_mpilup", threads=threads).run_bcftools_snp()
        log_coloredlogs("All snp files were saved in %s" %out5)
    else:
        out4 = out3
        out5 = Callvar(genome, markdup=False, conf="conf.txt", bcftools="bcftools", ind=out4, temp_outdir="temp_mpilup", threads=threads).run_bcftools_snp()
        log_coloredlogs("All snp files were saved in %s" %out5)



# @click.argument('rawreads_path', type=click.Path(exists=True), required=True)
# @click.argument('genome', type=click.Path(exists=True), required=True)
# @click.option('-T','--threads', type=int, default=4,  help="Threads used in fastp,hisat2,..... (default:4)")
# @click.option('--conf', default='conf.txt', type=click.Path(exists=True),  help="User-defined path to the configuration file (default:'conf.txt')")
# @click.option('--markdup', is_flag=True, help="Perform markduplication or not (default:True)")


    




