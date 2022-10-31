import logging
import re
import os
import string
from utils import log_coloredlogs
import sys

logging.basicConfig(
                #filename="./logfile",
                filemode='a',
                format='%(asctime)s: %(levelname)s\t%(message)s',
                datefmt='%H:%M:%S',
                level=logging.DEBUG
        )
def _snp_dp_filter(ind,dp = 5,del_ffile=False):
    cwd = os.getcwd()
    myindir = os.path.join(cwd,ind)
    myoutdir = os.path.join(cwd,"temp_DPfiltered")
    try:
        log_coloredlogs("Creating depth-filter folder in %s......" %myoutdir)
        os.mkdir(myoutdir)    
    except FileExistsError as e:
        #print("The output path existing, the file will be rewrited")
        log_coloredlogs("The output path of temp_DPfiltered existing, the old file will be rewrited")
    # except Exception as e:
    #     logging.exception(e)
        
    if os.path.isdir(myindir):
        logging.debug("File exist in %s, will use the existed data" %myindir)
    else:
        log_coloredlogs("The vcf file is not exist in %s, please check your data." %myindir,err=True)
        sys.exit()
        #logging.error("The vcf file is not exist, please check your data")
    in_list = os.listdir(myindir)
    #all_dict = {}
    log_coloredlogs("Filtering raw vcf file......")
    for inf in in_list:
        f_list = []
        with open(os.path.join(myindir,inf),"r") as in_f:
            with open(os.path.join(myoutdir,inf+".f"),"w") as out_f:           
                for lines in in_f:
                    if lines.startswith("#"):
                        f_list.append(lines)
                        out_f.write(lines)
                    else:
                        dp4 = re.search("DP4=(\d+),(\d+),(\d+),(\d+);",lines).groups()
                        len_dp4 = int(dp4[0]) + int(dp4[1]) +int(dp4[2]) + int(dp4[3])
                        if len_dp4 >= dp:
                            f_list.append(lines)
                            out_f.write(lines)
        log_coloredlogs("%-10s records retaind after DPfiltering for %s" %(str(len(f_list)),inf))
        log_coloredlogs("All file were saved in %s" %myoutdir)
        #all_dict[inf] = f_list

    return(myoutdir)
"""
if __name__ == '__main__':
    ind = "test"
    outpath = snp_dp_filter(ind)
    print(outpath)
"""

def _parse_gff(gff):
    gff_d = {}
    log_coloredlogs("Parsing your gff file......")
    with open(gff,"r") as f:
        for lines in f:
            if lines.startswith("#"):
                pass
            else:
                k = re.split("\t+",lines.strip())[0]
                v = lines.strip()
                if k in gff_d:
                    gff_d[k].append(v)    
                else:
                    gff_d[k] = [v]
    return(gff_d)

def vcf2csv(gff,ind,percent_alle,dp = 5,feature="CDS,three_prime_UTR,five_prime_UTR",del_ffile=False):
    #loading gff
    gff_f = _parse_gff(gff)
    #loading vcffilter file
    ffile = _snp_dp_filter(ind,dp)
    #creating folder
    cwd = os.getcwd()
    myoutdir = os.path.join(cwd,"snp_statistics")
    try:
        log_coloredlogs("Creating vcf2csv folder in %s......"%myoutdir)
        os.mkdir(myoutdir)
    except FileExistsError as e:
        #print("The output path existing, the file will be overwrited")
        log_coloredlogs("The output path of snp_statistics existing, the old file will be overwrited")
        #os.rmdir(myoutdir)
    
    ffile_l = os.listdir(ffile)
    log_coloredlogs("Transforming filtered vcf file to csv......\nThis may take some time.")
    for filtered_f in ffile_l:
        with open(os.path.join(myoutdir,filtered_f+".csv"),"w") as fo:
            fo.write("Chr\tFeature\tPosition\tDP_alt\tDP_ref\tAlt_percent\tEditing\tDescription\n")
            with open(os.path.join(ffile,filtered_f)) as f:
                f_list = f.readlines()
                for lines in f_list:
                    lines_list = re.split("\t+",lines.strip())

                    if lines.startswith("#"):
                        pass                
                    else:
                        dp4 = re.search("DP4=(\d+),(\d+),(\d+),(\d+);",lines).groups()
                        dp_ref = int(dp4[0]) + int(dp4[1])
                        dp_alle = int(dp4[2]) + int(dp4[3])
                        len_dp4 = dp_alle + dp_ref
                        perc_alle = dp_alle/len_dp4
                        if perc_alle >= percent_alle:
                            feature_list = feature.strip().split(",")
                            flag = 0
                            for gff_records in gff_f[lines_list[0]]:
                                if gff_records.startswith("#"):
                                    pass
                                else:
                                    gff_records_list = re.split("\t+",gff_records)
                                    if  gff_records_list[2] in feature_list:
                                        if gff_records_list[3] <= lines_list[1] and gff_records_list[4] >= lines_list[1]:
                                            flag = 1
                                            base_r = lines_list[3]
                                            base_a = lines_list[4]
                                            if gff_records_list[6] == "-":
                                                base_r = base_r.translate(str.maketrans("ACGTacgt","TGCAtgca"))
                                                base_a = base_a.translate(str.maketrans("ACGTacgt","TGCAtgca"))
                                            base_r = base_r.replace("T","U")
                                            base_a = base_a.replace("T","U")
                                            fo.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%s to %s\t%s\n" %(gff_records_list[0],gff_records_list[2],lines_list[1],str(dp_alle),str(dp_ref),perc_alle,base_r,base_a,gff_records_list[8]))
                            if flag == 0:
                                base_r = lines_list[3]
                                base_a = lines_list[4]
                                base_r = base_r.replace("T","U")
                                base_a = base_a.replace("T","U")
                                fo.write("%s\tMissMatch\t%s\t%s\t%s\t%.2f\t%s to %s\tMissMatch\n" %(lines_list[0],lines_list[1],str(dp_alle),str(dp_ref),perc_alle,base_r,base_a))
    log_coloredlogs("Completed !")
    if del_ffile:
        log_coloredlogs("Will delete the filtered vcf file(to save disk) %s" %ffile,warn=True)
        os.rmdir(ffile)
    
    return(myoutdir)


if __name__ == '__main__':
    gff = "/home/whosy/workstation/rnaediting/bcftools/microTom.gff3"
    ind = "/home/whosy/workstation/rnaediting/test3"
    percent_alle = 0.3
    test = vcf2csv(gff,ind,percent_alle,dp = 5,feature="mRNA,CDS,three_prime_UTR,five_prime_UTR",del_ffile=False)
    print(test)
