import os
import logging
import sys
from .utils import log_coloredlogs



class check_file:
    def __init__(self,ind):
        self.ind = ind
    
    def check_dir_nonexistent_or_empty(self):
        """
        Returns absolute path, exiting when the provided folder path doesn't exist or is empty.
        """
        input_dir = os.path.abspath(self.ind)
        if os.path.exists(input_dir):
            logging.debug("File existed in %s, will use the existed data" %self.ind)
            if len(os.listdir(input_dir)) == 0:
                log_coloredlogs("Empty folder!",err=True)
                sys.exit()
        else:
            log_coloredlogs("The file is not exist in %s, please check your data" %self.ind,err=True)
            sys.exit()
        
        return input_dir
    
    def check_file_nonexistent_or_empty(self):
        """
        Returns an error message when the provided file doesn't exist or is empty.
        """
        input_file = os.path.abspath(self.ind)
        if not os.path.exists(input_file):
        # print error message separate to ensure print out
            log_coloredlogs("File %s not exist, please check your data" %self.ind,err=True)
            sys.exit()
        #     return True
        # else:
        #     return False
        if os.path.getsize(input_file) == 0:
            # print error message separate to ensure print out
            log_coloredlogs("File %s empty, please check your data" %self.ind,err=True)
            sys.exit()
        #     return True
        # else:
        #     return False
        return input_file
    def check_file_nonexistent_or_empty_(self):
        """
        return: a boolean stating whether the file doesn't exist or is empty (True) 
        or whether it does exist and is not empty (False)
        """
        input_file = os.path.abspath(self.ind)
        if not os.path.exists(input_file):
        # print error message separate to ensure print out
            # log_coloredlogs("File %s not exist, please check your data" %self.ind,err=True)
            # sys.exit()
            return True
        else:
            return False
        if os.path.getsize(input_file) == 0:
            # print error message separate to ensure print out
            # log_coloredlogs("File %s empty, please check your data" %self.ind,err=True)
            # sys.exit()
            return True
        else:
            return False
    def check_file_from_lastrun(self,costom_name=None):
        """
        Checking the output file of certain function, if the same file exists (e.g. produced from last breaked run),
        it will be overwrited !
        """
        output_path = os.path.abspath(self.ind)
        try:
            logging.info("Creating %s folder in %s......" %(costom_name,output_path))
            os.mkdir(output_path)    
        except FileExistsError as e:
            #print("The output path existing, the file will be rewrited")
            logging.info("The output path of %s existing in %s, the old file will be overwrited" %(costom_name,output_path))
        
        return output_path
#    def check_gff(self):
