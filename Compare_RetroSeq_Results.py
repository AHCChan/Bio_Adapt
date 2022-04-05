HELP_DOC = """
RETROSEQ RESULTS COMPARER
(version 1.0)
by Angelo Chan

A tool to compare the RetroSeq Calls output files against a coordinates file to
see how many TEs were correctly called.



USAGE:
    
    python27 



MANDATORY:
    


OPTIONAL:



EXAMPLES SCENARIO EXPLANATION:

EXAMPLES:

USAGE:
"""

NAME = "Compare_RetroSeq_Results.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

DIRMOD__EDIT = "__INSERTED"
FILEMOD__COORDS = "__POST_INSERT_COORDS.tsv"
FILEMOD__SIZES = "__POST_INSERT_SIZES.tsv"
FILEMOD__FASTA = ".fa"

CONFIG__ignore_bad_slicing = False
CONFIG__mismatch_handling = 0



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__Mode = 4 # Perfect/Family/Skip/Fail
DEFAULT__Output = 1 # Flags



# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *
from RetroSeq_Calls_File_Reader import *



# Enums ########################################################################

class CLASSIFICATION:
    PASS_FAIL=1
    PASS_SKIP_FAIL=2
    PERFECT_FAM_FAIL=3
    PERFECT_FAM_SKIP_FAIL=4

class OUTPUT:
    FLAGS=1
    COUNT=2
    FLAGS_POSITVE=3

class TIEBREAKER:
    FULL=1
    PARTIAL=2
    NONE=3

class FLAG:
    SUCCESS=1
    SKIP=2
    FAIL=3
    FAMILY=4
    TIED_SUCCESS=5
    TIED_FAMILY=6




# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"Compare_RetroSeq_Results.py -h"

STR__VCF_Invalid = """
ERROR: Unable to scan the following folder for VCF files:
    {D}"""

STR__error_no_VCF = """
ERROR: No VCF files detected in:
    {D}"""

STR__invalid_classification = """
ERROR: Invalid classifcation mode specified:
    {S}"""

STR__invalid_output = """
ERROR: Invalid output mode specified:
    {S}"""

STR__invalid_tiebreaker = """
ERROR: Invalid tiebreaker method specified:
    {S}"""






STR__metrics = """
                       Total TE Sites: {A}
                    Total Calls Files: {B}

     Highest Perfect Matches Per File: {C}
     Average Perfect Matches Per File: {D}

    Highest Partial+ Matches Per File: {E}
    Average Partial+ Matches Per File: {F}

     Highest Perfect Matches Per Site: {G}
     Average Perfect Matches Per Site: {H}

    Highest Partial+ Matches Per Site: {I}
    Average Partial+ Matches Per Site: {J}
"""

STR__compare_begin = "\nRunning Compare_RetroSeq_Results..."

STR__compare_complete = "\nCompare_RetroSeq_Results successfully finished."



STR__unexpected_failure = "\nProgram exited with an unexpected error."



# Lists ########################################################################

LIST__vcf_1 = [".vcf", ".Vcf", ".VCF"]



# Dictionaries #################################################################

DICT__Values = {
    CLASSIFICATION.PASS_FAIL = {
        OUTPUT.FLAGS = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 3,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        OUTPUT.COUNT = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        OUTPUT.FLAGS_POSITIVE = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PASS_SKIP_FAIL = {
        OUTPUT.FLAGS = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 3,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        OUTPUT.COUNT = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        OUTPUT.FLAGS_POSITIVE = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PERFECT_FAM_FAIL = {
        OUTPUT.FLAGS = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: 3},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 2,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        OUTPUT.COUNT = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        OUTPUT.FLAGS_POSITIVE = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 7},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 5},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PERFECT_FAM_SKIP_FAIL = {
        OUTPUT.FLAGS = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: 3},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 2,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        OUTPUT.COUNT = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        OUTPUT.FLAGS_POSITIVE = {
            TIEBREAKER.FULL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 7},
            TIEBREAKER.PARTIAL = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 5},
            TIEBREAKER.NONE = {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}}}



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Compare_RetroSeq_Results__STR(baseline_file, calls_folder, chromosomes_file,
        output_file, classification, output, tiebreaker):
    """
    """
    vcf_files = Get_VCF_From_Folder(calls_folder)
    if not vcf_files: return 1
    list_of_chrs = Get_Chrs_From_File(chromosomes_file)
    exit_code = Compare_RetroSeq_Results__LIST(baseline_file, vcf_files,
        list_of_chrs, output_file, classification, output, tiebreaker)
    return 0

def Get_VCF_From_Folder(dirpath):
    """
    """
    try:
        files = os.listdir(dirpath)
        files = [string for string in files if Is_VCF(string)]
        if len(files) == 0:
            PRINT.printE(STR__error_no_VCF.format(S=dirpath))
            return []
        files = [dirpath + "\\" + string for string in files]
        return files
    except:
        PRINT.printE(STR__VCF_Invalid)
        return []

def Is_VCF(filepath):
    """
    """
    if filepath[:-4] in LIST__vcf_1: return True
    return False

def Get_Chrs_From_File(filepath):
    """
    """
    result = []
    f = Table_Reader(filepath)
    f.Set_Delimiter("\t")
    f.Open()
    while not f.End():
        f.Read()
        c = f.Get_Current()[0]
        result.append(c)
    f.Close()
    return result    
    
def Compare_RetroSeq_Results__LIST(baseline_file, calls_files, list_of_chrs,
        output_file, classification, output, tiebreaker):
    """
    """
    PRINT.printP(STR__compare_begin)
    
    # Pre-calculate
    length = len(calls_files)
    range_ = range(length)
    
    # Setup reporting
    sites = 0
    site_max_perfect = 0
    site_max_partial = 0
    files_perfect = [0]*length
    files_partial = [0]*length
    
    # Setup the I/O
    b = Table_Reader(baseline_file)
    b.Set_Delimiter("\t")

    files = [RetroSeq_Calls_File_Reader(f) for f in calls_files]
    for f in files:
        f.Create_Chr_Order(list_of_chrs)
    
    o = open(output_file ,"w")
    
    # Main loop
    b.Open()
    while not b.End():
        sites += 1
        # Read
        b.Read()
        values = b.Get_Current()
        chr_ = values[0]
        start = int(values[1])
        end = int(values[2])
        name = values[3]
        family_1 = values[4]
        family_2 = values[5]
        coords = [chr_, start, end]
        # Setup
        perfect = 0
        partial = 0
        # For all RetroSeq Calls files
        for i in range_:
            f = files[i]
            f.Read_Until(coords)
            best = f.Get_Best()
            calls = len(best)
            count = 0
            flag = FLAG.SKIP
            if calls == 0: pass
            elif calls == 1:
                c_chr_, c_start, c_count, c_name = best[0]
                count = c_count
                if c_name == name: # Perfect
                    # Stats
                    perfect += 1
                    partial += 1
                    files_perfect[i] += 1
                    files_partial[i] += 1
                    # Flag
                    flag = FLAG.SUCCESS                 
                elif Partial_Match(name, family_1, family_2, c_name): # Partial
                    # Stats
                    partial += 1
                    files_partial[i] += 1
                    # Flag
                    flag = FLAG.FAMILY
            else:
                for call in best:
                    c_chr_, c_start, c_count, c_name = call
                    if c_name == name:
                        flag = FLAG.TIED_SUCCESS
                        count = c_count
                    elif Partial_Match(name, family_1, family_2, c_name):
                        if flag != FLAG.TIED_SUCCESS:
                            flag = FLAG.TIED_FAMILY
                            count = c_count
            # Values
            value = DICT__Values[classification][output][tiebreaker][flag]
            if output == OUTPUT.COUNT: value = value*count
            values.append(value)
        # Max stats
        if perfect > site_max_perfect: site_max_perfect = perfect
        if partial > site_max_partial: site_max_partial = partial
        # Write
        sb = "\t".join(values) + "\n"
        o.write(sb)
    
    # Close up
    b.Close()
    o.close()
    for f in files:
        f.Close()
    
    # Reporting
    
    # Wrap up
    return 0



def Partial_Match(original, family_1, family_2, called):
    """
    """
    if((called in original) or (original in called) or
            (called in family_1) or (family_1 in called) or
            (called in family_2) or (family_2 in called)):
        return True
    return False



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Compare_RetroSeq_Results(raw_command_line_input):
    """
    Parse the command line input and call the Compare_RetroSeq_Results__STR
    function with appropriate arguments if the command line input is valid.
    """
    pass
    


# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Compare_RetroSeq_Results(sys.argv)


