HELP_DOC = """
RETROSEQ RESULTS COMPARER
(version 1.0)
by Angelo Chan

A tool to compare the RetroSeq Calls output files against a coordinates file to
see how many TEs were correctly called.



USAGE:
    
    python27 Compare_RetroSeq_Results.py <baseline_file> <calls_folder>
            <chromosomes_file> [-o <output_file>] [-c <classification>]
            [-v <values>] [-t <tiebreaker>]



MANDATORY:
    
    baseline_file
        
        The filepath of the input file containing the original TEs coordinates.
        The file is assumed to be a six-column TSV, with the columns containing
        the following:
            1:  Chromosome name (part of the genomic coordinates)
            2:  Start index (part of the genomic coordinates)
            3:  End index (part of the genomic coordinates)
            4:  TE name
            5:  TE family (1st)
            6:  TE famile (2nd)
    
    calls_folder
        
        The directory path of the folder containing the RetroSeq results. These
        VCF files are RetroSeq Call output files after undergoing through a
        Bedtools Sort.
    
    chromosomes_file
        
        The filepath of a chromosomes list file. This is a TSV containing at
        least one column, with the first column containing the names of the
        chromosomes in sorted order.



OPTIONAL:
    
    output_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file where the results will be output to.
    
    classification
        
        (DEFAULT: 4 - Pass/Family/Skip/Fail)
        
        The way the resuts are characterized. The modes are:
            1:  Pass/Fail
            2:  Pass/Skip/Fail
            3:  Pass/Family/Fail
            4:  Pass/Family/Skip/Fail
        
        The terms mean the following:
            Pass:   The correct TE subfamily was called.
            Family: The correct TE subfamily was not called, but the correct
                    TE family was called.
            Skip:   No calls were made.
            Fail:   None of the above. (Typically an incorect call. In Pass/Fail
                    mode, Skips are classified as Fails.
    
    values
        
        (DEFAULT: 1 - Flags)
        
        The values written to the output file. There 3 options are:
            
            1:  Flags - Set values indicating what the result was, based on
                    whether the best call was Pass/Family/Skip/Fail, and whether
                    or not there was a tie. These values are centered around 0.
                        
            2:  Counts - The number of reads identified by RetroSeq as
                    supporting a particular call)
                        
            3:  Positive Flags - Set values indicating what the result was, based on
                    whether the best call was PASS/FAMILY/SKIP/FAIL, and whether
                    or not there was a tie. These values are all positive.
        
        The flags, in ascending order, indicate:
            FAIL:       The most supported call(s) was/were all completely wrong.
            PASS:       No call was made which overlapped with the original TE's
                        position.
            TIE-FAMILY: Multiple calls were made, at least one of which was
                        the correct family of the TE at that position.
            TIE-PASS:   Multiple calls were made, at least one of which was
                        the correct TE at that position.
            FAMILY:     Multiple calls were made, at least one of which was the
                        correct family of the TE at that position.
            PASS:       Multiple calls were made, at least one of which was the
                        correct TE at that position.
        
        The flag values for Option 1 (Flags) are:
            -4, 0, 1, 2, 3, 4
        
        The flag values for Option 3 (Positive Flags) are:
            0, 4, 5, 6, 7, 8
    
    tiebreaker
        
        (DEFAULT: 2 - Partial)
        
        How multiple calls, with the same number of supporting reads, will be
        treated. There are 3 options:
            1:  (Full) If at least one of the calls is the correct Family or
                Subfamily, that TE will be given a full score.
            2:  (Partial) If at least one of the calls is the correct Family or
                Subfamily, that TE will be given a partial score.
            3:  (None) Ties will be treated as FAILs.

    

EXAMPLES SCENARIO EXPLANATION:
    
    1:
    A harsh pass/fail test. There are only two possible outcomes - successfully
    calling a TE, or failing to, either with a miscall, or with nothing called
    at the location of the original TE. Focuses exclusviely on the True Positive
    rate.
    
    2:
    A slightly less harsh pass/skip/fail test. Useful for gauging the extent of
    False Positive calling.
    
    3:
    A pass/family/skip/fail test with partial tiebreakers. Useful for gauging
    the exact degree of success/failure with imperfect calls.
    
    4:
    A lenient pass/family/skip/fail test. Full points are given in the event of
    a tie. At the time this software was written, no use-case scenario was
    envisioned for this kind of analysis, but the option is there.

EXAMPLES:
    
    python27 Compare_RetroSeq_Results.py Path/rmsk_rearranged.tsv Path/Calls
            Path/chr_sizes.tsv -o Path/Harsh_Compare.tsv -c 1 -v 1 -t 3
    
    python27 Compare_RetroSeq_Results.py Path/rmsk_rearranged.tsv Path/Calls
            Path/chr_sizes.tsv -o Path/Normal_Compare.tsv -c 2 -v 1 -t 3
    
    python27 Compare_RetroSeq_Results.py Path/rmsk_rearranged.tsv Path/Calls
            Path/chr_sizes.tsv -o Path/Detailed_Compare.tsv -c 4 -v 1 -t 2
    
    python27 Compare_RetroSeq_Results.py Path/rmsk_rearranged.tsv Path/Calls
            Path/chr_sizes.tsv -o Path/Lenient_Compare.tsv -c 4 -v 1 -t 1

USAGE:
    
    python27 Compare_RetroSeq_Results.py <baseline_file> <calls_folder>
            <chromosomes_file> [-o <output_file>] [-c <classification>]
            [-v <values>] [-t <tiebreaker>]
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

FILEMOD = "__RetroSeq_Compare.tsv"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__Classification = 4 # Perfect/Family/Skip/Fail
DEFAULT__Values = 1 # Flags
DEFAULT__Tiebreaker = 2 # Partial



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

class VALUES:
    FLAGS=1
    COUNT=2
    FLAGS_POSITIVE=3

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

STR__invalid_values = """
ERROR: Invalid output values mode specified:
    {S}"""

STR__invalid_tiebreaker = """
ERROR: Invalid tiebreaker method specified:
    {S}"""

STR__insufficient_columns = """
ERROR: Insufficient columns in coordinates file. The coordinates file should
contain at least 6 columns, the first 3 of which contain genome coordinates,
the 4th of which should contain the name of the TE, while the 5th and 6th should
contain the names of the TE family or grouping."""




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

LIST__PF = ["PF", "pf", "PASSFAIL", "PassFail", "passfail", "1"]
LIST__PSF = ["PSF", "psf", "PASSSKIPFAIL", "PassSkipFail", "passskipfail", "2"]
LIST__PFF = ["PFF", "pff", "PASSFAMILYFAIL", "PassFamilyFail", "passfamilyfail",
        "PASSFAMFAIL", "PassFamFail", "passfamfail", "3"]
LIST__PFSF = ["PFSF", "pfsf", "PASSFAMILYSKIPFAIL", "PassFamilySkipFail",
        "passfamilyskipfail", "PASSFAMSKIPFAIL", "PassFamSkipFail",
        "passfamskipfail", "4"]

LIST__FLAGS = ["F", "f", "FLAGS", "Flags", "flags", "FLAG", "Flag", "flag", "1"]
LIST__COUNTS = ["C", "c", "COUNTS", "Counts", "counts", "COUNT", "Count",
        "count", "2"]
LIST__P_FLAGS = ["P", "p", "PF", "pf", "POSITIVEFLAGS", "PositiveFlags",
        "positiveflags", "POSITIVEFLAG", "PositiveFlag", "positiveflag", "3",
        "+"]

LIST__FULL = ["F", "f", "FULL", "Full", "full", "1"]
LIST__PARTIAL = ["P", "p", "PARTIAL", "Partial", "partial", "2"]
LIST__NONE = ["N", "n", "NONE", "None", "none", "3"]



# Dictionaries #################################################################

DICT__Values = {
    CLASSIFICATION.PASS_FAIL: {
        VALUES.FLAGS: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 3,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        VALUES.COUNT: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        VALUES.FLAGS_POSITIVE: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PASS_SKIP_FAIL: {
        VALUES.FLAGS: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: 3,
                FLAG.TIED_FAMILY: -4},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: -4,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        VALUES.COUNT: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        VALUES.FLAGS_POSITIVE: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 0},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 0,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PERFECT_FAM_FAIL: {
        VALUES.FLAGS: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: 3},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 2,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: -4,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        VALUES.COUNT: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        VALUES.FLAGS_POSITIVE: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 7},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 5},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}}},
    CLASSIFICATION.PERFECT_FAM_SKIP_FAIL: {
        VALUES.FLAGS: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 4,
                FLAG.TIED_FAMILY: 3},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: 2,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 4,
                FLAG.SKIP: 0,
                FLAG.FAIL: -4,
                FLAG.FAMILY: 3,
                FLAG.TIED_SUCCESS: -4,
                FLAG.TIED_FAMILY: -4}},
        VALUES.COUNT: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 1,
                FLAG.TIED_FAMILY: 1},
            TIEBREAKER.NONE: {
                FLAG.SUCCESS: 1,
                FLAG.SKIP: 0,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 1,
                FLAG.TIED_SUCCESS: 0,
                FLAG.TIED_FAMILY: 0}},
        VALUES.FLAGS_POSITIVE: {
            TIEBREAKER.FULL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 8,
                FLAG.TIED_FAMILY: 7},
            TIEBREAKER.PARTIAL: {
                FLAG.SUCCESS: 8,
                FLAG.SKIP: 4,
                FLAG.FAIL: 0,
                FLAG.FAMILY: 7,
                FLAG.TIED_SUCCESS: 6,
                FLAG.TIED_FAMILY: 5},
            TIEBREAKER.NONE: {
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
        output_file, classification, output_values, tiebreaker):
    """
    Compare the results of multiple datasets, which have gone through the
    RetroSeq pipeline, against the original TE coordinates and output the
    comparison results into an output file.
    
    @baseline_file
            (str - filepath)
            The filepath of the input file containing the original TEs
            coordinates. The file is assumed to be a six-column TSV, with the
            columns containing the following:
                1:  Chromosome name (part of the genomic coordinates)
                2:  Start index (part of the genomic coordinates)
                3:  End index (part of the genomic coordinates)
                4:  TE name
                5:  TE family (1st)
                6:  TE famile (2nd)
    @calls_folder
            (str - dirpath)
            The directory path of the folder containing the RetroSeq results.
            These VCF files are RetroSeq Call output files after undergoing
            through a Bedtools Sort.
    @chromosomes_file
            (str - filepath)
            The filepath of a chromosomes list file. This is a TSV containing
            at least one column, with the first column containing the names of
            the chromosomes in sorted order.   
    @output_file
            (str - filepath)
            The filepath of the output file where the results will be output
            to.
    @classification
            (int) - Pseudo ENUM
            The way the resuts are characterized. The modes are:
                1:  Pass/Fail
                2:  Pass/Skip/Fail
                3:  Pass/Family/Fail
                4:  Pass/Family/Skip/Fail
    @output_values
            (int) - Pseudo ENUM
            The values written to the output file. There 3 options are:
                1:  Flags
                2:  Counts
                3:  Positive Flags
            The flags, in ascending order, indicate:
                FAIL, PASS, TIE-FAMILY, TIE-PASS, FAMILY, PASS
            The flag values for Option 1 (Flags) are:
                -4, 0, 1, 2, 3, 4
            The flag values for Option 3 (Positive Flags) are:
                0, 4, 5, 6, 7, 8
    @tiebreaker
            (int) - Pseudo ENUM
            How multiple calls, with the same number of supporting reads, will
            be treated. There are 3 options:
                1:  Full points (treated the same as a perfect call)
                2:  Partial points
                3:  No points (treated the same as a failed call)
    
    Compare_RetroSeq_Results__STR(str, str, str, str, int, int, int) -> int
    """
    vcf_files = Get_VCF_From_Folder(calls_folder)
    if not vcf_files: return 1
    list_of_chrs = Get_Chrs_From_File(chromosomes_file)
    exit_code = Compare_RetroSeq_Results__LIST(baseline_file, vcf_files,
        list_of_chrs, output_file, classification, output_values, tiebreaker)
    return 0

def Get_VCF_From_Folder(dirpath):
    """
    Return a list of VCF files (with complete filepathing) from the given
    directory path. This function relies on file extensions to judge whether or
    not a file is a VCF, rather than examining the files' contents.
    
    Get_VCF_From_Folder(str) -> list<str>
    """
    try:
        files = os.listdir(dirpath)
        files = [string for string in files if Is_VCF(string)]
        if len(files) == 0:
            PRINT.printE(STR__error_no_VCF.format(D=dirpath))
            return []
        files = [dirpath + "\\" + string for string in files]
        return files
    except:
        PRINT.printE(STR__VCF_Invalid.format(D=dirpath))
        return []

def Is_VCF(filepath):
    """
    Return whether or not a string (presumably a file name or filepath) ends in
    a known VCF file extension.
    
    Is_VCF(str) -> bool
    """
    if filepath[-4:] in LIST__vcf_1: return True
    return False

def Get_Chrs_From_File(filepath):
    """
    Return a list of chromosomes from a TSV file containing the names of said
    chromosomes in the first column.
    
    Get_Chrs_From_File(str) -> list<str>
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
        output_file, classification, output_values, tiebreaker):
    """
    Identical to the function {Compare_RetroSeq_Results__STR}, except that
    function takes a filepath for the chromosomes list file as an input, while
    this function takes a list of chromosomes instead for that input.
    
    Compare_RetroSeq_Results__LIST(str, str, list<str>, str, int, int, int)
            -> int
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
    
    # Setup the I/O (1)
    b = Table_Reader(baseline_file)
    b.Set_Delimiter("\t")

    # Header (1)
    headers = ["Chr", "Start", "End", "Name", "Family1", "Family2"]
    b.Open()
    b.Read()
    values = b.Get_Current()
    length = len(values)
    extras = length-6
    b.Close()
    if extras < 0:
        PRINT.printE(STR__insufficient_columns)
        return 1
    
    # Setup the I/O (2)
    files = [RetroSeq_Calls_Reader(f) for f in calls_files]
    for f in files:
        f.Create_Chr_Order(list_of_chrs)
    for f in files:
        f.Open()
    o = open(output_file ,"w")
    
    # Header (2)
    while extras > 0:
        headers.append("#")
        extras -= 1
    for f in calls_files:
        filename = Get_File_Name(f)
        headers.append(filename)
    sb = "\t".join(headers) + "\n"
    o.write(sb)
    
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
        name = values[3].upper()
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
                c_name = c_name.upper()
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
                    flag = FLAG.FAIL
            else:
                for call in best:
                    c_chr_, c_start, c_count, c_name = call
                    c_name = c_name.upper()
                    if c_name == name:
                        flag = FLAG.TIED_SUCCESS
                        count = c_count
                    elif Partial_Match(name, family_1, family_2, c_name):
                        if flag != FLAG.TIED_SUCCESS:
                            flag = FLAG.TIED_FAMILY
                            count = c_count
            # Values
            value=DICT__Values[classification][output_values][tiebreaker][flag]
            if output_values == VALUES.COUNT: value = value*count
            values.append(str(value))
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
    
    PRINT.printP(STR__compare_complete)
    
    # Reporting
    Report_Metrics(sites, site_max_perfect, site_max_partial, files_perfect,
            files_partial)
    sites = 0
    site_max_perfect = 0
    site_max_partial = 0
    files_perfect = [0]*length
    files_partial = [0]*length
    
    # Wrap up
    return 0



def Partial_Match(original, family_1, family_2, called):
    """
    Return whether or not there was a partial match, defined as either the
    query string (from the RetroSeq call) containing either of the other three
    strings (from the original coordinates file) as substrings, or either of
    the three original strings containing the query string as a substring.
    
    @original
            (str)
            The original TE element name.
    @family_1
            (str)
            The first of two strings from the original annotations for the TE
            element denoting its family.
    @family_2
            (str)
            The second of two strings from the original annotations for the TE
            element denoting its family.
    @called
            (str)
            The query string. Outputted by RetroSeq.
    
    Partial_Match(str, str, str, str) -> bool
    """
    if((called in original) or (original in called) or
            (called in family_1) or (family_1 in called) or
            (called in family_2) or (family_2 in called)):
        return True
    return False



def Report_Metrics(sites, site_max_perfect, site_max_partial, files_perfect,
            files_partial):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @sites
            (int)
            The number of TEs in the original file.
    @site_max_perfect
            (int)
            The number of files where a site was called correctly, for the site
            with the highest number of files where it was called correctly.
    @site_max_perfect
            (int)
            The number of files where a site was called either completely
            correct or partially correct, for the site with the highest number
            of files where it was called either completely correct or partially
            correct.
    @files_perfect
            (list<int>)
            A list of the number of perfect calls. Each integer corresponds to
            a file.
    @files_partial
            (list<int>)
            A list of the number of perfect or partial calls. Each integer
            corresponds to a file.
    
    Report_Metrics(int, int, int, list<int>, list<int>) -> None
    """
    # Basic numbers
    size = len(files_perfect)
    length = float(size)
    sum_perf = sum(files_perfect)
    sum_part = sum(files_partial)
    # Calculations
    max_perf_per_file = max(files_perfect)
    avg_perf_per_file = sum_perf/length
    max_part_per_file = max(files_partial)
    avg_part_per_file = sum_part/length
    avg_perf_per_site = float(sum_perf)/sites
    avg_part_per_site = float(sum_part)/sites
    # Strings
    sites = str(sites)
    size = str(size)
    max_perf_per_file = str(max_perf_per_file)
    avg_perf_per_file = str(avg_perf_per_file) + "0"
    max_part_per_file = str(max_part_per_file)
    avg_part_per_file = str(avg_part_per_file) + "0"
    site_max_perfect = str(site_max_perfect)
    avg_perf_per_site = str(avg_perf_per_site) + "0"
    site_max_partial = str(site_max_partial)
    avg_part_per_site = str(avg_part_per_site) + "0"
    # Pad
    max_size = Get_Max_Len([sites, size, max_perf_per_file, avg_perf_per_file,
            max_part_per_file, avg_part_per_file, site_max_perfect,
            avg_perf_per_site, site_max_partial, avg_part_per_site])
    sites = Pad_Str(sites, max_size, " ", 0)
    size = Pad_Str(size, max_size, " ", 0)
    max_perf_per_file = Pad_Str(max_perf_per_file, max_size, " ", 0)
    avg_perf_per_file = Pad_Str(avg_perf_per_file, max_size, " ", 0)
    max_part_per_file = Pad_Str(max_part_per_file, max_size, " ", 0)
    avg_part_per_file = Pad_Str(avg_part_per_file, max_size, " ", 0)
    site_max_perfect = Pad_Str(site_max_perfect, max_size, " ", 0)
    avg_perf_per_site = Pad_Str(avg_perf_per_site, max_size, " ", 0)
    site_max_partial = Pad_Str(site_max_partial, max_size, " ", 0)
    avg_part_per_site = Pad_Str(avg_part_per_site, max_size, " ", 0)
    # Print
    PRINT.printM(STR__metrics.format(A = sites, B = size, C = max_perf_per_file,
            D = avg_perf_per_file, E = max_part_per_file, F = avg_part_per_file,
            G = site_max_perfect, H = avg_perf_per_site, I = site_max_partial,
            J = avg_part_per_site))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Compare_RetroSeq_Results(raw_command_line_input):
    """
    Parse the command line input and call the Compare_RetroSeq_Results__STR
    function with appropriate arguments if the command line input is valid.
    """
    PRINT.printP(STR__parsing_args)
    # Remove the runtime environment variable and program name from the inputs
    inputs = Strip_Non_Inputs(raw_command_line_input, NAME)
    
    # No inputs
    if not inputs:
        PRINT.printE(STR__no_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Help option
    if inputs[0] in LIST__help:
        print(HELP_DOC)
        return 0
    
    # Initial validation
    if len(inputs) < 3:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_in_baseline = inputs.pop(0)
    valid = Validate_Read_Path(path_in_baseline)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in_baseline))
        PRINT.printE(STR__use_help)
        return 1
    b = Table_Reader(path_in_baseline)
    b.Set_Delimiter("\t")
    b.Open()
    b.Read()
    test_row = b.Get_Current()
    if len(test_row) < 6:
        PRINT.printE(STR__insufficient_columns)
        PRINT.printE(STR__use_help)
        return 1
    #
    path_in_folder = inputs.pop(0)
    valid = Validate_VCF_Folder(path_in_folder)
    if valid == 1:
        PRINT.printE(STR__IO_error_read_folder)
        PRINT.printE(STR__use_help)
        return 1
    elif valid == 2:
        PRINT.printE(STR__error_no_VCF.format(f = path_in_folder))
        PRINT.printE(STR__use_help)
        return 1
    #
    path_in_chromosomes = inputs.pop(0)
    valid = Validate_Read_Path(path_in_chromosomes)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in_chromosomes))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out = ""
    classification = DEFAULT__Classification
    values = DEFAULT__Values
    tiebreaker = DEFAULT__Tiebreaker
    
    # Initial validation
    while inputs:
        arg = inputs.pop(0)
        try: # Following arguments
            if arg in ["-o", "-c", "-v", "-t"]:
                arg2 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            print arg
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        # Flag-dependent response
        if arg == "-o": path_out = arg2
        elif arg == "-c":
            if arg2 in LIST__PF:
                classification = CLASSIFICATION.PASS_FAIL
            elif arg2 in LIST__PSF:
                classification = CLASSIFICATION.PASS_SKIP_FAIL
            elif arg2 in LIST__PFF:
                classification = CLASSIFICATION.PERFECT_FAM_FAIL
            elif arg2 in LIST__PFSF:
                classification = CLASSIFICATION.PERFECT_FAM_SKIP_FAIL
            else:
                PRINT.printE(STR__invalid_classification.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        elif arg == "-v":
            if arg2 in LIST__FLAGS: values = VALUES.FLAGS
            elif arg2 in LIST__COUNTS: values = VALUES.COUNT
            elif arg2 in LIST__P_FLAGS: values = VALUES.FLAGS_POSITIVE
            else:
                PRINT.printE(STR__invalid_values.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        else: #arg == "-t"
            if arg2 in LIST__FULL: tiebreaker = TIEBREAKER.FULL
            elif arg2 in LIST__PARTIAL: tiebreaker = TIEBREAKER.PARTIAL
            elif arg2 in LIST__NONE: tiebreaker = TIEBREAKER.NONE
            else:
                PRINT.printE(STR__invalid_tiebreaker.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
    
    # Automated output path generation
    if not path_out: path_out = path_in_folder + FILEMOD
    
    # Validate output path
    valid_out = Validate_Write_Path__FILE(path_out)
    if valid_out == 2: return 0
    if valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    
    # Run program
    exit_state = Compare_RetroSeq_Results__STR(path_in_baseline, path_in_folder,
            path_in_chromosomes, path_out, classification, values, tiebreaker)
    
    # Exit
    if exit_state == 0: return 0
    else:
        PRINT.printE(STR__use_help)
        return 1



def Validate_VCF_Folder(dirpath):
    """
    Validates the dirpath of the input folder as containing VCF files.
    Return 0 if the dirpath is valid and contains at least 1 VCF file.
    Return 1 if the dirpath is valid but contains no VCF files.
    Return 2 if the dirpath is invalid.
    
    Validate_Read_Path(str) -> int
    """
    try:
        os.listdir(dirpath)
        files = Get_Files_W_Extensions(dirpath, LIST__VCF)
        if len(files) > 0: return 0
        return 1
    except:
        return 2

def Validate_Write_Path__FILE(filepath):
    """
    Validates the filepath of the output file.
    Return 0 if the filepath is writtable.
    Return 1 if the user decides to overwrite an existing file.
    Return 2 if the user declines to overwrite an existing file.
    Return 3 if the file exists and the program is set to forbid overwriting.
    Return 4 if the program is unable to write to the filepath specified.
    
    Validate_Write_Path(str) -> int
    """
    try:
        f = open(filepath, "U")
        f.close()
    except: # File does not exist. 
        try:
            f = open(filepath, "w")
            f.close()
            return 0 # File does not exist and it is possible to write
        except:
            return 4 # File does not exist but it is not possible to write
    # File exists
    if WRITE_PREVENT: return 3
    if WRITE_CONFIRM:
        confirm = raw_input(STR__overwrite_confirm.format(f = filepath))
        if confirm not in LIST__yes: return 2
    # User is not prevented from overwritting and may have chosen to overwrite
    try:
        f = open(filepath, "w")
        f.close()
        if WRITE_CONFIRM: return 1 # User has chosen to overwrite existing file
        return 0 # Overwriting existing file is possible
    except:
        return 4 # Unable to write to specified filepath

def Validate_Write_Path__FOLDER(folder_path):
    """
    Validates the writepath of the output folder.
    Attempts to create the folder if it does not exist.
    
    Return 0 if the folder path is valid and empty and can be written into.
    Return 1 if the folder path is valid and the user decides to overwrite
            existing files.
    Return 2 if the folder path is valid and empty but cannot be written into.
    Return 3 if the folder path is valid and the user declines to overwrite
            existing files.
    Return 4 if the folder path is valid, but contains existing files and the
            program is set to forbid overwriting.
    Return 5 is the folder path does not exist and cannot be created.
    Return 6 for unexpected errors.
    
    Validate_Folder_Path(str, str) -> int
    """
    new_dir = False
    # Create folder if it does not exist
    if not os.path.isdir(folder_path):
        try:
            os.mkdir(folder_path)
        except:
            return 5
        new_dir = True
    
    # Create a random file for testing purposes
    random_name = str(Random.random())
    random_path = folder_path + "\\" + random_name
    while os.path.exists(random_path):
        random_name = str(Random.random())
        random_path = folder_path + "\\" + random_name
    # Attempt to write to the folder
    try:
        f = open(random_path, "w")
    except:
        return 2
    # Unexpected errors
    try:
        f.close()
        os.remove(random_path)
    except:
        return 6
    
    if new_dir: return 0
    
    # OVERWRITE TESTING
    FLAG_exists = 0
    if WRITE_PREVENT: return 4
    if WRITE_CONFIRM:
        confirm = raw_input(STR__overwrite_confirm_2.format(f=folder_path))
        if confirm not in LIST__yes: return 3        
    return 1



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Compare_RetroSeq_Results(sys.argv)


