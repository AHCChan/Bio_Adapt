HELP_DOC = """
BED POSTMERGE UNCOLLAPSE COMBINE REPLICATES
(version 1.0)
by Angelo Chan

Designed to run on the output of BED__Postmerge_Uncollapse.py. (Bio_Files_Tools)

This is a program for summarizing the count and absence of genetic elements in
multiple samples, aggregated by experimental group, with different replicates
of the same experimental group being combined.



USAGE:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py <input_path>
            <grouping_file_path> [-m <output_path_main_p> <output_path_main_a>
            <output_path_main_t> <output_path_main_m>] [-a <output_path_all_p>
            <output_path_all_a> <output_path_all_t> <output_path_all_m>]
            [-u <output_path_uniques_p> <output_path_uniques_a>
            <output_path_uniques_t> <output_path_uniques_m>]



MANDATORY:
    
    input_path
        
        The filepath of the input BED file. This input file should be the
        _SUMMARY file produced by BED__Postmerge_Uncollapse.py.

        This summary file will be a TSV with the following columns:
            1)  The first 3 columns will be the genomic coordinates of the
                merged regions. (1 region per row)
            2)  For each file which went into the pipeline, the number of loci
                from that file which went into merged region.
    
    grouping_file_path
        
        The two-column TSV denoting which columns belong to which experimental
        groups. The first column will contain the name of the experimental
        group, while second column contains comma-separated column numbers.
        
        The column numbers follow a 1-index system, and ignore the first 3
        columns of the BED file. (which contain the genomic coordinates) Thus,
        the very first data column is column "1".
        
        Ex.:
            Wildtype_ctrl   1,2,3,4
            Wildtype_drug   5,6,7
            Mutated_ctrl    8,9,10
            Mutated_drug    11,12,13,14

OPTIONAL:[-m <output_path_main>]
            [-u <output_path_uniques>] [-s <output_path_summary>]
    
    output_path_main_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output files which contain the results of combining
        all replicates from each experimental group.
    
    output_path_all_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output files which contain the genetic elements
        found in all experimental groups.
    
    output_path_uniques_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output folders which contain the group-specific
        files, for genetic elements found exclusively in one experimental
        group and none of the others. (None of the columns specified in the
        grouping file)
    
    (SUFFIXES)
        
        (present), (all), (total), (mean)
            
            [*_p] contains either a 1 or a 0. 1 if at least 1 replicate contains
                    at least 1 count. 0 otherwise.
            [*_a] contains either a 1 or a 0. 1 if all replicates contain at
                    least 1 count. 0 otherwise.
            [*_t] contains the total counts across all replicates.
            [*_m] contains the average number of counts per replicate.



EXAMPLES:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py
            data\summary.bed data\list_of_groups.tsv
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py
            data\summary.bed data\list_of_groups.tsv -m present.bed
            present_all_reps.bed main_total_counts.bed main_average_counts.bed
            -a present_all_groups.bed present_all_groups_all_reps.bed
            present_all_totals.bed present_all_mean.bed -u uniquely_present
            uniquely_present_all_reps uniquely_present_totals
            unique_present_averages

USAGE:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py <input_path>
            <grouping_file_path> [-m <output_path_main_p> <output_path_main_a>
            <output_path_main_t> <output_path_main_m>] [-a <output_path_all_p>
            <output_path_all_a> <output_path_all_t> <output_path_all_m>]
            [-u <output_path_uniques_p> <output_path_uniques_a>
            <output_path_uniques_t> <output_path_uniques_m>]
"""

NAME = "BED__Postmerge_Uncollapse.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD__MAIN = "__MAIN"
FILEMOD__ALL = "__ALLGROUPS"
FILEMOD__UNIQUE = "__UNIQUE"

FILEMOD__P = "_PRESENT"
FILEMOD__A = "_ALLREPS"
FILEMOD__T = "_TOTALS"
FILEMOD__M = "_AVERAGES"



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"BED__Postmerge_Uncollapse_Combine_Replicates.py -h"



STR__invalid_column_number = """
ERROR: Invalid column number: {s}
Please specify positive integers which correspond to column numbers which exist
in the main input file."""



STR__metrics = """
                       Rows in input file: {A}

    Universally present in all groups:

                In at least one replicate: {B}
                        In all replicates: {C}

    Uniquely present in one group:
        
                In at least one replicate: {D}
                        In all replicates: {E}
"""



STR__combine_begin = "\nRunning Combine_Replicates..."

STR__combine_complete = "\nCombine_Replicates successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Combine_Replicates(path_in, path_groups, paths_out):
    """
    Generate a series of FASTA files each containing a synthetic chromosome.
    
    @path_in
            (str - filepath)
            The filepath of the input BED file produced by
            BED__Postmerge_Uncollapse.py.
    @path_groups
            (str - filepath)
            The filepath of the TSV file specifying which columns belong to
            which experimental groups.
    @paths_out
            (list<list<str>> - 8x filepaths, 4x dirpaths)
            Nested lists, containing the output filepaths and dirpaths.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    Combine_Replicates(str, str, [[str, str, str, str], [str, str, str, str],
            [str, str, str, str]]) -> int
    """
    PRINT.printP(STR__combine_begin)
    
    # Setup reporting
    count_rows = 0
    count_universal_one = 0
    count_universal_all = 0
    count_unique_one = 0
    count_unique_all = 0
    
    # Get groups
    groups_dict = Get_Groups_Dict(path_groups)
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_in)
    f.Set_Delimiter("\t")
    f.Open()
    f.Close()
    outputs = Setup_Outputs(paths_out, groups_dict)
    
    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        count_rows += 1
        # Main loop logic
    
    # Finish
    f.Close()
    PRINT.printP(STR__combine_complete)
    
    # Reporting
    Report_Metrics([count_rows, count_universal_one, count_universal_all,
            count_unique_one, count_unique_all])
    
    # Wrap up
    return 0



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Combine_Uncollapsed(raw_command_line_input):
    """
    Parse the command line input and call the Generate_Synthetic_Genome function
    with appropriate arguments if the command line input is valid.
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
    if len(inputs) < 2:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_in = inputs.pop(0)
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    path_groups = inputs.pop(0)
    valid = Validate_Read_Path(path_groups)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_groups))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-m", "-a", "-u"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
                arg4 = inputs.pop(0)
                arg5 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        if arg == "-m":
            pass
        if arg == "-a":
            pass
        else: # arg == "-u"
            path_out_s = arg2
    
    # Validate output paths
    
    
    # Run program
    exit_state = Uncollapse_BED(path_in, path_groups, [])
    
    # Exit
    if exit_state == 0: return 0
    else: return 1


    
def Validate_Write_Path(filepath):
    """
    Validates the filepath of the input file.
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
        confirm = raw_input(STR__overwrite_confirm.format(f=filepath))
        if confirm not in LIST__yes: return 2
    # User is not prevented from overwritting and may have chosen to overwrite
    try:
        f = open(filepath, "w")
        f.close()
        if WRITE_CONFIRM: return 1 # User has chosen to overwrite existing file
        return 0 # Overwriting existing file is possible
    except:
        return 4 # Unable to write to specified filepath



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Combine_Uncollapsed(sys.argv)
