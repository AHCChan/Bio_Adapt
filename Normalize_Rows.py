HELP_DOC = """
ROW NORMALIZER
(version 1.0)
by Angelo Chan

This program takes a table file and normalizes each row so that the rows are
comparable.

Several strategies can be employed for this. A common strategy is to shift all
values up or down so that the average value of each row is 0, also known as
mean-shifting to 0. (Certain types of data may need to be log-transformed first
before this approach is appropriate) Median-shifting, is where all data values
are adjusted so that the median value of every row is the same. Quartile
shifting is when all values are shifted such that the 1st and 3rd quartile
values (the 3rd quartile value is the value greater than 75% of the values in
the data, while the 1st quartile value is the value greater than 25% of the
values in the data) are the same across all datasets.



USAGE:
    
    python27 Normalize_Rows.py <input_file> [-o <output_file> <output_format>]
            [-h <header>] [-i <id>] [-n MEAN|MEDIAN|QUARTILE <m1> [<m2>]]



MANDATORY:
    
    input_file
        
        The filepath of the input file to be normalized.
    
    input_format
        
        The file format of the input file. Acceptable options are:
            tsv - Tab-separated values
            csv - Comma-separated values
            ssv - Space-separated values



OPTIONAL:
    
    output_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file.

    output_format
        
        (DEFAULTS to the same format as the input file)
        
        The file format of the input file. Acceptable options are:
            tsv - Tab-separated values
            csv - Comma-separated values
            ssv - Space-separated values

    header
        
        (DEFAULT: Y)
        
        Whether or not the first row of the file is a header row and should be
        excluded from normalization.

    id
        
        (DEFAULT: Y)
        
        Whether or not the first column of the file contains ID numbers/words
        and should be omitted from the normalization process.
    
    MEAN|MEDIAN|QUARTILE
        
        (DEFAULT: MEAN)
        
        What kind of normalization to perform:
            
            MEAN:
                All rows will have the same mean.
            MEDIAN:
                All rows will have the same median.
            QUARTILE:
                All rows will have the same values for their 1st and 3rd
                quartiles.
    
    m1
        
        (DEFAULT: 0)
        
        The number which the mean, median, or 1st quartile will normalize to.

    m2
        
        (DEFAULT: none)
        
        The value which the 3rd quartile will normalize to. Only applicable when
        the QUARTILE normalization method has been specified.



EXAMPLES SCENARIO EXPLANATION:
    
    1:
    A minimal use case.
    
    2:
    Output file specified. Change the format from what was originally a CSV file
    to a TSV file.
    
    3:
    Dealing with a data table which contains no headers or row IDs.
    
    4:
    Change the median value of each row to 1000.
    
    5:
    Alter the data so that for each row, the first quartile is -10 and the third
    quartile is 10.

EXAMPLE:
    
    python27 Normalize_Rows.py Path/data.tab tsv
    
    python27 Normalize_Rows.py Path/data.csv csv -o Path/results.tsv tsv
    
    python27 Normalize_Rows.py Path/raw_data.tab tsv -h N -i N
    
    python27 Normalize_Rows.py Path/data.tab tsv -n MEDIAN 1000
    
    python27 Normalize_Rows.py Path/data.tab tsv -n QUARTILE -10 10

USAGE:
    
    python27 Normalize_Rows.py <input_file> [-o <output_file> <output_format>]
            [-h <header>] [-i <id>] [-n MEAN|MEDIAN|QUARTILE <m1> [<m2>]]
"""

NAME = "Normalize_Rows.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__NORMALIZED"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__headers = True
DEFAULT__IDs = True
DEFAULT__normalize = 1 # MEAN
DEFAULT__M1 = 0
DEFAULT__M2 = 0



# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Enums ########################################################################

class NORMALIZE:
    MEAN=1
    MEDIAN=2
    QUARTILE=3



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\tpython "\
"Row_Normalizer.py -h"



STR__invalid_normalization = """
ERROR: Invalid normalization method specified: {s}
Please specify one of the following:
    Mean
    Median
    Quartile"""

STR__specify_number = """
ERROR: Invalid parameter for normalization: {s}
Please specify a number."""

STR__need_two_for_quartile = "\nERROR: You need to specify two numbers when "\
"normalizing by quartile."



STR__metrics = """
       Total Rows: {A}
    Total Columns: {B}
    
            {C} {G} {K} {O}
     Lowest {D} {H} {L} {P}
    Average {E} {I} {M} {Q}
    Highest {F} {J} {N} {R}
"""
STR__metrics_h1 = "Mean"
STR__metrics_h1 = "1st Q"
STR__metrics_h1 = "Median"
STR__metrics_h1 = "3rd Q"

STR__report_begin = "\nRunning Normalize_Rows..."

STR__report_complete = "\nNormalize_Rows successfully finished."



# Lists ########################################################################

LIST__mean = ["MEAN", "Mean", "mean", "AVERAGE", "Average", "average"]
LIST__median = ["MED", "Med", "med", "MEDIAN", "Median", "median"]
LIST__quartile = ["Q", "q", "QUART", "Quart", "quart", "QUARTILE", "Quartile",
        "quartile"]



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Normalize_Rows(path_in, delim_in, path_out, delim_out, headers, IDs,
            normalize, normalize_params):
    """
    """
    PRINT.printP(STR__report_begin)
    
    # 
    
    PRINT.printP(STR__report_complete)
    
    # Reporting
    Report_Metrics(metrics)
    
    # Wrap up
    return 0



def Report_Metrics(metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The number of rows of data
                * The number of columns of data
                * The highest, average, and lowest means
                * The highest, average, and lowest 1st quartile values
                * The highest, average, and lowest medians
                * The highest, average, and lowest 3rd quartile values
    
    Report_Metrics(list<int>(14)) -> None
    """
    # Unpacking
    # Print
    pass



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Normalize_Rows(raw_command_line_input):
    """
    Parse the command line input and call the Normalize_Rows function with
    appropriate arguments if the command line input is valid.
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
    
    # Minimum inputs validation
    if len(inputs) < 2:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_in = inputs.pop(0) # Input file
    valid = Validate_FASTA_Folder(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    delim_in_str = inputs.pop(0) # Input delim
    input_delim = Validate_Table_Type(input_delim_str)
    if not input_delim:
        PRINT.printE(STR__invalid_table_format.format(s = input_delim_str))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    headers = DEFAULT__headers
    IDs = DEFAULT__IDs
    normalize = DEFAULT__normalize
    M1 = DEFAULT__M1
    M2 = DEFAULT__M2
    path_out = ""
    delim_out = delim_in
    
    # Initial validation
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-h", "-i", "-n"]:
                arg2 = inputs.pop(0)
            elif arg in ["-o"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        if arg in ["-h", "-i"]:
            boolean = Validate_Bool(arg2)
            if boolean == None:
                PRINT.printE(STR__invalid_bool.format(s = arg2))
                return 1
            if arg == "-h":
                headers = boolean
            else: # "-i"
                IDs = boolean
        elif arg == "-o":
            path_out = arg2
            delim_out = Validate_Table_Type(arg3)
            if not delim_out:
                PRINT.printE(STR__invalid_table_format.format(s = arg3))
                PRINT.printE(STR__use_help)
                return 1
        elif arg == "-n":
            normalize = Validate_Normalization(arg2)
            if not normalize:
                PRINT.printE(STR__invalid_normalization.format(s = arg2))
            M1 = Validate_Number(arg2)
            if M1 == None:
                PRINT.printE(STR__specify_number.format(s = arg2))
                return 1
            if normalize == NORMALIZE.QUARTILE:
                try:
                    arg3 = inputs.pop(0)
                except:
                    PRINT.printE(STR__need_two_for_quartile)
                M2 = Validate_Number(arg3)
                if M2 == None:
                    PRINT.printE(STR__specify_number.format(s = arg3))
                    return 1
    
    # Automated output path generation
    if not path_put:
        path_out = Generate_Default_Output_Filepath_Norm_Rows(path_in, delim_in)
    
    # Validate output path
    valid_out = Validate_Write_Path(path_out)
    if valid_out == 2: return 0
    if valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    
    # Run program
    Normalize_Rows(path_in, delim_in, path_out, delim_out, headers, IDs,
            normalize, [M1, M2])
    
    # Safe exit
    if exit_state == 0: return 0
    else:
        return 1



def Generate_Default_Output_Filepath_Norm_Rows(filepath, delim):
    """
    Generate output filepath based on the input filepath and delimiter.
    
    Generate_Default_Output_Filepath_Norm_Rows(str, str) -> str
    """
    abs_path = os.path.abspath(filepath)
    dir_path = os.path.dirname(dir_path)
    file_name = Get_File_Name(filepath)
    ext = DICT__table_delim_to_ext(delim)
    return dir_path + "\\" + file_name + "." + ext

def Validate_Normalization(string):
    """
    Validates the normalization method specified.
    Return 0 if it is invalid. Otherwise return the appropriate ENUM:
        1 - Mean-shift
        2 - Median-shift
        3 - Quartile-shift
    
    Validate_Normalization(str) -> int
    """
    if string in LIST__mean: return NORMALIZE.MEAN
    if string in LIST__median: return NORMALIZE.MEDIAN
    if string in LIST__quartile: return NORMALIZE.QUARTILE
    return 0



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



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Normalize_Rows(sys.argv)
