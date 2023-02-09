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
    
    python27 Normalize_Rows.py <input_file> <input_format> [-o <output_file>
            <output_format>] [-h <header>] [-i <id>] [-n MEAN|MEDIAN|QUARTILE
            <m1> [<m2>]]



MANDATORY:
    
    input_file
        
        The filepath of the input file to be normalized.
        
        NOTE: A table file with column headers is acceptable, a table file with
        header comments is not.
    
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
    
    python27 Normalize_Rows.py <input_file> <input_format> [-o <output_file>
            <output_format>] [-h <header>] [-i <id>] [-n MEAN|MEDIAN|QUARTILE
            <m1> [<m2>]]
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



STR__invalid_value = "\nERROR: Invalid value: {s}"



STR__empty_table_file = "\nERROR: Table file is empty: {s}"



STR__metrics = """
       Total Rows: {A}
    Total Columns: {B}
    
             {C}  {G}  {K}  {O}
     Lowest  {D}  {H}  {L}  {P}
    Average  {E}  {I}  {M}  {Q}
    Highest  {F}  {J}  {N}  {R}
"""
STR__metrics_h1 = "Mean"
STR__metrics_h2 = "1st Q"
STR__metrics_h3 = "Median"
STR__metrics_h4 = "3rd Q"

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
    Normalize the data so that every row has the same mean, median, or quartile
    values.
    
    @path_in
            (str - filepath)
            The filepath of the data table to be normalized.
            A table file with column headers is acceptable, a table file with
            header comments is not.
    @delim_in
            (str)
            The delimiter to be used for the input table file. File formats and
            their corresponding delimiters are as follows:
                TSV - "\t" (tab character)
                CSV - ","  (comma character)
                SSV - " "  (whitespace character)
    @path_out
            (str - filepath)
            The filepath of the file where the output will be written into.
    @delim_out
            (str)
            The delimiter to be used for the output table file. File formats and
            their corresponding delimiters are as follows:
                TSV - "\t" (tab character)
                CSV - ","  (comma character)
                SSV - " "  (whitespace character)
    @headers
            (bool)
            Whether or not there are column headers in the input file.
            If there are, they will be excluded from the normalization process
            and written into the output file as is.
    @IDs
            (bool)
            Whether or not the first column contains data entry IDs.
            If there are, they will be excluded from the normalization process
            and written into the output file as is.
    @normalize
            (int - ENUM)
            An integer denoting what kind of normalization should be performed.
                1 - Mean-shifting
                2 - Median-shifting
                3 - Quartile-shifting
    @normalize_params
            (list<float>)
            A list of parameters to be used for the normalization process. When
            mean-shifting or median-shifting, the first number in this list will
            be what the mean/median of the data will be shifted to. When
            quartile shifting, the first and second numbers in this list will be
            what the 1st and 3rd quartile of each row of data will be shifted
            to.
    
    Normalize_Rows(str, str, str, str, bool, bool, int, list<float>) -> int
    """
    PRINT.printP(STR__report_begin)
    
    # Setup reporting
    total_lines = 0
    total_columns = 0
    means = []
    q1s = []
    medians = []
    q3s = []
    
    # Setup file I/O
    o = open(path_out, "w")
    f = Table_Reader(path_in)
    f.Set_Delimiter(delim_in)
    f.Open()
    
    # Column count and indexes
    f.Read()
    values = f.Get()
    total_columns = len(values)
    if IDs: total_columns = total_columns - 1
    range_ = range(total_columns)
    indexes = Get_Indexes(total_columns)
    f.Close()
    f.Open()
    
    # Headers
    if headers:
        f.Read()
        values = f.Get()
        sb = delim_out.join(values) + "\n"
        o.write(sb)
    
    # Main loop
    while not f.End():
        total_lines += 1
        # Read
        f.Read()
        values = f.Get()
        # ID
        if IDs:
            ID = values[0]
            o.write(ID + delim_out)
            values = values[1:]
        # Convert to floats
        for i in range_:
            try:
                value = float(values[i])
            except:
                PRINT.printE(STR__invalid_value.format(s=values[i]))
                return 1
            values[i] = value
        # Metrics
        mean = sum(values)/total_columns
        q1, median, q3 = Get_Quartiles_And_Medians(values, indexes)
        means.append(mean)
        q1s.append(q1)
        medians.append(median)
        q3s.append(q3)
        # Normalize
        if normalize == NORMALIZE.MEAN:
            difference = normalize_params[0] - mean
            for i in range_:
                values[i] = str(values[i] + difference)
        elif normalize == NORMALIZE.MEDIAN:
            difference = normalize_params[0] - median
            for i in range_:
                values[i] = str(values[i] + difference)
        elif normalize == NORMALIZE.QUARTILE:
            gradient_goal = normalize_params[1] - normalize_params[0]
            gradient_real = q3-q1
            ratio = gradient_goal/gradient_real
            offset = normalize_params[0] - ratio * q1
            for i in range_:
                values[i] = str(ratio * values[i] + offset)
        # Write
        sb = delim_out.join(values) + "\n"
        o.write(sb)
    # Close up
    f.Close()
    o.close()
    PRINT.printP(STR__report_complete)
    
    # Empty file
    if total_lines == 0:
        PRINT.printE(STR__empty_table_file.format(s=path_in))
        return 1
    
    # Calculate
    table_metrics = Get_Metrics([means, q1s, medians, q3s])
    
    # Reporting
    Report_Metrics([total_lines, total_columns] + table_metrics)
    
    # Wrap up
    return 0



def Get_Indexes(length):
    """
    Return the indexes needed to get the 1st quartile value, the median value,
    and the 3rd quartile value.
    
    @length
            (int)
            The number of elements in the list.
    
    Get_Indexes(int) -> [int, int, int, int, int, int]
    """
    max_index = float(length-1)
    # Median
    median_1 = (length-1)/2
    median_2 = length/2
    # q1
    q1_index_float = max_index/4
    q1_index_int = int(q1_index_float)
    if q1_index_float == q1_index_int:
        q1_index_1 = q1_index_int
        q1_index_2 = q1_index_int
    else:
        q1_index_1 = q1_index_int
        q1_index_2 = q1_index_int+1
    # q3
    q3_index_float = (max_index*3)/4
    q3_index_int = int(q3_index_float)
    if q3_index_float == q3_index_int:
        q3_index_1 = q3_index_int
        q3_index_2 = q3_index_int
    else:
        q3_index_1 = q3_index_int
        q3_index_2 = q3_index_int+1
    # Return
    return [q1_index_1, q1_index_2, median_1, median_2, q3_index_1, q3_index_2]



def Get_Quartiles_And_Medians(values, indexes):
    """
    Return the 1st quartile value, the median, and the 3rd quartile value of a
    list of numbers.
    
    @values
            (list<float>)
            The numbers being analyzed.
    @indexes
            (list<int>(6))
            A precalculated set of indexes. The indexes are used to access the
            elements in a sorted list of values to calculate the 1st quartile
            value, the median, and the 3rd quartile value of the list of
            numbers.
    
    Get_Quartiles_And_Medians(list<float>, list<int>(6)) -> list<float>(3)
    """
    # Sort
    values = sorted(values)
    # q1
    q1_1 = values[indexes[0]]
    q1_2 = values[indexes[1]]
    q1 = (q1_1+q1_2)/2
    # Median
    median_1 = values[indexes[2]]
    median_2 = values[indexes[3]]
    median = (median_1+median_2)/2
    # q3
    q3_1 = values[indexes[4]]
    q3_2 = values[indexes[5]]
    q3 = (q3_1+q3_2)/2
    # Return
    return [q1, median, q3]



def Get_Metrics(lists):
    """
    Return the lowest, average, and highest value in every list in [lists].
    
    @lists
            (list<list<float>>(4))
            The lists of values to get the lowest, average, and highest values
            for.
    
    Get_Metrics(list<list<float>>(4)) -> list<float>(12)
    """
    length = len(lists[0])
    results = []
    for list_ in lists:
        # Setup
        lowest = list_[0]
        highest = list_[0]
        total = 0
        # Iterate
        for value in list_:
            if value < lowest: lowest = value
            if value > highest: highest = value
            total += value
        average = total/length
        # Append
        results.append(lowest)
        results.append(average)
        results.append(highest)
    return results



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
    # Setup
    length = len(metrics)
    range_ = range(length)
    # 
    for i in range_:
        value = metrics[i]
        string = str(value)
        if i > 1:
            string = Trim_Percentage_Str(string, 2)
        metrics[i] = string
    # Repacking
    row_and_column = metrics[:2]
    col_1 = [STR__metrics_h1] + metrics[2:5]
    col_2 = [STR__metrics_h2] + metrics[5:8]
    col_3 = [STR__metrics_h3] + metrics[8:11]
    col_4 = [STR__metrics_h4] + metrics[11:]
    # Pad
    row_and_column = Pad_Column(row_and_column, 0, 0, " ", 0)
    col_1 = Pad_Column(col_1, 0, 0, " ", 0)
    col_2 = Pad_Column(col_2, 0, 0, " ", 0)
    col_3 = Pad_Column(col_3, 0, 0, " ", 0)
    col_4 = Pad_Column(col_4, 0, 0, " ", 0)
    # Strings
    PRINT.printM(STR__metrics.format(
            A = row_and_column[0], B = row_and_column[1],
            C = col_1[0], D = col_1[1], E = col_1[2], F = col_1[3],
            G = col_2[0], H = col_2[1], I = col_2[2], J = col_2[3],
            K = col_3[0], L = col_3[1], M = col_3[2], N = col_3[3],
            O = col_4[0], P = col_4[1], Q = col_4[2], R = col_4[3]))  



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
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    delim_in_str = inputs.pop(0) # Input delim
    delim_in = Validate_Table_Type(delim_in_str)
    if not delim_in:
        PRINT.printE(STR__invalid_table_format.format(s = delim_in_str))
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
            if arg in ["-h", "-i"]:
                arg2 = inputs.pop(0)
            elif arg in ["-o", "-n"]:
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
            M1 = Validate_Number(arg3)
            if M1 == None:
                PRINT.printE(STR__specify_number.format(s = arg3))
                return 1
            if normalize == NORMALIZE.QUARTILE:
                try:
                    arg4 = inputs.pop(0)
                except:
                    PRINT.printE(STR__need_two_for_quartile)
                M2 = Validate_Number(arg4)
                if M2 == None:
                    PRINT.printE(STR__specify_number.format(s = arg4))
                    return 1
    
    # Automated output path generation
    if not path_out:
        path_out = Generate_Default_Output_Filepath_Norm_Rows(path_in, delim_in)
    
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
    exit_state = Normalize_Rows(path_in, delim_in, path_out, delim_out, headers,
            IDs, normalize, [M1, M2])
    
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
    dir_path = os.path.dirname(abs_path)
    file_name = Get_File_Name(filepath)
    ext = DICT__table_delim_to_ext[delim]
    return dir_path + "\\" + file_name + FILEMOD + "." + ext

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
