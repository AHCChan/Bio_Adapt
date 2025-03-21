HELP_DOC = """
SECONDARY ALIGNMENT MASKER
(version 1.1)
by Angelo Chan

Modify a SAM (sequence alignment map) file to hide secondary alignments.

The SAM file is presumed to contain secondary alignment info using the "XS" 
flag in the columns after the 11th. It is also presumed to have been 
generated such that the XS value will be strictly non-positive.



USAGE:
    
    python27 Mask_Secondary_Alignments <input_SAM> [-o <output_path>]
            [-t <threshold>]



MANDATORY:
    
    input_SAM
        
        The filepath of the input SAM file. No headers allowed.

OPTIONAL:
    
    output_path
        
        (DEFAULT path generation available)
        
        The filepath of the output file.
    
    threshold
        
        (DEFAULT: 0)
        
        The threshold at which the secondary alignments will be retained.
            * At 1, no secondary alignments will be retained.
            * At 0, only perfect secondary alignments will be retained.
            * For any negative integer, only secondary alignments with that
                    relative score or higher will be retained.



EXAMPLES:
    
    python27 Mask_Secondary_Alignments path\data.SAM

    python27 Mask_Secondary_Alignments path\data.SAM path\results.SAM -t -2

USAGE:
    
    python27 Mask_Secondary_Alignments <input_SAM> [-o <output_path>]
            [-t <threshold>]
"""

NAME = "Mask_Secondary_Alignments.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__NO_SECONDARY"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__threshold = 0



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import * # 1.9

from Table_File_Reader import * # 1.1



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"Mask_Secondary_Alignments.py -h"



STR__invalid_threshold = """
ERROR: Invalid threshold: {s}
Please specify either 1, 0, or a negative integer.
"""



STR__metrics = """
                      Reads in File: {A}

    Reads with Secondary Alignments: {B}

      Secondary Alignments Retained: {C} ({D}%)
       Secondary Alignments Removed: {E} ({F}%)
"""



STR__masking_begin = "\nRunning Mask_Secondary_Alignments..."

STR__masking_complete = "\nMask_Secondary_Alignments successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Mask_Secondary_Alignments(path_SAM, path_out, threshold):
    """
    Modify a SAM (sequence alignment map) file to hide secondary alignments.

    The SAM file is presumed to contain secondary alignment info using the "XS" 
    flag in the columns after the 11th. It is also presumed to have been 
    generated such that the XS value will be strictly non-positive.
    
    @path_SAM
            (str - filepath)
            The filepath of the input SAM file.
    @path_out
            (str - filepath)
            The filepath of the output file.
    @threshold
            (int)
            The threshold at which the secondary alignments will be retained.
                * At 1, no secondary alignments will be retained.
                * At 0, only perfect secondary alignments will be retained.
                * For any negative integer, only secondary alignments with that
                        relative score or higher will be retained.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    BED_Normalize_Chr(str, str, int) -> int
    """
    # Setup reporting
    row_count = 0
    secondary_count = 0
    retention_count = 0
    removal_count = 0
    
    PRINT.printP(STR__masking_begin)
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_SAM)
    f.Set_Delimiter("\t")
    f.Open()
    f.Close()
    o = open(path_out, "w")
    
    # Header scouting
    count = 0
    flag = True
    s = open(path_SAM, "U")
    while flag:
        line = s.readline()
        if line[0] == "@":
            count += 1
            o.write(line)
        else:
            flag = False
    s.close()
    
    # Header skipping
    f.Open()
    while count > 0:
        count -= 1
        f.Read()
    
    # Main loop
    while not f.EOF:
        f.Read()
        row_count += 1
        # Read and parse
        values = f.Get()
        non_flags = values[:12]
        flags = values[12:]
        # Check
        go_ahead = False
        for i in flags:
            if "XS:i:" in i: go_ahead = True
        if go_ahead:
            secondary_count += 1
            temp_flags = []
            for i in flags:
                if "XS:i:" in i:
                    if threshold != 1:
                        XS_str = i.split("XS:i:")[1]
                        XS = int(XS_str)
                        if XS >= threshold:
                            retention_count += 1
                            temp_flags.append(i)
                        else:
                            removal_count += 1
                    else:
                        removal_count += 1
                else:
                    temp_flags.append(i)
            # Combine and Write
            new_values = non_flags + temp_flags
            sb = "\t".join(new_values) + "\n"
            o.write(sb)
        else:
            raw = f.Get_Raw()
            o.write(raw)
    
    # Finish
    f.Close()
    o.close()
    PRINT.printP(STR__masking_complete)
    
    # Reporting
    Report_Metrics([row_count, secondary_count, retention_count, removal_count])
    
    # Wrap up
    return 0



def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The total number of rows
                * The total number of rows with secondary alignments
                * The total number of secondary alignments retained
                * The total number of secondary alignments removed
        
    Report_Metrics([int, int, int, int]) -> None
    """
    # Unpacking
    total, secondary, retained, removed = summary_metrics
    # Calculations
    if secondary:
        percentage_retained = (retained*100.00)/secondary
        percentage_removed = (removed*100.00)/secondary
    else:
        percentage_retained = 0.0
        percentage_removed = 0.0
    # Strings
    total = str(total)
    secondary = str(secondary)
    retained = str(retained)
    removed = str(removed)
    percentage_retained = str(percentage_retained)
    percentage_removed = str(percentage_removed)
    # Pad Column 1
    total = total
    secondary = secondary
    retained = retained
    removed = removed
    #
    col_1 = [total, secondary, retained, removed]
    col_1 = Pad_Column(col_1, 0, 0, " ", 0)
    total, secondary, retained, removed = col_1
    # Pad Column 2
    percentage_retained = Trim_Percentage_Str(percentage_retained, 2)
    percentage_removed = Trim_Percentage_Str(percentage_removed, 2)
    #
    col_2 = [percentage_retained, percentage_removed]
    col_2 = Pad_Column(col_2, 0, 0, " ", 0)
    percentage_retained, percentage_removed = col_2
    # Print
    PRINT.printM(STR__metrics.format(A = total, B = secondary,
            C = retained, D = percentage_retained,
            E = removed, F = percentage_removed))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Mask_Secondary_Alignments(raw_command_line_input):
    """
    Parse the command line input and call the Mask_Secondary_Alignments function
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
    if len(inputs) < 1:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_SAM = inputs.pop(0)
    valid = Validate_Read_Path(path_SAM)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_SAM))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out = Generate_Default_Output_File_Path_From_File(path_SAM, FILEMOD,
            True)
    threshold = DEFAULT__threshold
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o", "-t"]:
                arg2 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        if arg == "-o":
            path_out = arg2
        else: # arg == "-t"
            threshold = Validate_Int_Max(arg2, 1)
            if valid == None:
                PRINT.printE(STR__invalid_threshold.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
    
    # Validate output paths
    valid_out = Validate_Write_Path(path_out)
    if valid_out == 2: return 0
    if valid_out == 3:
        printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        printE(STR__In_error_write_unable)
        return 1
    
    # Run program
    exit_state = Mask_Secondary_Alignments(path_SAM, path_out, threshold)
    
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
    exit_code = Parse_Command_Line_Input__Mask_Secondary_Alignments(sys.argv)


