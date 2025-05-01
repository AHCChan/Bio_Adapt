HELP_DOC = """
SAM PAIRS TO XBED
(version 1.0)
by Angelo Chan

This is a program which takes an unsorted SAM file of aligned read pairs and
outputs the alignment data of the pairs. For each read pair, the following data
is outputted:
    
    1)  Read ID
    2*) Consensus chromosome name.
    3*) Fragment starting position.
    4*) Fragment ending position.
    5`) Distance between the two reads.
    6*) Fragment size.
    7)  Chromosome name of the first read's alignment.
    8)  Starting position of the first read's alignment.
    9)  Starting position of the first read's alignment.
    10')Chromosome name of the second read's alignment.
    11')Starting position of the second read's alignment.
    12')Starting position of the second read's alignment.
    
    * - Only applies if both reads aligned to the same chromosome. This column
        will have a placeholder if this is not the case.
    ` - Only applies if both reads aligned to the same chromosome and the two
        reads do not overlap.
    ' - Only applies if this is a read pair (based on read ID)

This data can then be refined and filtered using other tools, such as
Table2Table.py.

Placeholders:
    
    2)  .
    3)  -1
    4)  -1
    5)  -1
    6)  -1
    10) .
    11) -1
    12) -1



USAGE:
    
    python27 SAM_Pairs_To_xBED.py <SAM_file> [-o <output_file>]



MANDATORY:
    
    SAM_file
        
        (SAM file)
        The filepath of the input SAM file.

OPTIONAL:
    
    output_path
        
        (DEFAULT path generation available)
        
        The filepath of the output file.



EXAMPLES:
    
    python27 SAM_Pairs_To_xBED.py path\alignment.sam path\paired_pairs.tsv

USAGE:
    
    python27 SAM_Pairs_To_xBED.py <SAM_file> [-o <output_file>]
"""

NAME = "SAM_Pairs_To_xBED.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__PRE_BED"



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import * # 2.0

from Table_File_Reader import * #1.1



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"SAM_Pairs_To_xBED.py -h"



STR__metrics = """
                    Reads: {A}
                    Pairs: {B} ({C}%)
    Same-Chromosome Pairs: {D} ({E}% of Reads)
                           {F} ({G}% of Pairs)
         Average Gap Size: {H}
    Average Fragment Size: {I}
"""



STR__pair_begin = "\nRunning Pair_SAM_Reads..."

STR__pair_complete = "\nPair_SAM_Reads successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Pair_SAM_Reads(path_SAM, path_output):
    """
    Create a new set of values which would allow genomic coordinate data to be
    plotted linearly.
    
    @path_SAM
            (str - filepath)
            The filepath of the input SAM file.
    @path_output
            (str - filepath)
            The filepath of the output file.
    
    Return a value of 0 if the function runs successfully.
    Return a positive integer if there is a problem. The integer functions as an
    error code.
    
    Pair_SAM_Reads(str, str) -> int
    """
    # Setup reporting
    count_reads = 0
    count_pairs = 0
    count_pairs_same = 0
    total_gaps = 0
    total_frag_Ns = 0
    
    # Setup buffers
    stored = [] # [ID, chr, start, end]
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_SAM)
    f.Set_Delimiter("\t")
    f.Set_Header_Params(["@"]) # Optional
    o = open(path_output, "w")
    
    # Start and read header
    PRINT.printP(STR__pair_begin)
    f.Open()
    
    # Main loop
    while not f.EOF:
        f.Read()
        count_reads += 1
        ID = f[0]
        chrom = f[2]
        start = int(f[3])
        length = len(f[9])
        end = start + length - 1
        if (not stored):
            stored = [ID, chrom, start, end]
        elif ID != stored[0]:
            sb = (stored[0] + "\t.\t-1\t-1\t-1\t-1\t" + stored[1] + "\t" +
                    str(stored[2]) + "\t" + str(stored[3]) + "\t.\t-1\t-1\n")
            o.write(sb)
            stored = [ID, chrom, start, end]
        else: # ID == stored_ID
            count_pairs += 1
            if chrom != stored[1]: # Different chromosomes
                sb = (stored[0] + "\t.\t-1\t-1\t-1\t-1\t" + stored[1] + "\t" +
                        str(stored[2]) + "\t" + str(stored[3]) + "\t" + chrom +
                        "\t" + f[3] + "\t" + str(end) + "\n")
            else: # Same chromsomes
                count_pairs_same += 1
                gap = start - stored[3] - 1
                if gap < 0: gap = 0
                total_gaps += gap
                frag_size = (end - stored[2]) + 1
                if frag_size < 0: frag_size = 0
                total_frag_Ns += frag_size
                
                sb = (stored[0] + "\t" + chrom + "\t" + str(stored[2]) + "\t" +
                        str(end) + "\t" + str(gap) + "\t" + str(frag_size) +
                        "\t" + str(stored[2]) + "\t" + str(stored[3]) + "\t" +
                        chrom + "\t" + f[3] + "\t" + str(end) + "\n")
            o.write(sb)
            stored = []
    
    # Last line, if unpaired
    if stored:
        sb = (stored[0] + "\t.\t-1\t-1\t-1\t-1\t" + stored[1] + "\t" +
                str(stored[2]) + "\t" + str(stored[3]) + "\t.\t-1\t-1\n")
        o.write(sb)
    
    # Finish
    o.close()
    f.Close()
    PRINT.printP(STR__pair_complete)
    
    # Reporting
    metrics = [count_reads, count_pairs, count_pairs_same, total_gaps,
            total_frag_Ns]
    Report_Metrics(metrics)
    
    # Wrap up
    return 0

def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The number of reads
                * The number of read pairs
                * The number of read pairs for which both reads align to the
                    same chromosome
                * The sum of all the gaps between read pairs
                * The sum of all the fragment sizes
    
    Report_Metrics([int, int, int, int, int]) -> None
    """
    # Unpacking
    reads, pairs, pairs_same, gaps, frag_Ns = summary_metrics
    # Calculations
    if reads < 2:
        pct_paired = 0.0
        pcr_same_r = 0.0
    else:
        pct_paired = (100.0 * pairs)/(reads/2)
        pct_same_r = (100.0 * pairs_same)/(reads/2)
    if pairs == 0:
        pct_same_p = 0.0
    else:
        pct_same_p = (100.0 * pairs_same)/pairs
    if pairs_same == 0:
        avg_gaps = 0.0
        avg_frag = 0.0
    else:
        avg_gaps = (float(gaps))/pairs_same
        avg_frag = (float(frag_Ns))/pairs_same
    # Strings
    reads = str(reads) + "   "
    pairs = str(pairs) + "   "
    pairs_same = str(pairs_same) + "   "
    placeholder = "   "
    avg_gaps = str(avg_gaps)
    avg_frag = str(avg_frag)
    pct_paired = str(pct_paired)
    pct_same_r = str(pct_same_r)
    pct_same_p = str(pct_same_p)
    # Trim
    avg_gaps = Trim_Percentage_Str(avg_gaps, 2)
    avg_frag = Trim_Percentage_Str(avg_frag, 2)
    pct_paired = Trim_Percentage_Str(pct_paired, 2)
    pct_same_r = Trim_Percentage_Str(pct_same_r, 2)
    pct_same_p = Trim_Percentage_Str(pct_same_p, 2)
    # Pad column (1)
    col_1 = [reads, pairs, pairs_same, placeholder, avg_gaps, avg_frag]
    col_1 = Pad_Column(col_1, 0, 0, " ", 0)
    reads, pairs, pairs_same, placeholder, avg_gaps, avg_frag = col_1
    # Pad column (2)
    col_2 = [pct_paired, pct_same_r, pct_same_p]
    col_2 = Pad_Column(col_2, 0, 0, " ", 0)
    pct_paired, pct_same_r, pct_same_p = col_2
    # Print
    PRINT.printM(STR__metrics.format(A = reads, B = pairs, C = pct_paired,
            D = pairs_same, E = pct_same_r, F = placeholder, G = pct_same_p,
            H = avg_gaps, I = avg_frag))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Pair_SAM_Reads(raw_command_line_input):
    """
    Parse the command line input and call the Pair_SAM_Reads function with
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
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o"]:
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
        path_out = arg2
    
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
    exit_state = Pair_SAM_Reads(path_SAM, path_out)
    
    # Exit
    if exit_state == 0: return 0
    else: return exit_state



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
    exit_code = Parse_Command_Line_Input__Pair_SAM_Reads(sys.argv)


