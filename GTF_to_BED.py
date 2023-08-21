HELP_DOC = """
GTF TO BED CONVERTER
(version 1.0)
by Angelo Chan

This is a program for taking a particular kind of GTF file and converting it
into BED file. The input GTF file this program was designed for is the kind
which contains start codons, stop codons, and other elements, while the output
BED file will be a 5 column TSV containing the genomic coordinates (first three
columns), the name of the gene (4th column), and the directionality of the
sequence. (5th column)

The start and stop codon coordinates will be used to determine the start and end
of the coding sequence of the gene.

* THIS IS NOT THE SAME AS THE START/STOP OF THE MRNA SEQUENCE. *



USAGE:
    
    python27 GTF_to_BED.py <input_file> <field> [-o <output_file>] [-h
            <header_in> <header_out>] [-n <no_ID_file>]



MANDATORY:
    
    input_file
        
        The GTF file were each line is an entry for either an exon, a CDS, a
        start codon, or a stop codon. The GTF is assumed to be a TSV (tab
        separated values) where each row records the details of a segment of DNA
        with the following columns:
            
            1)  Chromosome name (genomic coordinates)
            2)  Biological function of the overall RNA sequence
                (protein coding, tRNA, ncRNA, etc)
                OR
                The name of the program which identified this element.
            3)  Biomolecular function of the specific segment of DNA
                (start codon, exon, CDS, stop codon, etc)
            4)  Upstream nucleotide number (genomic coordinates)
                  - Sometimes mistakenly referred to as the "start", or starting
                    nucleotide number. However, this is only true for with a "+"
                    directionality.
            5)  Downstream nucleotide number (genomic coordinates)
                  - Sometimes mistakenly referred to as the "end", or ending
                    nucleotide number. However, this is only true for with a "+"
                    directionality.
            6)  Placeholder column containing a full stop.
                OR
                A decimal number denoting the score which was calculated by the
                program which identified this element.
            7)  The directionality of the DNA element. "+" for an element on the
                forward strand of DNA (starts at a lower-numbered genomic
                coordinate and ends at a higher-numbered one), and "-" for an
                element on the backward strand of DNA (starts at a
                higher-numbered genomic coordinate and ends at a lower-numbered
                one).
                Also known as the strand.
            8)  The reading frame, relative to the first start codon:
                    0 -  In the same frame
                    1 -  Shifted forward one basepair
                    2 -  Shifted forward two basepairs
            9)  A series of tags, or semicolon-separated field-value entries.
                This program requires that all entries contain the "gene_id" tag
                and tolerates other tags being present.
                Other tags may include "transcript_if", "protein_id", etc.
    
    field
        
        The data field used for the name/ID. Examples include:
            gene_name
            gene_id
            transcript_id



OPTIONAL:
    
    output_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file.
    
    header_in
        
        (DEFAULT: N)
        
        Whether or not the GFT file contains headers.
    
    header_out
        
        (DEFAULT: N)
        
        Whether or not the output file should contain headers.
    
    no_ID_file
        
        (DEFAULT: None)
        
        The filepath of the output file for the entries which don't have a valid
        ID.
        If no filepath is specified, this data will be ignored.



EXAMPLE:
    
    python27 GTF_to_BED.py path/organism.gtf
    
    python27 GTF_to_BED.py path/organism.gtf -o path/coding_regions.bed -h N Y



USAGE:
    
    python27 GTF_to_BED.py <input_file> <field> [-o <output_file>] [-h
            <header_in> <header_out>] [-n <no_ID_file>]
"""

NAME = "GTF_to_BED.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__Coding_Regions.tsv"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__headers_in = False
DEFAULT__headers_out = False




# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from GTF_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"GTF_to_BED.py -h"



STR__column_header = "GENE_ID\tTOTAL_DATASETS\tTOTAL_SAMPLES"



STR__metrics = """
                        Total Elements: {A}
                                Coding: {B} ({C}%)
                            Non-Coding: {D} ({E}%)
                            No Name/ID: {F} ({G}%)

    Average Length of Coding Sequences: {H}
"""



STR__convert_begin = "\nRunning GTF_to_BED..."

STR__convert_complete = "\nWB_GE_Genes_Overlap_Report successfully finished."



# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Convert_GTF_to_BED(path_in, field, path_out, header_in, header_out,
            path_no_ID):
    """
    Convert a GTF file into a BED file, for GTF files which contain the
    coordinates of genetic elements like start codons, stop codons, and exons.
    
    Output BED file is 5-column TSV with the genomic coordinates in the first
    three columns, the gene name in the fourth, and the strand/directionality in
    the fifth.
    
    @path_in
            (str - filepath)
            The filepath of the input GTF file.
    @field
            (str)
            The name of the data field used as the name/ID of the resulting
            sequences.
    @path_out
            (str - filepath)
            The filepath of the output BED file.
    @header_in
            (bool)
            Whether or not there are headers in the input GTF file.
    @header_out
            (bool)
            Whether or not there should be headers in output BED file.
    @path_out
            (str - filepath)
            The filepath of the output BED file.
    
    Convert_GTF_to_BED(str, str, bool, bool) -> int
    """
    # Setup reporting
    total_entries = 0
    total_coding = 0
    total_span = 0
    total_NC = 0
    total_no_ID = 0
    
    # Setup the I/O
    f = GTF_Reader()
    f.Set_Grouping_Method("TAG")
    f.Set_Tag(field)
    f.Open(path_in)
    
    o = open(path_out, "w")

    n = None
    if path_no_ID: n = open(path_no_ID, "w")
    
    # Main loop
    PRINT.printP(STR__convert_begin)
    while not f.End():
        total_entries += 1
        # Read
        f.Read()
        gene = f.Get_Current_ID()
        coords = f.Get_Coords()
        # Non-coding
        if not coords: total_NC += 1
        # No ID
        elif not gene:
            total_no_ID += 1
            if n:
                raw = f.Get()
                for values in raw:
                    sb = "\t".join(values[:9]) + "\n"
                    n.write(sb)
        # Coding
        else:
            # Unpack
            chr_, start, end, strand, length = coords
            start = str(start)
            end = str(end)
            # Metrics
            total_coding += 1
            total_span += length
            # Write
            sb = (chr_ + "\t" + start + "\t" + end + "\t" + gene + "\t" +
                        strand + "\n")
            o.write(sb)
    PRINT.printP(STR__convert_complete)
    
    # Close up
    if n: n.close()
    o.close()
    f.Close()
    
    # Reporting
    Report_Metrics([total_entries, total_coding, total_span, total_NC,
            total_no_ID])

    # Wrap up
    return 0

def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The total number of DNA elements in the GTF file
                * The total number of coding sequences
                        (has at least one start and stop codon)
                * The total number of nucleotides spanned by the coding seqs
                * The total number of non-coding DNA elements
                * The total number of DNA elements without a valid ID
    
    Report_Metrics([int, int, int, int]) -> None
    """
    # Unpacking
    total, coding, span, NC, no_ID = summary_metrics
    # Calculations
    percentage_coding = (coding*100.0)/total
    percentage_NC = (NC*100.0)/total
    percentage_no_ID = (no_ID*100.0)/total
    avg_span = (float(span))/coding
    # Strings
    total = str(total)
    coding = str(coding)
    NC = str(NC)
    no_ID = str(no_ID)
    avg_span = str(avg_span)
    percentage_coding = str(percentage_coding)
    percentage_NC = str(percentage_NC)
    percentage_no_ID = str(percentage_no_ID)
    # Pad Column 1
    total = total + "   "
    coding = coding + "   "
    NC = NC + "   "
    no_ID = no_ID + "   "
    avg_span = Trim_Percentage_Str(avg_span, 2)
    # 
    col_1 = [total, coding, NC, no_ID, avg_span]
    col_1 = Pad_Column(col_1, 0, 0, " ", 0)
    total, coding, NC, no_ID, avg_span = col_1
    # Pad percentages
    percentage_coding = Trim_Percentage_Str(percentage_coding, 2)
    percentage_NC = Trim_Percentage_Str(percentage_NC, 2)
    percentage_no_ID = Trim_Percentage_Str(percentage_no_ID, 2)
    #
    percentages = [percentage_coding, percentage_NC, percentage_no_ID]
    percentages = Pad_Column(percentages, 0, 0, " ", 0)
    percentage_coding, percentage_NC, percentage_no_ID = percentages
    # Print
    PRINT.printM(STR__metrics.format(A = total, B = coding,
            C = percentage_coding, D = NC, E = percentage_NC, F = no_ID,
            G = percentage_no_ID, H = avg_span))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__GTF_to_BED(raw_command_line_input):
    """
    Parse the command line input and call the Convert_GTF_to_BED function
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
    field = inputs.pop(0)
    
    # Set up rest of the parsing
    path_out = ""
    header_in = DEFAULT__headers_in
    header_out = DEFAULT__headers_out
    path_no_ID = ""
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        try: # Following arguments
            if arg in ["-o", "-n"]:
                arg2 = inputs.pop(0)
            elif arg in ["-h"]:
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
        # Flag-dependent response
        if arg == "-o":
            path_out = arg2
        elif arg == "-n":
            path_no_ID = arg2
        else: # arg == "-h"
            header_in = Validate_Bool(arg2)
            if header_in == None:
                printE(STR__invalid_bool.format(s = arg2))
                return 1
            header_out = Validate_Bool(arg3)
            if header_in == None:
                printE(STR__invalid_bool.format(s = arg3))
                return 1
    
    # Automated output path generation
    if not path_out: Generate_Default_Output_File_Path_From_File(path_in,
            FILEMOD, False)
    
    # Validate output paths
    valid_out = Validate_Write_Path__FILE(path_out)
    if valid_out == 2: return 0
    if valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    if path_no_ID:
        valid_out = Validate_Write_Path__FILE(path_no_ID)
        if valid_out == 2: return 0
        if valid_out == 3:
            PRINT.printE(STR__IO_error_write_forbid)
            return 1
        if valid_out == 4:
            PRINT.printE(STR__IO_error_write_unable)
            return 1
    
    # Run program
    exit_state = Convert_GTF_to_BED(path_in, field, path_out, header_in,
            header_out, path_no_ID)
    
    # Exit
    if exit_state == 0: return 0
    else:
        if exit_state == 1: PRINT.printE(STR__unexpected_failure)
        PRINT.printE(STR__use_help)
        return 1



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
    exit_code = Parse_Command_Line_Input__GTF_to_BED(sys.argv)


