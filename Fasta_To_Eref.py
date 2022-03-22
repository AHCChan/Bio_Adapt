HELP_DOC = """
FASTA TO EREF
(version 1.0)
by Angelo Chan

This is a program for tranforming a FASTA file containing one or more reads,
into an EREF file and its corresponding folder.

The EREF file is a 2-column TSV containing the names of reads and their
corresponding file paths. The corresponding folder is a foldering containing
multiple FASTA files, each of which contains a single sequence.



USAGE:
    
    python27 Fasta_To_Eref.py <input_FASTA> [-o <output_EREF> <output_folder>]
            [-p <path_prefix>]



MANDATORY:
    
    input_FASTA
        
        The filepath of the input file.

OPTIONAL:
    
    output_EREF
        
        (DEFAULT path generation available)
        
        The name of the EREF file generated. Filepaths are accepted.
        
        An EREF file is a TSV with two columns. Each row contains the data for
        one sequence. The first column contains the name of the sequence and the
        second column contains the filepath of the FASTA file containing that
        sequence.
    
    output_folder
        
        (DEFAULT path generation available)
        
        The name of the output folder generated. Directory paths are accepted.
        By default, this will be the prefix of the filepaths in the EREF file.
    
    path_prefix
        
        (DEFAULT option available)
        
        If supplied, this is will be the prefix of the filepaths in the EREF
        file.



EXAMPLES:
    
    python27 Fasta_To_Eref.py RepbaseConsensuses.fa -o RetroSeq_Eref.tsv
            RetroSeq_Eref_Folder

USAGE:
    
    python27 Fasta_To_Eref.py <input_FASTA> [-o <output_EREF> <output_folder>]
            [-p <path_prefix>]
"""

NAME = "Fasta_To_Eref.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True

FILEMOD__EREF = "__EREF.tsv"
FILEMOD__FOLD = "__CON_SEQS"



# Minor Configurations #########################################################

DEFAULT__width = 80 # Width of the output FASTA files

FILEMOD__FASTA_EXTENSION = ".fa"



# Imported Modules #############################################################

import sys
import os

import random as Random



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from FASTA_File_Reader import *

from Width_File_Writer import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python Fasta_To_Eref.py -h"

STR__metrics = """
            Total reads: {A}
            Total bases: {B}
    Average Read Length: {C}"""

STR__parsing_args = "\nParsing arguments..."

STR__f2e_begin = "\nRunning Fasta_To_Eref..."

STR__f2e_complete = "\nFasta_To_Eref successfully finished."



STR__unexpected_failure = "\nProgram exited with an unexpected error."

STR__overwrite_confirm_2 = """
Files may exist in destination folder:
    {f}
Do you wish to overwrite them in the event of a naming clash? (y/n): """



# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Fasta_To_Eref(path_in, path_eref, path_folder, path_mod):
    """
    Generate an EREF file and corresponding folder from a FASTA file.
    
    @path_in
            (str - filepath)
            The filepath of the input FASTA file.
    @path_eref
            (str - filepath)
            The filepath of the output EREF file.
            An EREF file is a TSV with two columns. Each row contains the data
            for one sequence. The first column contains the name of the sequence
            and the second column contains the filepath of the FASTA file
            containing that sequence.
    @path_folder
            (str - dirpath)
            The filepath of the output folder where the individual sequences
            will be output.
    @path_mod
            (str - dirpath)
            The path appended to the filename in the reference column of the
            EREF file.
            Assumed to end in either "/" or "\".
    
    Fasta_To_Eref(str, str, str, str) -> int
    """
    # Sanitize paths
    if path_folder[-1] not in ["/", "\\"]:
        path_folder += "/"
    if path_mod[-1] not in ["/", "\\"]:
        path_mod += "/"
    
    # Set up reporting
    reads = 0
    bases = 0
    
    # Setup the I/O
    try:
        f = FASTA_Reader()
        f.Open(path_in)
    except:
        return 1
    try:
        e = open(path_eref, "w")
    except:
        return 1
    w = Width_File_Writer()
    w.Overwrite_Allow()
    w.Set_Width(DEFAULT__width)
    w.Set_Newline("\n")
    w.Toggle_Printing_M(True)
    
    # Main loop
    PRINT.printP(STR__f2e_begin)
    while not f.End():
        f.Read()
        name = f.Get_Name()
        name = name.split("\t")
        name = name[0]
        name_ = name + FILEMOD__FASTA_EXTENSION
        seq = f.Get_Seq()
        # Metrics
        reads += 1
        bases += len(seq)
        # Paths
        path1 = path_mod + name_ # String which is written to file
        path2 = path_folder + name_ # Actual path of the FASTA output file
        # Eref
        e.write(name + "\t" + path1 + "\n")
        # File
        w.Open(path2)
        w.Write_F(">" + name)
        w.Newline()
        w.Write(seq)
        w.Close()
    PRINT.printP(STR__f2e_complete)
    
    # Close up
    e.close()
    f.Close()
    
    # Reporting
    Report_Metrics(reads, bases)

    # Wrap up
    return 0

def Report_Metrics(reads, bases):
    """
    Print a report into the command line interface of the results of running
    this program.
    
    @reads
            (int)
            The total number of reads processed.
    @bases)
            (int)
            The total number of bases in the reads processed.
    
    Report_Metrics(int,int) -> None
    """
    # Calculations
    avg = float(bases)/reads
    # Strings
    str_reads = str(reads) + "   "
    str_bases = str(bases) + "   "
    str_avg = str(avg) + "0"
    # Padding and formatting
    str_avg = Trim_Percentage_Str(str_avg, 2)
    max_size = max([len(str_reads), len(str_bases), len(str_avg)])
    str_reads = Pad_Str(str_reads, max_size, " ", 0)
    str_bases = Pad_Str(str_bases, max_size, " ", 0)
    str_avg = Pad_Str(str_avg, max_size, " ", 0) 
    # Print
    PRINT.printM(STR__metrics.format(A = str_reads, B = str_bases, C = str_avg))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__F2E(raw_command_line_input):
    """
    Parse the command line input and call the Fasta_To_Eref function with
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
    
    # Initial validation (Redundant in current version)
    if len(inputs) < 1:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory input
    path_in = inputs.pop(0)
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_eref = ""
    path_folder = ""
    path_mod = ""
    
    # Parse arguments
    while inputs:
        arg = inputs.pop(0)
        
        # Confirm valid flag
        if arg in ["-p"]: # Second argument
            try:
                arg2 = inputs.pop(0)
            except:
                PRINT.printE(STR__insufficient_inputs)
                PRINT.printE(STR__use_help)
                return 1
        elif arg in ["-o"]: # Second and third arguments
            try:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
            except:
                PRINT.printE(STR__insufficient_inputs)
                PRINT.printE(STR__use_help)
                return 1
        else: # Invalid
            arg = Strip_X(arg)
            PRINT.printE(STR__invalid_argument.format(s = arg))
            PRINT.printE(STR__use_help)
            return 1
        
        # Individual validation
        if arg == "-o": # Output files - Actual validation done later
            path_eref = arg2
            path_folder = arg3
        elif arg == "-p":
            path_mod = arg2
    
    # Defaults
    if not path_eref:
        path_eref = Generate_Default_Output_File_Path_From_File(path_in,
                FILEMOD__EREF, False)
        path_folder = Generate_Default_Output_File_Path_From_File(path_in,
                FILEMOD__FOLD , False)
    if not path_mod:
        path_mod = path_folder
    
    # Validate output paths
    valid_eref = Validate_Write_Path(path_eref)
    if valid_eref == 1: PRINT.printM(STR__overwrite_accept)
    elif valid_eref == 2: return 0
    elif valid_eref == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    elif valid_eref == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    valid_folder = Validate_Write_Path__FOLDER(path_folder)
    if valid_folder == 0: pass
    elif valid_folder == 1:
        PRINT.printM(STR__overwrite_accept)
    else:
        if valid_folder == 2: PRINT.printE(STR__IO_error_write_folder_cannot)
        if valid_folder == 3: PRINT.printE(STR__overwrite_decline)
        if valid_folder == 4: PRINT.printE(STR__IO_error_write_folder_forbid)
        if valid_folder == 5:
            PRINT.printE(STR__IO_error_write_folder_nonexistent)
        if valid_folder == 6: PRINT.printE(STR__IO_error_write_unexpected)
        return 1
    
    # Run program
    Fasta_To_Eref(path_in, path_eref, path_folder, path_mod)
    
    # Safe exit
    return 0

def Validate_Write_Path(filepath):
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
    exit_code = Parse_Command_Line_Input__F2E(sys.argv)
