HELP_DOC = """
EXTRACT FLANKING
(version 2.0)
by Angelo Chan

(Modified from Sequence_Extractor.py, v1.0)

This is a program for copying the flanking sequences of genetic elements.

The sequences will be excised according to a table of genomic coordinates and
output into an output folder. The details of the extracts will be output in a
separate table file.



USAGE:
    
    python27 Extract_Flanking.py <genome_folder> <target_coordinates_table>
            [-o <extracted_sequences_folder> <coordinates_table>]
            [-w <window_size>] [-p <padding_size> <padding_basepair>]



MANDATORY:
    
    genome_folder
        
        The filepath of the input folder containing the template FASTA file(s).
        Each FASTA file is assumed to only contain one DNA sequence per file.
    
    target_coordinates_table
        
        The filepath of the input file containing the coordinates for extracting
        the target sequences from. The file needs to be in a TSV format and
        the first four columns need to contain the genomic coordinates:
            
            1) Chromosome name
            2) Start
            3) End
            4) Directionality* (+/-)
            
            * If the value in the Directionality column is neither "+" nor "-",
              the program will treat it as "+" by default. In cases where no
              directionality value is available, filling this column with random
              or irrelevant values will allow all entries to be treated as
              being on the forward strand of the "genome".

OPTIONAL:
    
    extracted_sequences_folder
        
        (DEFAULT path generation available)
        
        The filepath of the output folder where resultant excised sequences will
        be outputted to.
    
    coordinates_table
        
        (DEFAULT path generation available)
        
        The filepath of the output Coordinates Table file. The file will contain
        the details of the sequences extracted.
    
    window_size
        
        (DEFAULT: 50)
        
        The number of genomic nucleotides on either side of the element to
        extract.
    
    padding_size
        
        (DEFAULT: 100)
        
        The number of nucleotides to "pad" either side of the genomic
        nucleotides with.
    
    padding_basepair
        
        (DEFAULT: A)
        
        The nucleotide to pad the genomic nucleotides with.



EXAMPLES SCENARIO EXPLANATION:
    
    1:
    A modified RMSK file is used to provide coordinates for the Transposons in a
    genome, which are cut out of it.
    
    2:
    A modified RMSK file is used to provide coordinates for the Transposons in a
    genome, which are cut out of it. The data has shorter reads, so the window
    is smaller, as is the padding.

EXAMPLES:
    
    python27 Extract_Flanking.py Path/GenomeFolder rmsk__MOD.tsv
    
    python27 Extract_Flanking.py Path/GenomeFolder rmsk__MOD.tsv -o
            Path/Flanking_Regions -w 25 -p 50 A

USAGE:
    
    python27 Extract_Flanking.py <genome_folder> <target_coordinates_table>
            [-o <extracted_sequences_folder> <coordinates_table>]
            [-w <window_size>] [-p <padding_size> <padding_basepair>]
"""

NAME = "Extract_Flanking.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

DIRMOD = "__FLANKING"
FILEMOD__DETAILS = "__FLANKING_DETAILS.tsv"
FILEMOD__FLANKING = "__FLANKING.fa"

# For name string
ID_BASE = "TE_"
ID_SIZE = 9



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__width = 80

DEFAULT__window = 50
DEFAULT__pad_size = 75
DEFAULT__pad_str = "A"



# Imported Modules #############################################################

import sys
import os

import random as Random



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Chr_FASTA_File_Reader import *
from Table_File_Reader import *
from Width_File_Writer import *

from NSeq_Match import * #2.2



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"Extract_Flanking.py -h"

STR__error_no_FASTA = """
ERROR: No FASTA files detected in:
    {f}"""

STR__error_no_chr = """
ERROR: Unable to open chromosome FASTA file:
    {c}"""

STR__invalid_padding_size = """
ERROR: Invalid padding size:
    {s}
Please specify a non-negative integer.
"""

STR__invalid_padding_char = """
ERROR: Invalid padding char:
    {s}
Please specify a single character.
"""

STR__invalid_window_size = """
ERROR: Invalid window size:
    {s}
Please specify a non-negative integer.
"""



STR__metrics = """
         Chromosomes: {A}
    Basepairs copied: {B}
            Overlaps: {C}
    Sequences copied: {D}
"""

STR__extract_begin = "\nRunning Extract_Flanking..."

STR__extract_complete = "\nExtract_Flanking successfully finished."



STR__unexpected_failure = "\nProgram exited with an unexpected error."

STR__overwrite_confirm_2 = """
Files may exist in destination folder:
    {f}
Do you wish to overwrite them in the event of a naming clash? (y/n): """



# OS Strings ###################################################################

if sys.platform[:3] == "win":
    directory_spacer = "\\"
else:
    directory_spacer = "/"


# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Extract_Flanking(input_genome, input_coordinates, output_sequences,
            output_coordinates, window_size, padding_size, padding_char):
    """
    Copy the flanking sequences of genetic elements.

    The sequences will be excised according to a table of genomic coordinates
    and output into an output folder. The details of the extracts will be output
    in a separate table file.
    
    @input_genome
            (str - dirpath)
            The filepath of the folder containing the FASTA file(s) from which
            the sequences are to be extracted. Each file should only contain a
            single DNA sequence.
            The folder will typically be a "genome", with each of the files
            within being a chromosome, but does not absolutely have to be a
            "genome".
    @input_coordinates
            (str - filepath)
            The file containing the genomic coordinates and auxiliary
            information of the DNA sequences to be excised.
    @output_sequences
            (str - dirpath)
            A folder containing the sequences which were extracted. Each
            sequence will be stored in a separate FASTA file. Each sequence will
            also be given a new unique ID, which will also be used as the file
            name, and the sequence name. This ID will also be referenced in the
            resulting [output_coordinates] file, and subsequently derived files.
    @output_coordinates
            (str - filepath)
            The file containg the coordinates and details of the extracted
            sequences. It can be used to reinsert the extracted sequences back
            into the template. Derivatives of this file can be made which list
            genetic changes to the sequences, including substitutions,
            insertions, deletions, inversions, duplications, and fragmentations.
            These changes are listed in a way which is fairly human-readable,
            and without the need to create multiple folders containing highly
            similar information.
    @window_size
            (int)
            The number of genomic nucleotides on either side of the element to
            extract.
    @padding_size
            (int)
            The number of nucleotides to "pad" either side of the genomic
            nucleotides with.
    @padding_char
            (str)
            The character to pad the genomic nucleotides with.
            
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem accessing the data or if there are
            no valid FASTA files in the input genome folder.
    Return a value of 2 if there is a problem with the output file.
    Return a value of 3 if there is a problem during the sequence extraction
            process.
    Return a value of 4 if there are no elements specified in the coordinates
            file.
    
    Extract_Sequences(str, str, str, str, int, int, str) -> int
    """
    # Pad
    pad_str = padding_size * padding_char
    
    # Setup reporting
    chromosomes = 1
    basepairs_original = 0
    basepairs_copied = 0
    overlaps = 0
    seqs_copied = 0
    
    # Setup the I/O
    current_chr_name = ""
    seqs_current = []
    seq_next = []
    #
    f = Chr_FASTA_Reader()
    old_end = -1
    current_index = -1
    non_direction_flag = False # For when an entry has no +/-
    sb = ""
    #
    t = Table_Reader(input_coordinates)
    t.Set_Delimiter("\t")
    #
    c = open(output_coordinates, "w")
    o = Width_File_Writer()
    o.Overwrite_Allow()
    o.Set_Width(DEFAULT__width)
    o.Set_Newline("\n")
    o.Toggle_Printing_M(False)
    
    # Start
    PRINT.printP(STR__extract_begin)
    t.Open()
    remaining_flag = True
    
    # First entry
    t.Read()
    elements = t.Get_Current()
    if not elements:
        c.close()
        t.Close()
        return 4
    chr_name = elements[0]
    start = int(elements[1]) - window_size
    start_ = start - 1
    end = int(elements[2]) + window_size
    if start < 0: start = 0
    if start_ < 0: start_ = 0
    direction = elements[3]
    extras = elements[4:]
    #
    current_chr_name = chr_name
    seqs_current = [[chr_name, start, start_, end, direction, extras, "",
            elements]]
    
    # Open chromosome file
    chr_file_path = Get_Chr_File_Path(input_genome, chr_name)
    f.Open(chr_file_path)
    if f.End():
        c.close()
        f.Close()
        t.Close()
        PRINT.printE(STR__error_no_chr.format(c=chr_file_path))
        return 1
    
    # Second entry
    if not t.EOF:
        t.Read()
        elements = t.Get_Current()
        chr_name = elements[0]
        start = int(elements[1]) - window_size
        start_ = start - 1
        end = int(elements[2]) + window_size
        if start < 0: start = 0
        if start_ < 0: start_ = 0
        direction = elements[3]
        extras = elements[4:]
        #
        seq_next = [chr_name, start, start_, end, direction, extras, "",
                elements]
        
    # Main loop
    read_flag = True # In a rush, so this wasn't done via a special Class
    while read_flag:
        # Refill
        continue_flag = True
        if not seqs_current:
            seqs_current = [seq_next]
            t.Read()
            elements = t.Get_Current()
            if elements:
                chr_name = elements[0]
                start = int(elements[1]) - window_size
                start_ = start - 1
                end = int(elements[2]) + window_size
                if start < 0: start = 0
                if start_ < 0: start_ = 0
                direction = elements[3]
                extras = elements[4:]
                #
                seq_next = [chr_name, start, start_, end, direction, extras, "",
                        elements]
                #
                if current_chr_name != chr_name:
                    current_chr_name = chr_name
                    f.Close()
                    f.Open
                    chr_file_path = Get_Chr_File_Path(input_genome, chr_name)
                    f.Open(chr_file_path)
                    if f.End():
                        c.close()
                        f.Close()
                        t.Close()
                        PRINT.printE(STR__error_no_chr.format(c=chr_file_path))
                        return 1
                    current_index = -1
                    chromosomes += 1
            else:
                read_flag = False
                continue_flag = False
        # Latest
        earliest = seqs_current[0][2]
        latest = seqs_current[0][3]
        # Add coordinates until out of range
        while read_flag and continue_flag:
            # Read
            if seq_next:
                if seq_next[0] == current_chr_name:
                    if seq_next[1] <= latest:
                        seqs_current.append(seq_next)
                        if not t.EOF:
                            t.Read()
                            elements = t.Get_Current()
                            chr_name = elements[0]
                            start = int(elements[1]) - window_size
                            start_ = start - 1
                            end = int(elements[2]) + window_size
                            if start < 0: start = 0
                            if start_ < 0: start_ = 0
                            direction = elements[3]
                            extras = elements[4:]
                            #
                            current_chr_name = chr_name
                            seq_next = [chr_name, start, start_, end, direction,
                                    extras, "", elements]
                            # Update
                            if start_ < earliest: start = earliest
                            if end > latest: latest = end
                        else:
                            read_flag = False
                            continue_flag = False
                    else:
                        continue_flag = False
                else:
                    continue_flag = False
            else:
                read_flag = False
                remaining_flag = False
        # Read until start
        while current_index < earliest:
            f.Read()
            current_index += 1
        # Read until end
        while current_index < latest:
            f.Read()
            current_index += 1
            nuc = f.Get()
            overlap_temp = -1
            for seq in seqs_current:
                if current_index >= seq[1] and current_index < seq[3]:
                    overlap_temp += 1
                    seq[6] += nuc
            if overlap_temp > -1:
                if overlap_temp > 0:
                    overlaps += overlap_temp
        # Process elements
        for seq in seqs_current:
            seqs_copied += 1
            ID = Generate_Seq_ID(seqs_copied)
            if seq[4] == "-":
                seq[6] = Get_Complement(seq[6], True)
            sb = "\t".join([ID, seq[0], str(seq[1]), str(seq[3]), seq[4]] +
                    seq[5]) + "\n"
            c.write(sb)
            #
            extracted_seq = (seq[6][:window_size] + seq[6][-window_size:])
            padded_seq = (pad_str + extracted_seq + pad_str)
            path = output_sequences + "/" + ID + FILEMOD__FLANKING
            elements = seq[7]
            o.Open(path)
            o.Write_F(">" + ID + "\t" + "\t".join(elements))
            o.Newline()
            o.Write(padded_seq)
            o.Close_Newline()
            #
            basepairs_copied += len(extracted_seq)
        seqs_current = []
    
    PRINT.printP(STR__extract_complete)
    
    # Close up
    c.close()
    f.Close()
    t.Close()
    
    # Reporting
    Report_Metrics([chromosomes, basepairs_copied, overlaps, seqs_copied])

    # Wrap up
    return 0

def Get_Chr_File_Path(genome_folder_path, chr_name):
    """
    Return the file path to the Chromosomal FASTA file with [chr_name] as its
    name from the directory [genome_folder_path].
    Return an empty string if no matching file name is found.
    """
    names = os.listdir(genome_folder_path)
    for name in names:
        first = name.split(".")[0]
        if first == chr_name:
            filepath = genome_folder_path + directory_spacer + name
            return filepath
    return ""
    
def Generate_Seq_ID(counter):
    """
    Generate a sequence ID for a DNA sequence based on how many sequences have
    already been generated.
    """
    sb = Pad_Str(str(counter), ID_SIZE, "0", 0)
    sb = ID_BASE + sb
    return sb

def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The number of chromosomes from which sequeces were extracted
                * The total number of unique basepairs copied
                * The total number of overlaps
                * The total number of sequences copied
    
    Report_Metrics([int, int, int, int]) -> None
    """
    # Unpacking
    chromosomes, basepairs, overlaps, sequences = summary_metrics
    # Strings
    chromosomes = str(chromosomes)
    basepairs = str(basepairs)
    overlaps = str(overlaps)
    sequences = str(sequences)
    # Pad Column 1
    col_1 = [chromosomes, basepairs, overlaps, sequences]
    col_1 = Pad_Column(col_1, 0, 0, " ", 0)
    chromosomes, basepairs, overlaps, sequences = col_1
    # Print
    PRINT.printM(STR__metrics.format(A = chromosomes, B = basepairs,
            C = overlaps, D = sequences))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Extract_Flanking(raw_command_line_input):
    """
    Parse the command line input and call the Extract_Sequences function
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
    path_in_folder = inputs.pop(0)
    valid = Validate_FASTA_Folder(path_in_folder)
    if valid == 1:
        PRINT.printE(STR__IO_error_read_folder)
        PRINT.printE(STR__use_help)
        return 1
    elif valid == 2:
        PRINT.printE(STR__error_no_FASTA.format(f = path_in_folder))
        PRINT.printE(STR__use_help)
        return 1
    path_in_file = inputs.pop(0)
    valid = Validate_Read_Path(path_in_file)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in_file))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out_seqs = Generate_Default_Output_File_Path_From_File(path_in_file,
            DIRMOD, False)
    path_out_coords = Generate_Default_Output_File_Path_From_File(path_in_file,
            FILEMOD__DETAILS, False)
    window_size = DEFAULT__window
    padding_size = DEFAULT__pad_size
    padding_char = DEFAULT__pad_str
        
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        try: # Following arguments
            if arg in ["-o", "-p"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
            elif arg in ["-w"]:
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
        # Flag-dependent response
        if arg == "-o":
            path_out_seqs = arg2
            path_out_coords = arg3
        elif arg == "-p":
            padding_size = Validate_Int_NonNeg(arg2)
            padding_char = arg3
            if padding_size == -1:
                PRINT.printE(STR__invalid_padding_size.format(s = arg2))
                return 1
            if len(padding_char) != 1:
                PRINT.printE(STR__invalid_padding_char.format(s = arg3))
                return 1
        else: # arg == "-w"Validate_Int_NonNeg(arg2)
            window_size = Validate_Int_NonNeg(arg2)
            if window_size == -1:
                PRINT.printE(STR__invalid_window_size.format(s = arg2))
                return 1

    # Validate output paths
    valid_out = Validate_Write_Path__FOLDER(path_out_seqs)
    if valid_out == 0: pass
    elif valid_out == 1: PRINT.printM(STR__overwrite_accept)
    else:
        if valid_out == 2:
            PRINT.printE(STR__IO_error_write_folder_cannot)
        elif valid_out == 3: PRINT.printE(STR__overwrite_decline)
        elif valid_out == 4:
            PRINT.printE(STR__IO_error_write_folder_forbid)
        elif valid_out == 5:
            PRINT.printE(STR__IO_error_write_folder_nonexistent)
        elif valid_out == 6: PRINT.printE(STR__read_file_invalid)
        elif valid_out == 7:
            PRINT.printE(STR__IO_error_write_unexpected)
        return 1
    valid_out = Validate_Write_Path__FILE(path_out_coords)
    if valid_out == 2: return 0
    elif valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    elif valid_out == 4:
        PRINT.printE(STR__In_error_write_unable)
        return 1
    
    # Run program
    exit_state = Extract_Flanking(path_in_folder, path_in_file, path_out_seqs,
            path_out_coords, window_size, padding_size, padding_char)
    
    # Exit
    if exit_state == 0: return 0
    else:
        if exit_state == 1: PRINT.printE(STR__unexpected_failure)
        PRINT.printE(STR__use_help)
        return 1
    
    

def Validate_FASTA_Folder(dirpath):
    """
    Validates the dirpath of the input file as containing FASTA files.
    Return 0 if the dirpath is valid and contains at least 1 FASTA file.
    Return 1 if the dirpath is invalid.
    Return 2 if the dirpath is valid but contains no FASTA files.
    
    Validate_Read_Path(str) -> int
    """
    try:
        files = Get_Files_W_Extensions(dirpath, LIST__FASTA)
        if len(files) > 0: return 0
        return 2
    except:
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
    random_path = folder_path + directory_spacer + random_name
    while os.path.exists(random_path):
        random_name = str(Random.random())
        random_path = folder_path + directory_spacer + random_name
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
    exit_code = Parse_Command_Line_Input__Extract_Flanking(sys.argv)


