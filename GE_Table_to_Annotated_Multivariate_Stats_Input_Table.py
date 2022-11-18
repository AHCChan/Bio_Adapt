HELP_DOC = """
GENE EXPRESSION TABLE TO ANNOTATED MULTIVARIATE STATISTICS INPUT TABLE
(version 1.0)
by Angelo Chan
This is a program for taking a folder of Gene Expression table files, combining
them, transposing it (swapping the X and Y axes), optionally adding annotations
from an annotation file.

It is assumed that the Gene Expression table files have the samples across the
X axis, and the genes across the Y axis. That is to say, the first row is a
header row containing the names of the samples, and each row afterwards contains
the values for one gene.

The annotation file is assumed for each row to containing the annotation data for
one sample. It will work with or without a header row.




USAGE:
    
    python27 GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py
            <input_folder> <file_identifier> <input_format> <shortlist_genes>
            [-a <annotation_file>] [-o <output_file>] [-h <header>]
            [-b <annotation_only>]



MANDATORY:
    
    input_folder
        
        The directory path of the folder containing the Gene Expression files.
    
    file_identifier
        
        The substring (i.e., the text) found in the file names and/or file
        extensions of the Gene Expression files.

        For example, if the folder contains files downloaded directly from
        Wormbase, without changing the file names, all files will begin with
        "WBPaper" and have a ".csv" file extension. Either of those could be
        used as the file identifier.
    
    input_format
        
        The file format of the input file. Acceptable options are:
            tsv - Tab-separated values
            csv - Comma-separated values
            ssv - Space-separated values

    shortlist_genes
        
        The file path of the file containing a list of all the genes which are
        universally found in all the files in the input folder, and which are
        desirable in the output file.
        
        The file should contain one gene per line, and nothing else.
        
        WB_GE_Genes_Overlap_Report.py can be used to obtain such a file for any
        given folder.



OPTIONAL:
    
    annotation_file
        
        The filepath of the annotation file containing the sample annotations.
        This annotation file needs to be a TSV, with the sample IDs in the
        first column. The file can have any number of columns, as well as a
        header row, although a header row is not required.
        
        If an annotation file contains all the desired data, but is in the wrong
        table format, contains undesired columns, or has the columns in the
        wrong order, Table_to_Table.py can be used to adjust it.
        (https://github.com/AHCChan/Table_To_Table)
    
    output_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file.

    header
        
        (DEFAULT: Y)
        
        Whether or not to include headers in the output file.

    annotated_only
        
        (DEFAULT: N)
        
        Whether or not to only write samples into the output file which have
        annotations.



EXAMPLE:
    
    python27 GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py Genes.txt
            Folder/Of/WB/Files WBPaper tsv -a Sample_Annotations.tsv -o
            Consolidated_WB_GEs.tsv -h Y -b Y
    
    python27 GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py Genes.txt
            Folder/Of/GE/Files .tsv tsv -o Simple_Consolidated_GEs.tsv -h N



USAGE:
    
    python27 GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py
            <input_folder> <file_identifier> <input_format> <shortlist_genes>
            [-a <annotation_file>] [-o <output_file>] [-h <header>]
            [-b <annotation_only>]
"""

NAME = "GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD_1 = "__Consolidated_GEs.tsv"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__Header = True
DEFAULT__Annotated_Only = False



# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *




# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\tpython "\
"GE_Table_to_Annotated_Multivariate_Stats_Input_Table.py -h"



STR__error_no_GE = "ERROR: No files containing the file idenfier string: \n\t "\
"{i}\n... were found in the folder: \n\t{f}"

STR__error_annotate_file_inconsistent = "ERROR: Annotation file does not "\
"contain a consistent number of columns."



STR__metrics = """
         Total Datasets: {A}
          Total Samples: {B}
            
      Shortlisted Genes: {C}

      Annotated Samples: {D}
    Unannotated Samples: {E}

        Written to File: {F}
"""

STR__report_begin = "\nRunning Consolidate_GE_Files..."

STR__report_complete = "\nConsolidate_GE_Files successfully finished."



STR__unexpected_failure = "\nProgram exited with an unexpected error."


# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Consolidate_GE_Files(input_folder, file_identifier, file_delimiter,
            shortlist_genes, annotation_file, output_file, header,
            annotated_only):
    """
    Consolidate the gene expression data in [input_folder] into a single file.
    
    Returns exit codes:
        0: No errors detected
        1: Annotation file has an inconsistent number of columns
    
    @input_folder
            (str - dirpath)
            The directory path of the folder containing the Gene Expression
            files.
    @file_identifier
            (str)
            The substring (i.e., the text) found in the file names and/or file
            extensions of the Gene Expression files.
            For example, if the folder contains files downloaded directly from
            Wormbase, without changing the file names, all files will begin with
            "WBPaper" and have a ".csv" file extension. Either of those could be
            used as the file identifier.
    @file_delimiter
            (str)
            The table delimiter. (Tabs for TSVs, Commas for CSVs, etc)
    @shortlist_genes
            (str - filepath)
            The file path of the file containing a list of all the genes which
            are universally found in all the files in the input folder, and
            which are desirable in the output file.
            The file should contain one gene per line, and nothing else.
    @annotation_file
            (str - filepath)
            The filepath of the annotation file containing the sample
            annotations. This annotation file needs to be a TSV, with the sample
            IDs in the first column. The file can have any number of columns, as
            well as a header row, although a header row is not required.
            An empty string can be supplied if there are no annotations.
    @output_file
            (str - filepath)
            The filepath of the output file.
    @header
            (bool)
            Whether or not to include headers in the output file.
    @annotated_only
            (bool)
            Whether or not to only write samples into the output file which have
            annotations.
            
    
    Consolidate_GE_Files(str, str, str, str, str, bool, bool) -> int
    """
    PRINT.printP(STR__report_begin)
    
    # Setup - Paths and Lists
    file_list = Get_Files_W_Substring(input_folder, file_identifier)
    
    gene_list = Col_No_To_Set([shortlist_genes], "TSV", 0)
    gene_list = sorted(gene_list)
    
    # Get data
    annotations = Get_Annotations(annotation_file)
    if annotations == 1: return 1 # Rows are if unequal length
    metrics = Write_GE_Data(file_list, file_delimiter, gene_list, annotations,
            output_file, header, annotated_only)
        
    PRINT.printP(STR__report_complete)
    
    # Reporting
    Report_Metrics(metrics)
    
    # Wrap up
    return 0



def Get_Samples(file_list):
    """
    Return a list of all the sample IDs contained in the Wormbase Gene
    Expression files in [file_list].
    
    @file_list
            (list<str - filepath>)
            A list containing the filepaths of all the files to be scanned for
            sample IDs.
    
    Get_Samples(list<str>) -> list<str>
    """
    results = []
    #
    f = Table_Reader()
    f.Set_Delimiter("\t")
    #
    for file_ in file_list:
        f.Set_New_Path(file_)
        f.Open()
        f.Read()
        headers = f.Get()
        samples = headers[1:]
        for sample in samples:
            results.append(sample)
        f.Close()
    #
    results = list(results)
    return results

def Get_Annotations(annotation_file):
    """
    Return a dictionary containing the annotation string for the samples named
    in the annotation file. The dictionary will also contain a value string
    corresponding to the key:int(0) which can be used for unannotated samples,
    while the column headers in the first line of the annotation file will also
    correspond to the key:int(1).
    
    @annotation_file
            (str - filepath)
            A TSV file containing the annotations.
    
    Get_Annotations(str) -> dict
    """
    if not annotation_file: return {0: "", 1: ""}
    # Setup
    results = {}
    f = open(annotation_file, "U")
    count = 0
    # First row
    line = f.readline()
    line = line.strip("\n")
    if not line: return {0: "", 1: ""}
    values = line.split("\t", 1)
    if len(values) != 2: # Single value
        count = 0
        results[values[0]] = ""
    else: # Contains annotations
        key, annotations = values
        count = annotations.count("\t")
        results[key] = annotations
        results[0] = line
        results[1] = (count+1)*("\t")
    # Rest of the file
    while line:
        values = line.split("\t", 1)
        if len(values) != 2: # Single value
            temp_count = 0
            results[values[0]] = ""
        else: # Contains annotations
            key, annotations = values
            temp_count = annotations.count("\t")
            results[key] = annotations
        if temp_count != count: return 1
        # Read next
        line = f.readline()
        line = line.strip("\n")
    # Finish up and return
    f.close()
    return results



def Write_GE_Data(file_list, file_delimiter, gene_list, annotations,
            output_file, header, annotated_only):
    """
    Write the gene expression data in [file_list] and write it to file.
    
    Return a list of integers which represent the following metrics:
        * The total number of datasets
        * The total number of samples
        * The total number of genes in the shortlist file
        * The number of samples with annotations
        * The number of samples without annotations
        * The number of samples written to file
    
    @file_list
            (list<str - filepath>)
            A list containing the filepaths of all the files to be scanned for
            sample IDs.
    @file_delimiter
            (str)
            The table delimiter. (Tabs for TSVs, Commas for CSVs, etc)
    @gene_list
            (list<str>)
            A list of the genes to be written to file.
    @annotations
            (dict<str:str>)
            A dictionary containing the annotations for the samples.
            The key:int(0) corresponds to the header row values.
            The key:int(1) corresponds to the placeholder string to be used if
            a particular sample cannot be found in this dictionary.
    @output_file
            (str - filepath)
            The filepath of the output file.
    @header
            (bool)
            Whether or not to include headers in the output file.
    @annotated_only
            (bool)
            Whether or not to only write samples into the output file which have
            annotations.
    
    Write_GE_Data(list<str>, list<str>, dict, str, bool, bool) -> list<int>
    """
    # Metrics
    n_datasets = 0
    n_samples = 0
    n_genes = len(gene_list)
    n_annotated = 0
    n_unannotated = 0
    n_written = 0    
    # Setup
    data = {}
    o = open(output_file, "w")
    f = Table_Reader()
    f.Set_Delimiter(file_delimiter)
    gene_set = set(gene_list)
    # Write - header
    if header:
        o.write(annotations[0])
        sb = ""
        for gene in gene_list:
            sb += "\t" + gene
        o.write(sb)
        o.write("\n")
    # Go through all files
    for file_ in file_list:
        # Open
        f.Set_New_Path(file_)
        f.Open()
        # Header
        f.Read()
        headers = f.Get()
        headers = headers[1:]
        # Setup
        length = len(headers)
        range_ = range(length)
        mapper = {}
        for i in range_:
            data[i] = {}
            mapper[i] = headers[i]
        # Populate with data
        while not f.EOF:
            f.Read()
            columns = f.Get()
            gene = columns[0]
            values = columns[1:]
            if gene in gene_set:
                for i in range_:
                    data[i][gene] = values[i]
        # Close
        f.Close()
        # Write - body
        for i in range_:
            flag = True
            sample = mapper[i]
            if sample in annotations:
                n_annotated += 1
                annotation = "\t" + annotations[sample]
            else:
                n_unannotated += 1
                annotation = annotations[1]
                if annotated_only: flag = False
            if flag:
                n_written += 1
                o.write(sample)
                o.write(annotation)
                sb = ""
                for gene in gene_list:
                    sb += "\t" + data[i][gene]
                o.write(sb)
                o.write("\n")
        # Other metrics
        n_datasets += 1
        n_samples += length
    # Return metrics
    return [n_datasets, n_samples, n_genes, n_annotated, n_unannotated,
            n_written]



def Report_Metrics(metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The total number of datasets
                * The total number of samples
                * The total number of genes in the shortlist file
                * The number of samples with annotations
                * The number of samples without annotations
                * The number of samples written to file
    
    Report_Metrics(list<int>(6)) -> None
    """
    # Unpacking
    datasets = metrics[0]
    samples = metrics[1]
    genes = metrics[2]
    annotated = metrics[3]
    unannotated = metrics[4]
    written = metrics[5]
    # Print
    PRINT.printM(STR__metrics.format(A = datasets, B = samples, C = genes,
            D = annotated, E = unannotated, F = written))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Consolidate_GE_Files(raw_command_line_input):
    """
    Parse the command line input and call the Consolidate_GE_Files function with
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
    if len(inputs) < 4:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    
    # Set up rest of the parsing
    annotation_file = ""
    path_out = ""
    header = DEFAULT__Header
    annotation_only = DEFAULT__Annotated_Only
    
    # Initial validation
    
    # Automated output path generation
    
    # Validate output path
    
    # Run program
    exit_state = Consolidate_GE_Files(file_list, file_delimiter, gene_list,
            shortlist_file, annotation_file, path_out, header, annotation_only)
    
    # Safe exit
    if exit_state == 0: return 0
    else:
        if exit_state == 1:
            PRINT.printE(STR__error_annotate_file_inconsistent)
        PRINT.printE(STR__use_help)
        return 1



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Consolidate_GE_Files(sys.argv)
