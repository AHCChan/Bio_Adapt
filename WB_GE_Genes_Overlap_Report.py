HELP_DOC = """
WORMBASE GENE EXPRESSION GENES OVERLAP REPORT
(version 1.0)
by Angelo Chan

This is a program for taking a folder of Wormbase Gene Expression files and
generating a report on the genes listed within those files.

The report will state the number of datasets, the number of samples, the number
of genes found in all samples, some other information, and for each gene, state
the number of datasets and samples which contain said gene.



USAGE:
    
    python27 WB_GE_Genes_Overlap_Report.py <input_folder> [-o <report_file>
            <universal_file>] [-a <abridge>]



MANDATORY:
    
    input_folder
        
        The directory path of the folder containing the raw Wormbase Gene
        Expression file downloads. It is assumed that all valid data files
        within said folder are prefixed by "WBPaper" and have a ".csv" file
        extension, but are formatted as TSVs.
        (Consistent with the download results from October, 2022)
        
        For the most reliable metrics, it should be ensured that no unique
        biological sample occurs in more than one dataset.



OPTIONAL:
    
    report_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file where the overall report will be output
        to. The report file will contain a summary of the results at the top
        of the file, followed by a table showing which genes are present in
        which datasets.
    
    shortlist_file
        
        (DEFAULT path generation available)
        
        The filepath of the output file where a shortlist of genes will be
        output to. These are the genes present in all datasets.
    
    abridge
        
        (DEFAULT: N)
        
        Whether or not to abridge the results to, for each gene, only how many
        datasets and samples a gene occurs in.



EXAMPLE:
    
    python27 WB_GE_Genes_Overlap_Report.py Path/WB_GE_Results -o
            Genes_Summary.tsv



USAGE:
    
    python27 WB_GE_Genes_Overlap_Report.py <input_folder> [-o <output_file>]
"""

NAME = "WB_GE_Genes_Overlap_Report.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD_1 = "__Overlap_Summary.tsv"
FILEMOD_2 = "__Universal_Genes.tsv"



# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"WB_GE_Genes_Report.py -h"



STR__header_1 = "# Total number of datasets: "
STR__header_2 = "# Total number of samples: "
STR__header_3 = "# Number of universal genes: "
STR__header_4 = "# Total number of genes: "
STR__header_5 = "# The best performing gene is found in this many datasets: "
STR__header_6 = "# The best performing gene is found in this many samples: "

STR__column_header = "SAMPLE_ID\tTOTAL_DATASETS\tTOTAL_SAMPLES"


STR__metrics = """
                    Total Datasets: {A}
                     Total Samples: {B}

       No. Of Genes In All Samples: {C}
    
                       Total Genes: {D}
                    
    Best Performing Gene(s) Present in:
                          Datasets: {E} ({F}%)
                           Samples: {G} ({H}%)
"""

STR__report_begin = "\nRunning WB_GE_Genes_Overlap_Report..."

STR__report_complete = "\nWB_GE_Genes_Overlap_Report successfully finished."



STR__unexpected_failure = "\nProgram exited with an unexpected error."


# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def WB_GE_Genes_Overlap_Report(input_folder, report_file, universal_file,
            abridge_report):
    """
    Summarize the coverage overlap between the different datasets in
    [input_folder]. That is to say, the degree to which genes appear in most or
    all datasets.
    
    This is useful for determining the feasibility of a meta-analysis combining
    the results of both datasets.
    
    This overlap reporter was developed specifically for Wormbase Gene
    Expression datasets downloaded from the ".cgi" browser tool, where all genes
    have been converted to Wormbase Gene IDs. This script is not applicable to
    datasets where the IDs used have not been normalized.
    
    @input_folder
            (str - dirpath)
            The directory path of the folder containing the raw Wormbase Gene
            Expression file downloads. It is assumed that all valid data files
            within said folder are prefixed by "WBPaper" and have a ".csv" file
            extension, but are formatted as TSVs.
            (Consistent with the download results from October, 2022)
            For the most reliable metrics, it should be ensured that no unique
            biological sample occurs in more than one dataset.
    @report_file
            (str - filepath)
            The filepath of the output file where the overall report will be
            output to. The report file will contain a summary of the results at
            the top of the file, followed by a table showing which genes are
            present in how many datasets, how many samples, which datasets and
            which samples.
    @universal_file
            (str - filepath)
            The filepath of the output file where a shortlist of genes will be
            output to. These are the genes present in all datasets.
    @abridge_report
            (bool)
            Whether or not to abridge the report to just the number of datasets
            and samples a gene is in.
    
    WB_GE_Genes_Overlap_Report(str, str, str) -> int
    """
    PRINT.printP(STR__report_begin)
    
    # Setup - Paths and Lists
    file_list = Get_WB_GE_Files(input_folder)
    
    gene_list = Get_Genes(file_list)
    gene_list = sorted(gene_list)
    dataset_list = Get_Datasets(file_list)
    sample_list = Get_Samples(file_list)
    
    # Setup - Dict
    template = {}
    for sample in sample_list: template[sample] = 0
    template[0] = 0 # 0 is the key for Total count
    
    data = {} # Per sample data
    for gene in gene_list: data[gene] = dict(template)

    #
    template_2 = {}
    for dataset in dataset_list: template_2[dataset] = 0
    template_2[0] = 0 # 0 is the key for Total count
    
    data_2 = {} # Per dataset data
    for gene in gene_list: data_2[gene] = dict(template_2)
    
    
    # Populate with data
    Populate_With_WB_GE_Data(data, data_2, file_list)
    
    # Assess results
    summary_metrics = Get_Summary_Metrics(data, data_2, gene_list, dataset_list,
            sample_list)
    total_datasets = len(dataset_list)
    total_samples = len(sample_list)
    summary_metrics = [total_datasets, total_samples] + summary_metrics
    
    # Write outcomes to file
    Write_Summary_to_File(data, data_2, gene_list, dataset_list, sample_list,
            summary_metrics, report_file, abridge_report)
    Write_Shortlist_to_File(data, gene_list, total_samples, universal_file)
        
    PRINT.printP(STR__report_complete)
    
    # Reporting
    Report_Metrics(summary_metrics)
    
    # Wrap up
    return 0



def Get_WB_GE_Files(input_folder):
    """
    Return a list of the full filepaths of all the files in the [input_folder]
    folder which are prefixed with "WBPaper" and have a ".csv" file extension.
    
    @input_folder
            (str - filepath)
            The directory path of the folder containing the raw Wormbase Gene
            Expression file downloads.
    
    Get_WB_GE_Files(str) -> list<str>
    """
    results = []
    try:
        files = os.listdir(input_folder)
    except:
        return results
    for file_ in files:
        if (file_[:7] == "WBPaper") and (file_[-4:] == ".csv"):
            full_path = input_folder + "/" + file_
            results.append(full_path)
    return results
    
def Get_Genes(file_list):
    """
    Return a list of all the genes contained in the Wormbase Gene Expression
    files in [file_list].
    
    @file_list
            (str - filepath)
            A list containing the filepaths of all the files to be scanned for
            gene names.
    
    Get_Genes(list<str>) -> list<str>
    """
    results = set([])
    #
    f = Table_Reader()
    f.Set_Delimiter("\t")
    #
    for file_ in file_list:
        f.Set_New_Path(file_)
        f.Open()
        f.Read() # Header
        while not f.EOF:
            f.Read()
            gene = f[0]
            if gene not in results: results.add(gene)
        f.Close()
    #
    results = list(results)
    return results
    
def Get_Datasets(file_list):
    """
    Return a list of all the datasets contained in the Wormbase Gene Expression
    files in [file_list].
    
    @file_list
            (str - filepath)
            A list containing the filepaths of all the files.
    
    Get_Datasets(list<str>) -> list<str>
    """
    results = []
    for file_ in file_list:
        file_name = Get_File_Name(file_)
        results.append(file_name)
    return results

def Get_Samples(file_list):
    """
    Return a list of all the sample IDs contained in the Wormbase Gene
    Expression files in [file_list].
    
    @file_list
            (str - filepath)
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



def Populate_With_WB_GE_Data(data, data_2, file_list):
    """
    Populate the given dictionaries [data] and [data_2] with the data in the
    files in [file_list].

    For every gene in [file_list], there will be two corresponding dictionaries:

    In the first dictionary:
        All sample IDs in [file_list] will have an entry, where the key is the
        sample ID, and the value is a 1 or a 0. 1 if the gene is found in that
        sample, and a 0 otherwise.

    In the second dictionary:
        All datasets in [file_list] will have an entry, where the key is the
        dataset file name, and the value is a 1 or a 0. 1 if the gene is found
        in that dataset, and a 0 otherwise.
    
    @data
            (dict<str:dict<str:int>>)
            The "blank" dictionary containing every gene and for every gene, a
            dictionary containing every sample ID corresponding to a value of 0.
    @data_2
            (dict<str:dict<str:int>>)
            The "blank" dictionary containing every gene and for every gene, a
            dictionary containing every dataset name corresponding to a value of
            0.
    @file_list
            (str - filepath)
            A list containing the filepaths of all the data files.
    
    Get_Samples(list<str>) -> None
    """
    f = Table_Reader()
    f.Set_Delimiter("\t")
    #
    for file_ in file_list:
        # Dataset
        dataset = Get_File_Name(file_)
        # Samples
        f.Set_New_Path(file_)
        f.Open()
        ### Header
        f.Read()
        headers = f.Get()
        samples = headers[1:]
        length = len(samples)
        ### Body
        while not f.EOF:
            f.Read()
            gene = f[0]
            # Populate
            data[gene][0] += length # Total count - samples
            for sample in samples:
                data[gene][sample] = 1 # Individual samples
            data_2[gene][0] += 1 # Total count - datasets
            data_2[gene][dataset] = 1 # Individual datasets
        f.Close()



def Get_Summary_Metrics(data, data_2, gene_list, dataset_list, sample_list):
    """
    Return a list of summary metrics for the data, including:
        * The number of genes found universally in all samples
        * The total number of genes
        * For the best performing gene(s), how many datasets were they in
        * For the best performing gene(s), how many samples were they in
    
    @data
            (dict<str:dict<str:int>>)
            The dictionary containing every gene and for every gene, a
            dictionary containing every sample ID corresponding to a value of
            either 1 or 0. A 1 if the gene is in that sample, and a 0 otherwise.
    @data_2
            (dict<str:dict<str:int>>)
            The dictionary containing every gene and for every gene, a
            dictionary containing every dataset name corresponding to a value of
            0. A 1 if the gene is in that sample, and a 0 otherwise.
    @gene_list
            (list<str>)
            A list of all the genes in the data.
    @dataset_list
            (list<str>)
            A list of all the datasets.
    @sample_list
            (list<str>)
            A list of all the samples in the data.
    
    Get_Summary_Metrics(dict<str:dict<str:int>>, dict<str:dict<str:int>>) ->
            [int, int, int, int]
    """
    # Setup
    gene_count = len(gene_list)
    sample_count = len(sample_list)
    perfect_count = 0
    maximum_datasets = 0
    maximum_samples = 0
    # Loop through all genes
    for gene in gene_list:
        gene_in_datasets = data_2[gene][0]
        gene_in_samples = data[gene][0]
        # Adjust max
        if gene_in_datasets > maximum_datasets:
            maximum_datasets = gene_in_datasets
        if gene_in_samples > maximum_samples:
            maximum_samples = gene_in_samples
        # Adjust count
        if data[gene][0] == sample_count: perfect_count += 1
    #
    return [perfect_count, gene_count, maximum_datasets, maximum_samples]



def Write_Summary_to_File(data, data_2, gene_list, dataset_list, sample_list,
            summary_metrics, report_file, abridge_report):
    """
    Write a report of the genes to [report_file].
    The first lines of the output file will begin with "#", and contain the
    following information:
        * The total number of datasets
        * The total number of samples
        * The number of genes in all samples
        * The total number of genes
        * For the best performing gene(s), how many datasets were they in
        * For the best performing gene(s), how many samples were they in
    This will be followed by a table detailing which dataset and which samples
    each gene is found in.
    
    @data
            (dict<str:dict<str:int>>)
            The dictionary containing every gene and for every gene, a
            dictionary containing every sample ID corresponding to a value of
            either 1 or 0. A 1 if the gene is in that sample, and a 0 otherwise.
    @data_2
            (dict<str:dict<str:int>>)
            The dictionary containing every gene and for every gene, a
            dictionary containing every dataset name corresponding to a value of
            0. A 1 if the gene is in that sample, and a 0 otherwise.
    @gene_list
            (list<str>)
            A list of all the genes in the data.
    @dataset_list
            (list<str>)
            A list of all the datasets.
    @sample_list
            (list<str>)
            A list of all the samples in the data.
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The total number of datasets
                * The total number of samples
                * The number of genes in all samples
                * The total number of genes
                * For the best performing gene(s), how many datasets were they
                    in
                * For the best performing gene(s), how many samples were they in
    @report_file
            (str - filepath)
            The filepath to which the results will be written into.
    @abridge_report
            (bool)
            Whether or not to abridge the report to just the number of datasets
            and samples a gene is in.
    
    Write_Summary_to_File(dict<str:dict<str:int>>, dict<str:dict<str:int>>int,
            list<int>, str) -> None
    """
    # Setup
    w = open(report_file, "w")
    # Headers
    w.write(STR__header_1 + str(summary_metrics[0]) + "\n")
    w.write(STR__header_2 + str(summary_metrics[1]) + "\n")
    w.write(STR__header_3 + str(summary_metrics[2]) + "\n")
    w.write(STR__header_4 + str(summary_metrics[3]) + "\n")
    w.write(STR__header_5 + str(summary_metrics[4]) + "\n")
    w.write(STR__header_6 + str(summary_metrics[5]) + "\n")
    # Column headers
    w.write(STR__column_header)
    if not abridge_report:
        for dataset in dataset_list:
            w.write("\t" + dataset)
        for sample in sample_list:
            w.write("\t" + sample)
    w.write("\n")
    # Body
    for gene in gene_list:
        sb = gene
        sb += "\t" + str(data_2[gene][0])
        sb += "\t" + str(data[gene][0])
        if not abridge_report:
            for dataset in dataset_list:
                sb += "\t" + str(data_2[gene][dataset])
            for sample in sample_list:
                sb += "\t" + str(data[gene][sample])
        sb += "\n"
        w.write(sb)
    # Finish
    w.close()

def Write_Shortlist_to_File(data, gene_list, total_samples, universal_file):
    """
    Write a file of the genes which are found in all samples.
    
    @data
            (dict<str:dict<str:int>>)
            The dictionary containing every gene and for every gene, a
            dictionary containing every sample ID corresponding to a value of
            either 1 or 0. A 1 if the gene is in that sample, and a 0 otherwise.
    @gene_list
            (list<str>)
            A list of all the genes in the data.
    @total_samples
            (int)
            The total number of samples.
    @universal_file
            (str - filepath)
            The filepath to which the genes are to be written into.
    
    Write_Shortlist_to_File(dict<str:dict<str:int>>, int, str) -> None
    """
    # Setup
    w = open(universal_file, "w")
    # Body
    for gene in gene_list:
        if data[gene][0] == total_samples:
            w.write(gene + "\n")
    # Finish
    w.close()



def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The total number of datasets
                * The total number of samples
                * The number of genes in all samples
                * The total number of genes
                * For the best performing gene(s), how many datasets were they
                    in
                * For the best performing gene(s), how many samples were they in
    
    Report_Metrics(int, int, int, list<int>, list<int>) -> None
    """
    # Unpacking
    datasets, samples, u_genes, genes, in_datasets, in_samples = summary_metrics
    # Calculations
    percentage_datasets = (in_datasets*100.0)/datasets
    percentage_samples = (in_samples*100.0)/samples
    # Strings
    datasets = str(datasets)
    samples = str(samples)
    u_genes = str(u_genes)
    genes = str(genes)
    in_datasets = str(in_datasets)
    in_samples = str(in_samples)
    percentage_datasets = str(percentage_datasets)
    percentage_samples = str(percentage_samples)
    # Pad percentages
    percentage_datasets = Trim_Percentage_Str(percentage_datasets, 2)
    percentage_samples = Trim_Percentage_Str(percentage_samples, 2)
    # Pad ints
    max_size = Get_Max_Len([datasets, samples, u_genes, genes, in_datasets,
            in_samples])
    datasets = Pad_Str(datasets, max_size, " ", 0)
    samples = Pad_Str(samples, max_size, " ", 0)
    u_genes = Pad_Str(u_genes, max_size, " ", 0)
    genes = Pad_Str(genes, max_size, " ", 0)
    in_datasets = Pad_Str(in_datasets, max_size, " ", 0)
    in_samples = Pad_Str(in_samples, max_size, " ", 0)
    # Print
    PRINT.printM(STR__metrics.format(A = datasets, B = samples, C = u_genes,
            D = genes, E = in_datasets, F = percentage_datasets, G = in_samples,
            H = percentage_samples))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__WB_GE_Genes_Overlap_Report(
            raw_command_line_input):
    """
    Parse the command line input and call the WB_GE_Genes_Overlap_Report
    function with appropriate arguments if the command line input is valid.
    """
    PRINT.printP(STR__parsing_args)
    
    # Safe exit
    return 0



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__WB_GE_Genes_Overlap_Report(sys.argv)
