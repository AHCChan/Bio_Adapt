HELP_DOC = """
WORMBASE GENE EXPRESSION GENES OVERLAP REPORT
(version 2.0)
by Angelo Chan
This is a program for taking a folder of Wormbase Gene Expression files and
generating a report on the genes listed within those files.

The report will state the number of datasets, the number of samples, the number
of genes found in all samples, some other information, and for each gene, state
the number of datasets and samples which contain said gene.

NOTE: In some datasets, a gene is not detected in some samples, but the gene is
still included in the results. If the results are absolute values like copy
number or signal strength, then such samples will simply have a value of 0.
However, it is common practice to log-transform gene expression data for a
variety of reasons. However, 0 cannot be log-transformed. A common solution to
this problem is to set the "log-tansformed value" to a very low number, such
as -33.2192809488736. (Under a log2 transform system, this would be the result
if the original, non-log2-trasnformed value was 0.0000000001.) Such an approach
may work for some analyses, but greatly distort the results for others. More
importantly, genes with such values need to be excluded from certain analyses,
despite possibly being "present" in all datasets.

Optionally, commonality scores for each dataset can also be generated.
Commonality scores are calculated as follows:
        
  * For each sample, the sample's commonality score is the sum of the
    commonality scores of all the genes which are present in that sample.
    
  * The commonality score for a gene is the number of samples which that gene is
    present in.



USAGE:
    
    python27 WB_GE_Genes_Overlap_Report.py <input_folder> [-o <report_file>
            <universal_file>] [-a <abridge>] [-c <commonality_file>]



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
    
    commonality_file
        
        The filepath of the output file where the commonality scores are output
        to.



EXAMPLE:
    
    python27 WB_GE_Genes_Overlap_Report.py Path/WB_GE_Results -o
            Genes_Summary.tsv Genes_Universal.tsv -a Y



USAGE:
    
    python27 WB_GE_Genes_Overlap_Report.py <input_folder> [-o <report_file>
            <universal_file>] [-a <abridge>]
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



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__Abridge = False



# Imported Modules #############################################################

import sys
import os



import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"WB_GE_Genes_Report.py -h"



STR__error_no_WB_GE = "ERROR: No Wormbase gene expression data files detected "\
" in:\n\t{f}"



STR__header_1 = "# Total number of datasets: "
STR__header_2 = "# Total number of samples: "
STR__header_3 = "# Total number of genes: "
STR__header_4 = "# Number of universal genes: "
STR__header_5 = "# The best performing gene is found in this many datasets: "
STR__header_6 = "# The best performing gene is found in this many samples: "

STR__column_header = "GENE_ID\tTOTAL_DATASETS\tTOTAL_SAMPLES"

STR__commonality_header = "DATASET\tSAMPLE_SCORE\tDATASET_SCORE"



STR__metrics = """
                    Total Datasets: {A}
                     Total Samples: {B}
                       Total Genes: {C}
                       
       No. Of Genes In All Samples: {D} ({E}%)
                    
    Best Performing Gene(s) Present in:
                          Datasets: {F} ({G}%)
                           Samples: {H} ({I}%)
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
            abridge_report, commonality_file):
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
    @commonality_file
            (None/str - filepath)
            The filepath of the output file where the commonality scores are
            output to.
    
    WB_GE_Genes_Overlap_Report(str, str, str, bool, str) -> int
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
    total_genes = len(gene_list)
    summary_metrics = ([total_datasets, total_samples, total_genes] +
            summary_metrics)
    
    # Write outcomes to file
    Write_Summary_to_File(data, data_2, gene_list, dataset_list, sample_list,
            summary_metrics, report_file, abridge_report)
    Write_Shortlist_to_File(data, gene_list, total_samples, universal_file)
    if commonality_file:
        Write_Commonality_to_File(data, data_2, gene_list, dataset_list,
                commonality_file)
    
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
        if Validate_WB_CE_File(file_):
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
    
    NOTE: In some datasets, a gene is not detected in some samples, but the gene
    is still included in the results. If the results are absolute values like
    copy number or signal strength, then such samples will simply have a value
    of 0. However, it is common practice to log-transform gene expression data
    for a variety of reasons. However, 0 cannot be log-transformed. A common
    solution to this problem is to set the "log-tansformed value" to a very low
    number, such as -33.2192809488736. (Under a log2 transform system, this
    would be the result if the original, non-log2-trasnformed value was
    0.0000000001.) Such an approach may work for some analyses, but greatly
    distort the results for others. More importantly, genes with such values
    need to be excluded from certain analyses, despite possibly being "present"
    in all datasets.
    
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
    
    Populate_With_WB_GE_Data(dict, dict, list<str>) -> None
    """
    f = Table_Reader()
    f.Set_Delimiter("\t")
    invalid_1 = "-33"
    invalid_2 = "\\N"
    invalid_3 = "Inf"
    invalids = [invalid_1, invalid_2, invalid_3]
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
        range_ = range(length)
        ### Body
        while not f.EOF:
            f.Read()
            gene = f[0]
            # Populate
            flag_any = 0
            for i in range_:
                sample = samples[i]
                value = f[i+1]
                if value[:3] not in invalids:
                    flag_any += 1
                    data[gene][sample] = 1 # Individual samples
            if flag_any:
                data[gene][0] += flag_any # Total count - samples
                data_2[gene][0] += 1 # Total count - datasets
                data_2[gene][dataset] = 1 # Individual datasets
        f.Close()



def Write_Commonality_to_File(data, data_2, gene_list, dataset_list,
            commonality_file):
    """
    Generate a commonality score for each dataset, and write those scores to
    [commonality_file].

    For every gene in, there will be two corresponding dictionaries:

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
    @commonality_file
            (None/str - filepath)
            The filepath of the output file where the commonality scores are
            output to.
    
    Write_Commonality_to_File(dict, dict, list<str>, list<str>, str) -> None
    """
    # Setup
    results = {}
    # Iterate over data
    for dataset in dataset_list:
        score_1 = 0
        score_2 = 0
        for gene in gene_list:
            if data_2[gene][dataset] == 1:
                score_1 += data[gene][0]
                score_2 += data_2[gene][0]
        results[dataset] = [score_1, score_2]
    # Sort
    order = sorted(results.keys(), key=lambda kv: results[kv][1])
    # Write - initial
    o = open(commonality_file, "w")
    o.write(STR__commonality_header)
    o.write("\n")
    # Write - body
    for dataset in order:
        score_1, score_2 = results[dataset]
        o.write(dataset + "\t" + str(score_1) + "\t" + str(score_2) + "\n")
    # Close
    o.close()
    return



def Get_Summary_Metrics(data, data_2, gene_list, dataset_list, sample_list):
    """
    Return a list of summary metrics for the data, including:
        * The number of genes found universally in all samples
        * The percentage of genes found universally in all samples (string)
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
            [int, str, int, int]
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
    # Perfect string
    percentage = (perfect_count*100.0)/gene_count
    percentage = str(percentage)
    percentage = Trim_Percentage_Str(percentage, 2)
    # Return
    return [perfect_count, percentage, maximum_datasets, maximum_samples]



def Write_Summary_to_File(data, data_2, gene_list, dataset_list, sample_list,
            summary_metrics, report_file, abridge_report):
    """
    Write a report of the genes to [report_file].
    The first lines of the output file will begin with "#", and contain the
    following information:
        * The total number of datasets
        * The total number of samples
        * The total number of genes
        * The number of genes in all samples
        * The percentage of genes found universally in all samples (string)
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
            (list<int/str>)
            A list of summary metrics for the data, including:
                * The total number of datasets
                * The total number of samples
                * The total number of genes
                * The number of genes in all samples
                * The percentage of genes found universally in all samples
                        (string)
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
    w.write(STR__header_4 + str(summary_metrics[3]) + " (" +
            summary_metrics[4] + "%)\n")
    w.write(STR__header_5 + str(summary_metrics[5]) + "\n")
    w.write(STR__header_6 + str(summary_metrics[6]) + "\n")
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
                * The total number of genes
                * The number of genes in all samples
                * The percentage of genes found universally in all samples
                        (string)
                * For the best performing gene(s), how many datasets were they
                    in
                * For the best performing gene(s), how many samples were they in
    
    Report_Metrics([int, int, int, int, str, list<int>, list<int>]) -> None
    """
    # Unpacking
    datasets = summary_metrics[0]
    samples = summary_metrics[1]
    genes = summary_metrics[2]
    u_genes = summary_metrics[3]
    u_percent = summary_metrics[4]
    in_datasets = summary_metrics[5]
    in_samples = summary_metrics[6]
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
    PRINT.printM(STR__metrics.format(A = datasets, B = samples, C = genes,
            D = u_genes, E = u_percent, F = in_datasets,
            G = percentage_datasets, H = in_samples, I = percentage_samples))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__WB_GE_Genes_Overlap_Report(
            raw_command_line_input):
    """
    Parse the command line input and call the WB_GE_Genes_Overlap_Report
    function with appropriate arguments if the command line input is valid.
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
    path_in = inputs.pop(0)
    valid = Validate_WB_GE_Data_Folder(path_in)
    if valid == 1:
        PRINT.printE(STR__error_no_WB_GE.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    elif valid == 2:
        PRINT.printE(STR__IO_error_read_folder)
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    abridge = DEFAULT__Abridge
    path_out_report = ""
    path_out_universal = ""
    path_out_commonality = None
    
    # Initial validation
    while inputs:
        arg = inputs.pop(0)
        try: # Following arguments
            if arg in ["-o"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
            elif arg in ["-a", "-c"]:
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
            path_out_report = arg2
            path_out_universal = arg3
        elif arg == "-a":
            abridge = Validate_Bool(arg2)
            if abridge == None:
                printE(STR__invalid_bool)
                return 1
        else: # arg == "-c"
            path_out_commonality = arg2
    
    # Automated output path generation
    if not path_out_report: path_out = path_in_folder + FILEMOD_1
    if not path_out_universal: path_out = path_in_folder + FILEMOD_2
    
    # Validate output path
    valid_out = Validate_Write_Path__FILE(path_out_report)
    if valid_out == 2: return 0
    if valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    valid_out = Validate_Write_Path__FILE(path_out_universal)
    if valid_out == 2: return 0
    if valid_out == 3:
        PRINT.printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        PRINT.printE(STR__IO_error_write_unable)
        return 1
    if path_out_commonality:
        valid_out = Validate_Write_Path__FILE(path_out_commonality)
        if valid_out == 2: return 0
        if valid_out == 3:
            PRINT.printE(STR__IO_error_write_forbid)
            return 1
        if valid_out == 4:
            PRINT.printE(STR__IO_error_write_unable)
            return 1
    
    # Run program
    exit_state = WB_GE_Genes_Overlap_Report(path_in, path_out_report,
            path_out_universal, abridge, path_out_commonality)
    
    # Safe exit
    if exit_state == 0: return 0
    else:
        PRINT.printE(STR__use_help)
        return 1
    
def Validate_WB_GE_Data_Folder(dirpath):
    """
    Validates the dirpath of the input folder as containing raw downloaded
    Wormbase Gene Expression data files.
    
    A raw downloaded Wormbase Gene Expression data files is defined as a file
    suffixed by "WBPaper" and has a ".csv" extension.
    
    Return 0 if the dirpath is valid and contains at least 1 valid file.
    Return 1 if the dirpath is valid but contains no valid files.
    Return 2 if the dirpath is invalid.
    
    Validate_Read_Path(str) -> int
    """
    try:
        os.listdir(dirpath)
        files = os.listdir(dirpath)
        files = [f for f in files if Validate_WB_CE_File(f)]
        if len(files) > 0: return 0
        return 1
    except:
        return 2

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
    
def Validate_WB_CE_File(filepath):
    """
    Validates a file as a raw Wormbase gene expression data file.
    
    A raw downloaded Wormbase Gene Expression data files is defined as a file
    suffixed by "WBPaper" and has a ".csv" extension.
    
    Return True if the file is such a file.
    Return False otherwise.
    """
    if filepath[:7] == "WBPaper" and filepath[-4:] == ".csv": return True
    return False



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__WB_GE_Genes_Overlap_Report(sys.argv)
