HELP_DOC = """
BED POSTMERGE UNCOLLAPSE COMBINE REPLICATES
(version 1.0)
by Angelo Chan

Designed to run on the output of BED__Postmerge_Uncollapse.py. (Bio_Files_Tools)

This is a program for summarizing the count and absence of genetic elements in
multiple samples, aggregated by experimental group, with different replicates
of the same experimental group being combined.



USAGE:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py <input_path>
            <grouping_file_path> [-m <output_path_main_p> <output_path_main_a>
            <output_path_main_t> <output_path_main_m>] [-a <output_path_all_p>
            <output_path_all_a> <output_path_all_t> <output_path_all_m>]
            [-u <output_path_uniques_p> <output_path_uniques_a>
            <output_path_uniques_t> <output_path_uniques_m>]



MANDATORY:
    
    input_path
        
        The filepath of the input BED file. This input file should be the
        _SUMMARY file produced by BED__Postmerge_Uncollapse.py.

        This summary file will be a TSV with the following columns:
            1)  The first 3 columns will be the genomic coordinates of the
                merged regions. (1 region per row)
            2)  For each file which went into the pipeline, the number of loci
                from that file which went into merged region.
    
    grouping_file_path
        
        The two-column TSV denoting which columns belong to which experimental
        groups. The first column will contain the name of the experimental
        group, while second column contains comma-separated column numbers.
        
        The column numbers follow a 1-index system, and ignore the first 3
        columns of the BED file. (which contain the genomic coordinates) Thus,
        the very first data column is column "1".
        
        Ex.:
            Wildtype_ctrl   1,2,3,4
            Wildtype_drug   5,6,7
            Mutated_ctrl    8,9,10
            Mutated_drug    11,12,13,14

OPTIONAL:[-m <output_path_main>]
            [-u <output_path_uniques>] [-s <output_path_summary>]
    
    output_path_main_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output files which contain the results of combining
        all replicates from each experimental group.
    
    output_path_all_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output files which contain the genetic elements
        found in all experimental groups.
    
    output_path_uniques_{*}
        
        (DEFAULT path generation available)
        
        The filepath of the output folders which contain the group-specific
        files, for genetic elements found exclusively in one experimental
        group and none of the others. (None of the columns specified in the
        grouping file)
    
    (SUFFIXES)
        
        (present), (all), (total), (mean)
            
            [*_p] contains either a 1 or a 0. 1 if at least 1 replicate contains
                    at least 1 count. 0 otherwise.
            [*_a] contains either a 1 or a 0. 1 if all replicates contain at
                    least 1 count. 0 otherwise.
            [*_t] contains the total counts across all replicates.
            [*_m] contains the average number of counts per replicate.



EXAMPLES:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py
            data\summary.bed data\list_of_groups.tsv
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py
            data\summary.bed data\list_of_groups.tsv -m present.bed
            present_all_reps.bed main_total_counts.bed main_average_counts.bed
            -a present_all_groups.bed present_all_groups_all_reps.bed
            present_all_totals.bed present_all_mean.bed -u uniquely_present
            uniquely_present_all_reps uniquely_present_totals
            unique_present_averages

USAGE:
    
    python27 BED__Postmerge_Uncollapse_Combine_Replicates.py <input_path>
            <grouping_file_path> [-m <output_path_main_p> <output_path_main_a>
            <output_path_main_t> <output_path_main_m>] [-a <output_path_all_p>
            <output_path_all_a> <output_path_all_t> <output_path_all_m>]
            [-u <output_path_uniques_p> <output_path_uniques_a>
            <output_path_uniques_t> <output_path_uniques_m>]
"""

NAME = "BED__Postmerge_Uncollapse_Combine_Replicates.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = False # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD__MAIN = "__MAIN"
FILEMOD__ALL = "__ALLGROUPS"
FILEMOD__UNIQUE = "__UNIQUE"

FILEMOD__P = "_PRESENT"
FILEMOD__A = "_ALLREPS"
FILEMOD__T = "_TOTALS"
FILEMOD__M = "_AVERAGES"



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *

import random as Random



# Strings ######################################################################

STR__BED_header = "chr\tstart\tend"



STR__use_help = "\nUse the -h option for help:\n\t python "\
"BED__Postmerge_Uncollapse_Combine_Replicates.py -h"



STR__invalid_column_number = """
ERROR: Invalid column number: {s}
Raw text: {t}
Please specify positive integers which correspond to column numbers which exist
in the main input file."""



STR__metrics = """
                       Rows in input file: {A}

    Universally present in all groups:

                In at least one replicate: {B}
                        In all replicates: {C}

    Uniquely present in one group:
        
                In at least one replicate: {D}
                        In all replicates: {E}
"""



STR__combine_begin = "\nRunning Combine_Replicates..."

STR__combine_complete = "\nCombine_Replicates successfully finished."



STR__overwrite_confirm_2 = """
Files may exist in destination folder:
    {f}
Do you wish to overwrite them in the event of a naming clash? (y/n): """



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Combine_Replicates(path_in, path_groups, paths_out):
    """
    Generate a series of FASTA files each containing a synthetic chromosome.
    
    @path_in
            (str - filepath)
            The filepath of the input BED file produced by
            BED__Postmerge_Uncollapse.py.
    @path_groups
            (str - filepath)
            The filepath of the TSV file specifying which columns belong to
            which experimental groups.
    @paths_out
            (list<list<str>> - 8x filepaths, 4x dirpaths)
            Nested lists, containing the output filepaths and dirpaths.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    Combine_Replicates(str, str, [[str, str, str, str], [str, str, str, str],
            [str, str, str, str]]) -> int
    """
    PRINT.printP(STR__combine_begin)
    
    # Setup reporting
    count_rows = 0
    count_universal_one = 0
    count_universal_all = 0
    count_unique_one = 0
    count_unique_all = 0
    
    # Get groups
    processed_groups = Process_Groups(path_groups)
    if processed_groups == None: return 1
    groups_list, groups_dict, empty_dicts_dict, lengths_dict = processed_groups
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_in)
    f.Set_Delimiter("\t")
    f.Close()
            
    header_str = STR__BED_header
    for group in groups_list:
        header_str += "\t" + group
    header_str += "\n"
    outputs = Setup_Outputs(paths_out, groups_dict, header_str)
    
    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        values = f.Get()
        count_rows += 1
        # Flags, totals and tallies, for ALL
        flag_all_p = True
        flag_all_a = True
        dict_present = {}
        dict_all = {}
        dict_total = {}
        dict_avgs = {}
        groups_present = 0
        last_group_present = ""
        last_str = ""
        # Per group
        for group in groups_list:
            # Flags, totals and tallies, per group
            flag_p = False
            flag_a = True
            total = 0
            col_nos = groups_dict[group]
            # Sb
            sb = ""
            # Go over all replicates
            for col in col_nos:
                string = values[col]
                # Stringbuilder
                sb += "\t" + string
                # Counts and flags
                number = int(string)
                if number:
                    flag_p = True
                    total += number
                else:
                    flag_a = False                
            # For individuals
            dict_total[group] = total
            avg = total/lengths_dict[group]
            dict_avgs[group] = avg
            # Update dicts
            if flag_p:
                dict_present[group] = 1
                #
                groups_present += 1
                last_group_present = group
                last_str = sb
            else:
                flag_all_p = False
                #
                dict_present[group] = 0
            if flag_a:
                dict_all[group] = 1
            else:
                flag_all_a = False
                #
                dict_all[group] = 0
        # Stringbuilder
        coords = values[:3]
        sb = "\t".join(coords)
        sb_present = sb
        sb_all = sb
        sb_total = sb
        sb_avg = sb
        for group in groups_list:
            present = dict_present[group]
            sb_present += "\t" + str(present)
            all_ = dict_all[group]
            sb_all += "\t" + str(all_)
            total = dict_total[group]
            sb_total += "\t" + str(total)
            avg = dict_avgs[group]
            sb_avg += "\t" + str(avg)
        sb_present += "\n"
        sb_all += "\n"
        sb_total += "\n"
        sb_avg += "\n"
        # Write
        outputs[0][0].write(sb_present)
        outputs[0][1].write(sb_all)
        outputs[0][2].write(sb_total)
        outputs[0][3].write(sb_avg)
        if flag_all_p: # Main
            count_universal_one += 1
            outputs[1][0].write(sb_present)
            outputs[1][2].write(sb_total)
            outputs[1][3].write(sb_avg)
        if flag_all_a: # All
            count_universal_all += 1
            outputs[1][1].write(sb_all)
        if groups_present == 1: # Unique
            count_unique_one += 1
            group = last_group_present # The last one to be present
            # Reiterating the data
            sb_unique = sb + last_str + "\n"
            outputs[2][0][group].write(sb_unique)
            if dict_all[group] == 1:
                count_unique_all += 1
                outputs[2][1][group].write(sb_unique)
            # Individual totals and averages
            total = dict_total[group]
            sb_unique_total = sb + "\t" + str(total) + "\n"
            avg = dict_avgs[group]
            sb_unique_avg = sb + "\t" + str(avg) + "\n"
            outputs[2][2][group].write(sb_unique_total)
            outputs[2][3][group].write(sb_unique_avg)
    
    # Finish
    f.Close()
    for i in outputs[:2]:
        for j in i:
            j.close()
    for j in outputs[2]:
        for group in groups_list:
            j[group].close()
    PRINT.printP(STR__combine_complete)
    
    # Reporting
    Report_Metrics([count_rows, count_universal_one, count_universal_all,
            count_unique_one, count_unique_all])
    
    # Wrap up
    return 0



def Process_Groups(path_groups):
    """
    Process a "groups" file to obtain the relevant data for the analysis.
    Return a list of all the groups, a dictionary of all the column numbers
    (computational, 0-indexed column numbers) for each group, and a dictionary
    which contains, for each group, a dictionary with all the column numbers as
    keys, and 0 as values. Also return a dictionary giving the number of columns
    for each group.
    
    Return None if there is an invalid column number.

    Ex.
        A groups file will the following:
            Control 1,2,5
            Treatment   3,4,6
        Will produce:
            ["Control", "Treatment"]
            {"Control": [3,4,7], "Treatment": [5,6,8]}
            {"Control": {3:0, 4:0, 7:0}, "Treatment": {5:0, 6:0, 8:0}}
            {"Control": 3.0, "Treatment": 3.0}

    @paths_groups
            (str - filepath)
            The filepath of the TSV file specifying which columns belong to
            which experimental groups.
    
    Process_Groups(str) -> [[list<str>, dict<str:<list<int>>>,
            dict<str<dict<int:int>>>, dict<str:float>]
    Process_Groups(str) -> None
    """
    # Output setup
    groups_list = []
    groups_dict = {}
    empty_dicts_dict = {}
    lengths_dict = {}
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_groups)
    f.Set_Delimiter("\t")

    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        # Process raw
        name = f[0]
        cols_raw = f[1]
        cols_str = cols_raw.split(",")
        cols_str = [s.strip(" ") for s in cols_str]
        cols = []
        for s in cols_str:
            try:
                i = int(s)
            except:
                PRINT.printE(STR__invalid_column_number.format(s=i, t=cols_raw))
                return None
            cols.append(i)
        cols = [i+2 for i in cols]
        # Output
        lengths_dict[name] = float(len(cols))
        groups_list.append(name)
        groups_dict[name] = cols
        temp = {}
        for i in cols:
            temp[i] = 0
        empty_dicts_dict[name] = temp
    f.Close()
    
    # Return
    return [groups_list, groups_dict, empty_dicts_dict, lengths_dict]

def Setup_Outputs(paths_out, groups_dict, header_str):
    """
    Setup all the necessary output files.
    
    TODO
    """
    result = []
    # MAIN and ALL
    for i in paths_out[:2]:
        temp = []
        for j in i:
            file_writer = open(j, "w")
            file_writer.write(header_str)
            temp.append(file_writer)
        result.append(temp)
    # UNIQUE
    uniques = []
    for j in paths_out[2]:
        temp = {}
        for group in groups_dict:
            file_path = j + "\\" + group + ".bed"
            file_writer = open(file_path, "w")
            temp[group] = file_writer
        uniques.append(temp)
    result.append(uniques)
    # Return
    return result

def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
            TODO
    
    Report_Metrics() -> None
    """
    pass



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Combine_Uncollapsed(raw_command_line_input):
    """
    Parse the command line input and call the Generate_Synthetic_Genome function
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
    
    # Validate mandatory inputs (file I/O)
    path_in = inputs.pop(0)
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    path_groups = inputs.pop(0)
    valid = Validate_Read_Path(path_groups)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_groups))
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs (file contents)
    groups = Process_Groups(path_groups)
    if groups == None:
        PRINT.printE(STR__use_help)
        return 1
    names, void, void, void = groups
    
    # Set up rest of the parsing
    paths_out = []
    for i in [FILEMOD__MAIN, FILEMOD__ALL]:
        temp = []
        for j in [FILEMOD__P, FILEMOD__A, FILEMOD__T, FILEMOD__M]:
            mod = i + j
            path_temp = Generate_Default_Output_File_Path_From_File(path_in,
                    mod, True)
            temp.append(path_temp)
        paths_out.append(temp)
    for i in [FILEMOD__UNIQUE]:
        temp = []
        base_path = Generate_Default_Output_Folder_Path(path_in)
        for j in [FILEMOD__P, FILEMOD__A, FILEMOD__T, FILEMOD__M]:
            path_temp = base_path + i + j
            temp.append(path_temp)
        paths_out.append(temp)
    Generate_Default_Output_Folder_Path
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-m", "-a", "-u"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
                arg4 = inputs.pop(0)
                arg5 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        if arg == "-m":
            paths_out[0] = [arg2, arg3, arg4, arg5]
        elif arg == "-a":
            paths_out[1] = [arg2, arg3, arg4, arg5]
        else: # arg == "-u"
            paths_out[2] = [arg2, arg3, arg4, arg5]
    
    # Validate output paths
    for i in paths_out[:2]:
        for j in i:
            valid_out = Validate_Write_Path__FILE(j)
            if valid_out == 2: return 0
            if valid_out == 3:
                PRINT.printE(STR__IO_error_write_forbid)
                return 1
            if valid_out == 4:
                PRINT.printE(STR__In_error_write_unable)
                return 1
    for i in [paths_out[2]]:
        for j in i:
            valid_folder = Validate_Write_Path__FOLDER(j)
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
            for name in names:
                file_path = j + "\\" + name + ".bed"
                valid_out = Validate_Write_Path__FILE(file_path)
                if valid_out == 2: return 0
                if valid_out == 3:
                    PRINT.printE(STR__IO_error_write_forbid)
                    return 1
                if valid_out == 4:
                    PRINT.printE(STR__In_error_write_unable)
                    return 1
    
    # Run program
    exit_state = Combine_Replicates(path_in, path_groups, paths_out)
    
    # Exit
    if exit_state == 0: return 0
    else: return 1


    
def Validate_Write_Path__FILE(filepath):
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
    exit_code = Parse_Command_Line_Input__Combine_Uncollapsed(sys.argv)
