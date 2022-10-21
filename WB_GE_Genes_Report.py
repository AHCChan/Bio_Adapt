HELP_DOC = """
WORMBASE GENE EXPRESSION GENES REPORT
(version 1.0)
by Angelo Chan

This is a program for taking a folder of Wormbase Gene Expression files and
generating a report on the genes listed within those files.

The report will state the number of datasets, the number of samples, the number
of genes found in all samples, some other information, and for each gene, state
the number of datasets and samples which contain said gene.
"""

NAME = "WB_GE_Genes_Report.py"



# Configurations ###############################################################



# Minor Configurations #########################################################



# Imported Modules #############################################################



# Strings ######################################################################



# Lists ########################################################################



# Dictionaries #################################################################



# Apply Globals ################################################################



# Functions ####################################################################

def WB_GE_Genes_Report():
    """
    """
    pass



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__WB_GE_Genes_Report(raw_command_line_input):
    """
    Parse the command line input and call the WB_GE_Genes_Report function with
    appropriate arguments if the command line input is valid.
    """
    PRINT.printP(STR__parsing_args)
    
    # Safe exit
    return 0



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__WB_GE_Genes_Report(sys.argv)
