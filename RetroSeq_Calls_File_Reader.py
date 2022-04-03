"""
RETROSEQ CALLS FILE READER
(version 1.0)
by Angelo Chan

This module contains a Class capable of reading and interpretting a RetroSeq
Calls VCF file.
"""

# Imported Modules #############################################################

from File_Reader import *
from Table_File_Reader import *



# Lists ########################################################################



# Dictionaries #################################################################



# Classes ######################################################################

class RetroSeq_Calls_Reader(Table_Reader):
    """
    The RetroSeq Calls Reader is a file reader designed specifically to work
    with RetroSeq Calls files, which are VCF files plus a few unique
    stipulations.
    
    VCF files themselves are a type of TSV file.
    
    Aside from reading and parsing the data, row by row, the RetroSeq Calls
    reader was specifically made to look at data across multiple rows which
    fall within the same "region" and pick the best "call".

    Designed for the following use:
    
    f = RetroSeq_Calls_Reader()
    f.Set_New_Path("F:/Filepath.vcf")
    f.Open()
    
    x = Other_File_Reader()
    # Intialize other file reader
    
    while not x.EOF:
        coordinates = x.Get_Coordinates() # Example function
        f.Read_Until(coordinates)
        best = f.Get_Best()
        # Your code - You may access buffered elements in best
    
    x.Close()
    f.Close()
    """
    
    # Minor Configurations #####################################################
    
    empty_element = []
    
    # Minor Configurations #####################################################
    
    _CONFIG__print_errors = True
    _CONFIG__print_progress = False
    _CONFIG__print_metrics = True
    
    
    
    # Strings ##################################################################
    
    _MSG__object_type = "RetroSeq Calls File Reader"
    _MSG__units_of_measure = "Rows"
    
    _MSG__faulty_line = "Faulty line of data:\n\t{S}"
    
    

    # Constructor & Destructor #################################################
    
    def __init__(self, file_path="", auto_open=False, delimiter="",
                enclosers=[], header_params=[], keep_enclosers=True):
        """
        Creates a File Reader object. The filepath will be tested if a
        filepath is supplied.
        """
        Table_Reader.__init__(self, file_path, False, "\t", [], [], True)
        self.current_raw = self.next_raw = ""
        self.header_text = ""
        self.buffer = []
    
    
    
    # Property Methods #########################################################
    
    # File I/O Methods #########################################################
            
    # File Reading Methods #####################################################
    
    def Get_Buffer(self):
        """
        Return a copy of the current element.
        """
        result = []
        for list_ in self.buffer:
            result.append(list(list_))
        return result
    
    def Read_Until(self, coords):
        """
        Read until either the end of the file, or until the next entry in the
        file is outside the specified range.

        @coords
                ([str, int, int])
                The coordinates to read until the end of. The values are,
                respectively, the chromosome name, the start of the region of
                interest (bp) and the end of the region of interest (bp).
        """
        self._clear_buffer(coords)
        flag = False
        while self._behind_range(coords) and not self.EOF:
            self.Read(1)
        while self._within_range(coords) and not self.EOF:
            self.Read(1)
            current = self.Get_Current()
            self.buffer.append(current)
    
    def Get_Best(self):
        """
        Return the data entry or entries with the highest read count, in the
        current buffer of data entries.
        
        If there are two or more data entries with the same and highest read
        counts, return them all.
        
        If the buffer is empty, return an empty list.
        """
        if not self.buffer: return [] # Empty buffer
        best = [self.buffer[0]] # First element in buffer
        for entry in self.buffer[1:]: # Look at the rest
            if entry[2] > best[0][2]: # New winner
                best = [entry]
            elif entry[2] == best[0][2]: # Tied
                best.append(entry)
        return best
    
    def _clear_buffer(self, coords):
        """
        Removes all data in the buffer which does not fall within the given
        coordinates.
        """
        flag = True
        while flag and self.buffer:
            if coords[0] != self.buffer[0][0] or coords[1] > self.buffer[0][1]:
                self.buffer.pop(0)
            else:
                flag = False
    
    def _behind_range(self, coords):
        """
        Check to see if the "next" entry is behind the [coords] specified.
        """
        future = self.next_element
        if not future: return False
        if future[0] != coords[0]: return True
        if future[1] < coords[1]: return True
        return False
    
    def _within_range(self, coords):
        """
        Check to see if the "next" entry is within the [coords] specified.
        """
        future = self.next_element
        if not future: return False
        if future[0] != coords[0]: return False
        if future[1] < coords[1]: return False
        if future[1] > coords[2]: return False
        return True
    
    def _get_next_element(self):
        """
        Read in the next row and process it.
        
        Return an empty string if the end of the file has been reached.
        """
        raw = self.file.readline()
        data = self._process_line(raw)
        return data
    
    def _process_line(self, line):
        """
        Process a line of data.
        
        Returns a list containing the chromsome name, a basepair coordinate,
        the number of reads which support the called TE, and the called TE.
        
        _proces_line(str) -> [str, int, int, str]
        """
        if not line: return []
        values = line.split(self.delimiter)
        if not values: return []
        if len(values) != 10:
            self.printE(self._MSG__faulty_line.format(S = str([line])[1:-1]))
            return []
        chromosome = values[0]
        coord = int(values[1])
        count = int(values[5])
        string = values[7]
        substrings = string.split(",")
        subsubstrings = substrings[0].split("=")
        called = subsubstrings[1]
        return [chromosome, coord, count, called]


