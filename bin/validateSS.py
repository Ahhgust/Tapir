#!/usr/bin/env python3
# 
# Author: August Woerner
# This script does a rudimentary sanity test of a sample sheet.
# It checks for duplicate sample names, nonalphanumeric characters in the sample names (a requirement of bcl2fastq2)
# as well as blank sample names.
# Note that the name is drawn from Sample_Name ; if that's blank, Sample_ID is used
# (this is also what bcl2fastq2 does)
# This script may or may not be appropriate if bclconvert is used instead.

# Benjamin Crysup: added whodunargs support

import os
import sys

import whodunargs


def sanitizeString(s):
    """
    strips out non-alphanumeric symbols (and _, -)  from a string s (input), and returns it
    """
    return "".join(c for c in s if c.isalnum() or c == "_" or c == '-')


def parseSamplesheet(fh):
    """
    Takes in a sample sheet (SampleSheet.csv is most common)
    and does sanity testing on it. 
    Checks for nonalphanumeric (- also allowed) characters in the sample names
    duplicate sample names
    and returns two dictionaries;
    the first emulates stdout, the second emulates stderr
    key, value => name of error, what the error is)
    Note that all keys are strings, and all values are lists of strings (for multiple types of the same class of error)
    """
    getData=False
    getHeader=False
    idCol = -1
    nameCol=-1
    duptest= set()
    experiment=""
    
    # returned
    errors={}
    out={}
    
    
    for line in fh:
        line=line.rstrip()
        if line.startswith("[Data]"):
            getHeader=True
        elif not getHeader and experiment=="":
            if line.startswith("Experiment Name"):
                # this is a CSV element. grab the first non-empty element.
                s = line.split(",")[1:]
                for elem in s:
                    if elem != "":
                        experiment=sanitizeString(elem)
                        break
            
        elif getHeader:
            getData=True
            getHeader=False
            sp = line.split(",")
            i=0
            for colid in sp:
                if colid== "Sample_ID":
                    idCol=i
                elif colid=="Sample_Name":
                    nameCol=i
                i += 1

            if idCol< 0 and nameCol < 0:
                errors["Corrupt"] = ["Error: Your sample sheet appears to be corrupt. Are you sure you grabbed the right file?"]
                return out, errors
                 

                
        elif getData:
            sp = line.split(",")
            sid=sname=""
            if idCol >= 0:
                sid = sp[idCol]
            
            if nameCol>=0:
                sname=sp[nameCol]
            
            # if no sample name is given, the sample name is the sid
            name = sid
            if sname != "":
                name = sname
            
            if name == "":
                if "Sample" not in errors:
                    errors["Sample"] = []
                errors["Sample"].append("Error: Blank sample name detected.")
                
            elif name in duptest:
                if "Sample" not in errors:
                    errors["Sample"] = []
                errors["Sample"].append("Error: Duplicate sample name detected: " + name)
            else:
                for c in name:
                    if not c.isalnum() and c != '-':
                        if "Sample" not in errors:
                            errors["Sample"] = []
                        errors["Sample"].append("Error: sample name: " + name + " is malformed.(Only alphanumeric and '-' are allowed)")
                        break
                    
            duptest.add(name)
    
    if experiment == "":
        errors["Experiment"]= ["Failed to parse the name of the experiment. Are you sure you gave the script a sample sheet (and not some other file)?"]
    else:
        out["Experiment"] = ["'" + experiment + "' is the name of the experiment directory"]
 
    return out, errors


def parseSamplesheetByPath(path):
    '''
    Parse a sample sheet at the given path.
    @param path: The path to the sheet.
    @return: Tuple with two dictionaries: general information and actual errors.
    '''
    fh = open(path, "r")
    toR = parseSamplesheet(fh)
    fh.close()
    return toR


class ValidateSampleSheet(whodunargs.StandardProgram):
    '''Definition of the command line arguments.'''
    def __init__(self):
        '''Set up the program.'''
        whodunargs.StandardProgram.__init__(self)
        self.name = "ValidateSampleSheet"
        self.summary = "Run checks on the sample sheet."
        self.usage = "python3 validateSS.py --src path.csv --gui"
        self.version = "validateSS 0.0\nCopyright (C) 2024 August Woerner\nLicense LGPLv3: GNU LGPL version 3\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n"
        # --src
        self.srcOpt = whodunargs.ArgumentOptionFileRead("--src","Sample Sheet","The sample sheet to verify.","--src sheet.csv")
        self.srcOpt.validExts.add(".csv")
        self.options.append(self.srcOpt)
        # --log
        self.logOpt = whodunargs.ArgumentOptionFileWrite("--log","Output File","The file to write any information to.","--log results.txt")
        self.options.append(self.logOpt)
        # --gui
        self.guiOpt = whodunargs.ArgumentOptionFlag("--no-gui","No GUI","Do not display a gui.")
        #self.options.append(self.guiOpt)
    def baseRun(self):
        # figure out the sample sheet to check
        ssStream = sys.stdin
        if len(self.srcOpt.value) > 0:
            ssStream = open(self.srcOpt.value, "r")
        # get any problems
        hotOE = parseSamplesheet(ssStream)
        ssStream.close()
        # if there is a log file, write to it
        logStream = sys.stdout
        if len(self.logOpt.value) > 0:
            logStream = open(self.logOpt.value, "w")
        for e, vals in hotOE[1].items():
            for val in vals:
                logStream.write("ERROR " + e + " : " + val + "\n")
        for e, vals in hotOE[0].items():
            for val in vals:
                logStream.write("INFO " + e + " : " + val + "\n")
        logStream.close()
        # if any error entries, throw
        if len(hotOE[1]) > 0:
            raise ValueError("Sample Sheet has problems, check log.")


if __name__ == "__main__":
    mainPro = ValidateSampleSheet()
    mainPro.parse(sys.argv[1:])
    mainPro.run()

