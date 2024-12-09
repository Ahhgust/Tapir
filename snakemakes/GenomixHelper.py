import os, sys

class GenomixHelper:

    def __init__(self, current_basedir, config):
        """
        Takes in the result of: workflow.current_basedir
        and config (dictionary! not filename)
        """
        self.ROOT=os.path.abspath( os.path.join(current_basedir, "..") )
        self.config=config


    def getResources(self, arg):
        """
        Requests a genomic resource. Pulled from the relevant config
        and returns as an absolute path to said resource
        Note that the resource definition in the config file can use relative paths
        """
        config=self.config
        
        if arg in config["refs"]:
            ref = config["refs"][arg]
        else:
            print("Failed to get the resource: ", arg, file=sys.stderr)
            exit(1)

        if ref.startswith("/"):
            return ref
            
        return os.path.join(self.ROOT, ref)
        
    def getParam(self, arg, key):
        """
        takes in the name of a genomic tool (eg, gatk)
        and returns a callable executable (by default; when key==binary)
        note that tools ending in .py are prepended with python3
        and that tools ending in .R are prepended with Rscript
        """
        
        config=self.config
        # the config file is structured as:
        # gatkParams:
        #    binary: 
        if not arg.endswith("Params"):
            arg += "Params"
        
        # looks for gatk
        if arg not in config:
            print("Failed to get the resource: ", arg, file=sys.stderr)
            exit(1)
        
        # looks for the binary associated w/ gatk
        if key not in config[arg]:
            print("Failed to get the resource binary: ", arg, key, file=sys.stderr)
            exit(1)
            
        prog = config[arg][key]

        exe=""
        # if there's no whitespace (eg, Rscript foo.R)
        # then 
        i=prog.find(" ")
        if i < 0:
            if prog.endswith(".R"):
                exe = "Rscript "
            elif prog.endswith(".pl"):
                exe = "perl "
            elif prog.endswith(".py"):
                exe="python3 " # sorry yo. no python2 support.
            else: # we have to assume that the user knows what they're doing
                1
        elif i==0:
            print("Leading whitespace found with: ", prog, "Please tix this", sep="\n", file=sys.stderr)
            exit(1)


        apath=""
        # paths may be relative to the parent directory of the snakemake file 
        # in which case, they must start with bin/
        if not prog.startswith("/") and prog.find("/") >= 0:
            apath=self.ROOT


        if exe == "":
            return os.path.join(apath,prog)
        return exe + " " + os.path.join(apath,prog)
        
    
    def getBinary(self, arg):
        return self.getParam(arg, key="binary")
    
    def getSummarizer(self, arg):
        """
        Some tools require a post-hoc script to summarize the output. 
        It's embedded as the "summarizer" key on the config (eg, config[gatkParams][summarizer])
        """
        return self.getParam(arg, key='summarizer')
