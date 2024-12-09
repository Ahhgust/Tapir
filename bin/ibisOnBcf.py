#!/usr/bin/env python3

# Runs ibis on a BCF file...
# assumes GRCh38

import sys
import os

#gmap = "/eva/staging/augustw/public/genome_simulations/genetic_maps/ibis.sexaveraged.GRCh38.gmap"
gmap="/eva/edatums/reference_materials/genome_simulations/genetic_maps/ibis.sexaveraged.GRCh38.gmap"

args = sys.argv[1:]

if len(args) < 1:
    print("Gimme a BCF file!", file=sys.stderr)




for bcf in args:
    
    if not os.path.exists(bcf + ".csi"):
        print(bcf , " is not indexed. Please fix this!", file=sys.stderr)
        continue

    bn = os.path.splitext(os.path.basename(bcf))[0]

    if not os.path.exists(bn):
        os.mkdir(bn)

    # convert to plink format.
    if not os.path.exists(os.path.join(bn, bn + ".bim")):
        plinkCommand = "plink --memory 20000 --bcf " + bcf + " --snps-only --set-missing-var-ids @:#_\$1_\$2 --allow-extra-chr --biallelic-only --make-bed --double-id --out " + os.path.join(bn, bn)
        os.system(plinkCommand)

    # update the bim file w/ the genetic map
    if not os.path.exists(os.path.join(bn, bn + ".bim.cm")):
        gmapCommand = "add-map-plink.pl -noheader -cm " + os.path.join(bn, bn + ".bim") + " " + gmap + " > " + os.path.join(bn, bn + '.bim.cm')
        os.system(gmapCommand)


    # extra arguments are passed thru
    ibisCommand = "ibis " + os.path.join(bn, bn + ".bed") + " "  + os.path.join(bn, bn + ".bim.cm") + ' '  + os.path.join(bn, bn + ".fam") + " -2  -printCoef -gzip -noConvert -o " + os.path.join(bn, "ibis")
    #print(ibisCommand + " " + " ".join(args[1:]) + " &> /dev/null ")
    os.system(ibisCommand + " " + " ".join(args[1:]) + " &> /dev/null ")

    # okay. that's ugly.
    break

