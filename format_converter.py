#!/usr/bin/python
"""A script for converting between the binary and pdb formats."""

import argparse, sys
import bin_to_pdb_conversion
import pdb_to_bin_conversion

def pars_inp():
    """Return the results of parsing the input."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", required=True,
                        help="The input file.")
    parser.add_argument("-o", "--outputFile", required=False, default="",
                        help="The output file.")
    parser.add_argument("--inputLaminFile", required=False, default="", help="The input lamin file.")
    parser.add_argument("--outputLaminFile", required=False, default="", help="The output lamin file.")
    parser.add_argument('--bin_to_pdb', action='store_true', help="Convert binary format to pdb.")
    parser.add_argument('--pdb_to_bin', action='store_true', help="Convert pdb format to binary.")
    parser.add_argument("--fpkf", type=int, default=10, required=False,
                        help="""Number of frames per keyframe, equal to the number of delta frames + 1.
                                Used only when converting from pdb to binary (defaul 10).""")
    parser.add_argument('-s', '--step', type=int, default=0,
                        help="""The step between consecutive binary frames converted to the pdb format.
                                The (default) value zero leads to converting all the frames in which the energy changed.
                                Used only when converting from binary to pdb.""")
    args = parser.parse_args()

    if (args.bin_to_pdb and args.pdb_to_bin):
        print "Error, you cannot provide both the --bin_to_pdb and --pdb_to_bin options"
        sys.exit(1)

    spl = args.inputFile.split(".")
    ext = spl[-1] # Extension of the input file 
    file_pref = "".join(spl[:-1]) # Part of the input file without extension 

    if args.bin_to_pdb:
        bin_to_pdb = True
    elif args.pdb_to_bin:
        bin_to_pdb = False
    else: # guess the conversion direction from the file extension
        if ext == 'bin':
            bin_to_pdb = True
        elif ext == 'pdb':
            bin_to_pdb = False
        else:
            print "Error, you must provide the --bin_to_pdb or --pdb_to_bin option or the input file should have bin or pdb extension"
            sys.exit(1)

    if args.outputFile:
        out_name = args.outputFile
    elif bin_to_pdb and ext == "bin":
        out_name = file_pref + ".pdb"
    elif not bin_to_pdb and ext == "pdb":
        out_name = file_pref + ".bin"
    else:
        print """Error, you must provide the output file name, or
        the input file should have a bin extension in the case of bin to pdb conversion,
        or the input file should have a pdb extendion in the case of pdb to bin conversion.
        In the latter two cases the extension for output will be changed to pdb and bin 
        respectively."""
        sys.exit(1)
    if args.inputLaminFile:
        if args.outputLaminFile:
            out_lamin = args.outputLaminFile
        elif bin_to_pdb:
            out_lamin = args.inputLaminFile.split('.bin')[0] + '.pdb' 
        else:
            out_lamin = args.inputLaminFile.split('.pdb')[0] + '.bin'
    else: # do not convert the lamin file
        out_lamin = "" 
    
    return args.inputFile, out_name, bin_to_pdb, args.fpkf, args.step, args.inputLaminFile, out_lamin

input_name, out_name, bin_to_pdb, fpkf, step, in_lamin, out_lamin = pars_inp()
if bin_to_pdb:
    assert step >= 0
    bin_to_pdb_conversion.convert_bin_to_pdb(input_name, out_name, step, in_lamin, out_lamin)
else:
    pdb_to_bin_conversion.convert_pdb_to_bin(input_name, out_name, fpkf, in_lamin, out_lamin)
