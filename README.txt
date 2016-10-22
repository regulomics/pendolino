This is a project for running simulations of chromatin chains and for converting between various simulation outputs. 

1. Usage

The project contains two executable python scripts: pendolino.py and format_converter.py.
For details about their arguments see their help, by running "./pendolino.py -h" or "./format_converter.py -h" respectively.
Here we shall describe only some of their simplest functionalities.

1.1 pendolino.py

The script pendolino.py is used for running the simulation. Example usage:

python pendolino.py --regularBSitesFile dros_genom_regular.txt --laminBSitesFile dros_genom_lamin.txt -p all_genom2 -s 1000 -l 1151,1058,1228,1396,68,1122 -r 1151,1228 -b 10000 --fpkf 10 --binary

It will output the results of a 1000-step simulation into two binary files: 
all_genom2.bin - containing information about the positions of chromatin chains and binders after each fpkf - 1 moves. 
all_genom2_lamin.bin - containing the information about positions of lamins.

The formats of the binary files are specified by the message.proto protocol buffers file and in the Bachelor thesis 
"Efficient chromatin movement simulation" (see its Appendix B).
Some of the additional parameters of pendolino.py include --simulation_name SIMULATION_NAME 
and --simulation_description SIMULATION_DESCRIPTION, which can be used for setting the simulation name and description 
respectively. 
Without the option --binary analogous files but in the pdb format will be outputted. 

1.2 format_converter.py

The script format_converter.py can be used for converting between the binary and pdb files.
Example of conversion from binary to pdb:

python format_converter.py -i all_genom2.bin --inputLaminFile all_genom2_lamin.bin -s 10

It will convert all_genom2.bin to all_genom2.pdb and all_genom2_lamin.bin to all_genom2_lamin.pdb.
The optional parameter -s denotes the step between consecutive binary frames converted to the pdb format. 
The (default) value zero leads to converting all the frames in which the energy changed.

Example of conversion from pdb to binary format:

python format_converter.py -i all_genom2.pdb --inputLaminFile all_genom2_lamin.pdb --fpkf 10

It will convert all_genom2.pdb to all_genom2.bin and all_genom2_lamin.pdb to all_genom2_lamin.bin. 
The optional parameter -fpkf (default 10) denotes the number of frames per keyframe in the binary file all_genom2.bin,
so that the number of delta frames per keyframe is fpkf -1. 
The program tries to find a minimum number of movements transforming one pdb frame into another
and produces a binary frame with such moves in delta-frames after each fpkf - 1 moves found (if the number of 
moves is not divisible by fpkf - 1 then the remaining delta-frames are filled with zeros).

The parameter --inputLaminFile for both orders of conversion is optional. If it is not provided then only the 
simulation files are converted. 
The conversion direction in the script is inferred from the extension of 
the parameter after -i, but it can be also specified using the --bin_to_pdb or --pdb_to_bin flags. 
The name of the output simulation and lamin files can also be specified using the -o and --outputLaminFile flags.

2. File structure

Below we provide brief descriptions of all the files in the repository (apart from this README.txt file).

pendolino.py - contains the implementation of simulation of chromatin movement, as well as initialization of chromatin positions.
writing_binary.py - contains functions for outputting to the binary format. 
writing_pdb.py - contains functions for outputting to the pdb format.
constants.py - contains some constants used throughout th project.
format_converter.py - a script for converting between the binary and pdb formats, calls functions from the below two modules.
bin_to_pdb_conversion.py - contains functions for converting from the binary to the pdb format.
pdb_to_bin_conversion.py - contains functions for converting from the pdb to the binary format.
message.proto - the protostream format specification. 
message_pb2.py - protostream library generated automatically from message.proto.
dros_genom_regular.txt - contains sample positions of regular binding sites.
dros_genom_lamin.txt - contains sample positions of binding sites of both lamins and regular binders.


