"""A module containing constants used by other modules."""

chain_list = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T"]

EMPTY = 0
BINDER = 1
LAMIN = 2
BSITE_R = 3 # interacts only with binders 
BSITE_L = 4 # interacts both with lamins and binders
REGDNA = 5 # does not interact
DIST = 3
    
binding_forces = { # vectors of binding forces with lamins and binders of atoms in chains in different states
    BSITE_R: (0, 1), 
    BSITE_L: (1, 1), 
    REGDNA: (0, 0)
}

pdb_reg_names = { # names of regions of DNA used in PDB files
    REGDNA: "UNB",
    BSITE_R: "BOU",
    BSITE_L: "LAM"
}                 
