#!/usr/bin/env python3

"""
Given the original (full) rigid receptor and the output of --out_flex,
which contains side chain coordinates from flexible docking with smina/gnina,
output the full restored receptors with the flexible residues re-inserted.
   
This is a bit fragile and may break if atoms do not have standard names.
"""

import prody, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Assemble full receptor from flexible docking results."
)
parser.add_argument("rigid", type=str, help="Rigid receptor (pdb)")
parser.add_argument("flex", type=str, help="Flexible sidechains from docking (pdb)")
parser.add_argument("out", type=str, help="Output file name (pdb)")
args = parser.parse_args()

rigidname = args.rigid
flexname = args.flex
outfile = args.out

out = open(outfile, "w")

flex = prody.parsePDB(flexname)
flexres = set(zip(flex.getChids(), flex.getResnums()))
backbone = {
    "N",
    "O",
    "H",
    "HN",
}  # C and CA are included in the flex part, but don't move
PDBLINE = "%s%-4s%s%8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"

# Print flexible residues
print("Flexres:", flexres)

def atype_perception(atype, aname):
    """
    Automatically perceive atom type from atom name, when the atom type is not
    included in the PDB file.

    Assumpsions:
        * Atom names start either with a number or the correct element name
        * Atom names follow the periodic table (i.e. HG is a type of hydrogen,
          Hg is mercury)
    """

    atype = atype.strip()
    
    if atype == "": # PDB file is missing atom types
        # Remove numbers from atom name
        atype = "".join(c for c in aname if not c.isdigit())
        atype = atype[:2]

        # Check if two-letter name is element (upper + lowercase) or custom
        # atom name (e.g. HG is hydrogen of carbon CG, Hg is mercury)
        if atype[-1].isupper():
            atype = atype[0]

    return atype

for ci in range(flex.numCoordsets()):
    which = defaultdict(int)
    out.write("MODEL %d\n" % ci)
    for line in open(rigidname):
        if line.startswith("ATOM"):
            chain = line[21]
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()
            aname = line[12:16].strip()
            atype = atype_perception(line[76:].strip(), aname)

            if (chain, resnum) in flexres and aname not in backbone:

                if atype != "H":
                    resatoms = flex[chain].select("resnum %d and not name H" % resnum)
                    w = which[(chain, resnum)]
                    which[(chain, resnum)] += 1  # update to next index
                    atom = resatoms[w]  # this is the atom to replace this line with
                    c = atom.getCoordsets(ci)
                    line = PDBLINE % (
                        line[:13],
                        aname,
                        line[17:30],
                        c[0],
                        c[1],
                        c[2],
                        1.0,
                        0.0,
                        atype,
                    )
                else:
                    line = ""
        elif line.startswith("END"):
            line = ""

        out.write(line)
    
    out.write("ENDMDL\n")
out.write("END\n")
