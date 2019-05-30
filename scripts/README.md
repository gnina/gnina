# Scripts

## MakeFlex

The script `makeflex.py` can be used after flexible docking to build a full receptor. Given the original receptor and the flexible side chains obtained via flexible docking (with `smina` or `gnina`), the script builds the full receptor with the flexible residues re-inserted.

Usage:
```
python makeflex.py RIGID.pdb FLEXIBLE.pdb OUT.pdb
```
where `RIGID.pdb` contains the original receptor, `RIGID.pdb` contains the flexible side chains (output of flexible docking) and `OUT.pdb` is the output containing the full receptor with the flexible resiudes re-inserted.

`makeflex.py` supports multiple models (`MODEL`/`ENDMDL`) in `FLEXIBLE.pdb`.