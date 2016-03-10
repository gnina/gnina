# BIN='~/git/gnina_clone/build/linux-gpu/debug/smina.gpu'
# ARGS='-r receptor.pdb -l MP5.pdb --autobox_ligand MES.pdb --gpu'

# $BIN $ARGS

../build/linux-gpu/release/smina.gpu -r receptor.pdb -l MP5.pdb --autobox_ligand MES.pdb --gpu --minimize

