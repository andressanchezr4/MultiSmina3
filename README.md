# MultiSmina3
MultiSmina is a new wrapper script for Python 3, updated from the version found [here](https://sourceforge.net/projects/smina/files/).

- It works by splitting the input SDF file into an arbitrary number of parts, each of which is processed simultaneously by a separate core.
  
- In order to provide feedback about the task progress, it can use the progressbar module if installed.

A flag ```-s```/```--smina_path``` was explicitly added to specify the path of smina executable.

Example of use:
```
python3 multi_smina3.py --smina_path /path/to/smina.static --ligand /path/to/molecules.sdf --receptor /path/to/receptor.pdbqt --autobox_ligand /path/to/ligand.pdb --out /path/to/output.sdf --scoring vinardo 
```
Other smina flags may also be included.

### Disclaimer:
If after running multismina the output is an empty file, make sure the ```--ligand``` sdf file (that should contain all the ligands to be docked) is properly built. 
