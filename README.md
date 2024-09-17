# MultiSmina3
MultiSmina is a new wrapper script for Python 3, updated from the version found [here](https://sourceforge.net/projects/smina/files/).

- It works by splitting the input SDF file into an arbitrary number of parts, each of which is processed simultaneously by a separate core.
  
- In order to provide feedback about the task progress, it can use the progressbar module if installed.

A flag -s/--smina_path was added to include the path of smina.static.

