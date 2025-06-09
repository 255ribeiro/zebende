from pathlib import Path
import shutil

# Determine the directory of the current script
script_dir = Path(__file__).parents[1]

print(script_dir)

# Define the relative path to the source folder
source_folder = script_dir / 'zig-out/'

# Define the relative path to the destination folder
destination_folder = script_dir / 'zebende/zig_libs'

ignore = shutil.ignore_patterns('lib', '*.pdb')
# Copy the folder
shutil.copytree(source_folder, destination_folder, dirs_exist_ok=True, ignore=ignore)

print("Folder copied successfully from \n{} \nto \n{}.".format(source_folder, destination_folder))