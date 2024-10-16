import os

# Define the directory and output file
directory = "output/"
output_file = "Design_list.txt"

# Get all file names in the specified directory
file_names = os.listdir(directory)

# Write the file names to the output file, one name per line
with open(output_file, 'w') as f:
    for file_name in file_names:
        f.write("output/"+file_name+"\n")

