
# -*- coding: utf-8 -*-

import ntpath 
import os




def get_filename_without_extension(file_path):
    file_basename = os.path.basename(file_path)
    filename_without_extension = file_basename.split('.')[0]
    return filename_without_extension
    
def files_to_PPIFDT (files_paths):
    files = []
    files_to_PPIFDT = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(files_paths):
        for file in f:
            if '.fasta' in file:
                files.append(os.path.join(r, file))

    for f in files:
       files_to_PPIFDT.append(get_filename_without_extension(f))
    return files_to_PPIFDT