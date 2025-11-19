# Given a multifasta, this script will create a folder named job1...jobN with subfolders named data_inputs and inference_inputs.
# This script will create a custom json file for each fasta sequence, and place them in batches of batch size (e.g. 10) under each job1 to jobN folders
# the batch size can be provided by the user through a command line argument.
# usage: set_up_folders.py [-h] fasta_file batch_size

# for example, let test.fasta look like this, and batch_size = 5:
# >sequence1
# ABCDEFG
# >sequence2
# HIJKLMNOP
# >sequence3
# QRSTUV
# >sequence4
# ABCDEFG
# >sequence5
# HIJKLMNOP
# >sequence6
# QRSTUV
# This script will create these 2 job folders, with 5 json files in each.
# job1
# ├── data_inputs
# │   └── fold_input_1.json
# │   └── fold_input_2.json
# │   └── fold_input_3.json
# │   └── fold_input_4.json
# │   └── fold_input_5.json
# └── inference_inputs
# job2
# ├── data_inputs
# │   └── fold_input_6.json
# └── inference_inputs

## USAGE ##
# chmod +x set_up_directory.py
# python set_up_directory test.fasta batch_size

# Import libraries
import json
from Bio import SeqIO
import os
import argparse

def fasta_to_json(fasta_file, batch_size):
    # read in the fasta sequence from the user
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    # count how many sequences there are in the multi-fasta file
    count = len(sequences)

    
    if count >= batch_size:
        # Determine how many batchs of size batch_size there should be.
        # In the example able (6 + 5 - 1) // 5 = 2 --> 2 batches --> expected folder: job1 , job2
        num_batches = (count + batch_size - 1) // batch_size
        for batch_index in range(num_batches):
            # Make a jobN and data_inputs and inference_inputs for each folder
            job_dir = f"job{batch_index + 1}"
            data_inputs_dir = os.path.join(job_dir, "data_inputs")
            inference_inputs_dir = os.path.join(job_dir, "inference_inputs")
            
            os.makedirs(data_inputs_dir, exist_ok=True)
            os.makedirs(inference_inputs_dir, exist_ok=True)

            # Create a json file for that batch of sequences
            for i, record in enumerate(sequences[batch_index * batch_size : (batch_index + 1) * batch_size], start=1):
                data = {
                    "name": record.id,
                    "sequences": [
                        {
                            "protein": {
                                "id": ["A"],
                                "sequence": str(record.seq)
                            }
                        }
                    ],
                    "modelSeeds": [1],
                    "dialect": "alphafold3",
                    "version": 1
                }
                
                json_filename = os.path.join(data_inputs_dir, f"fold_input_{(batch_index * batch_size) + i}.json")
                with open(json_filename, "w") as json_file:
                    json.dump(data, json_file, indent=2)
    else:
        # If there are fewer sequences than the batch size, then there should only be one folder created, i.e. job1.
        job_dir = "job1"
        data_inputs_dir = os.path.join(job_dir, "data_inputs")
        inference_inputs_dir = os.path.join(job_dir, "inference_inputs")
        
        os.makedirs(data_inputs_dir, exist_ok=True)
        os.makedirs(inference_inputs_dir, exist_ok=True)

        # create the json files 
        for i, record in enumerate(sequences, start=1):
            data = {
                "name": record.id,
                "sequences": [
                    {
                        "protein": {
                            "id": ["A"],
                            "sequence": str(record.seq)
                        }
                    }
                ],
                "modelSeeds": [1],
                "dialect": "alphafold3",
                "version": 1
            }
            
            json_filename = os.path.join(data_inputs_dir, f"fold_input_{i}.json")
            with open(json_filename, "w") as json_file:
                json.dump(data, json_file, indent=2)
    
    print("JSON files saved under respective job folders in the current directory")

# Help page to let the user know how to use this script.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA file to JSON format with structured job directories.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file")
    parser.add_argument("batch_size", type=int, help="Number of JSON files per job folder")
    
    args = parser.parse_args()
    fasta_to_json(args.fasta_file, args.batch_size)
