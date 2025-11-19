#!/usr/bin/env python3
"""
Generate AlphaFold3 fold_input.json files for multiple Piwi protein sequences,
including a 3' 2'-O-methyl modification on the final nucleotide of the piRNA.

Each output JSON is placed under:
    job#_name/data_inputs/fold_input.json
with a companion empty inference_inputs/ directory created beside it.

Modification is applied using PDB CCD codes:
    A → MAD
    C → OMC
    G → OMG
    U → OMU
"""

import os
import csv
import json
import argparse

def get_3prime_methyl_modifications(seq):
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("piRNA sequence is empty.")

    last_base = seq[-1]
    pos = len(seq)

    ccd_map = {"A": "MAD", "C": "OMC", "G": "OMG", "U": "OMU"}

    if last_base not in ccd_map:
        raise ValueError(f"Unexpected base for methylation: {last_base}")

    return [{"modificationType": ccd_map[last_base], "basePosition": pos}]

def generate_fold_input(name, protein_seq, piRNA_seq, target_RNA_seq):
    sequences = [
        {
            "protein": {
                "id": "A",
                "sequence": protein_seq.strip().replace("\n", "")
            }
        },
        {
            "rna": {
                "id": "R",
                "sequence": piRNA_seq.strip().replace("\n", "").upper(),
                "modifications": get_3prime_methyl_modifications(piRNA_seq)
            }
        }
    ]

    if target_RNA_seq and target_RNA_seq.strip() != "":
        sequences.append({
            "rna": {
                "id": "T",
                "sequence": target_RNA_seq.strip().replace("\n", "").upper()
            }
        })

    return {
        "name": name,
        "sequences": sequences,
        "modelSeeds": [1, 2, 3],
        "dialect": "alphafold3",
        "version": 1
    }

def main():
    parser = argparse.ArgumentParser(description="Generate fold_input.json files for AlphaFold3 jobs.")
    parser.add_argument("--csv", required=True, help="Input CSV (columns: name,fasta_protein_sequence)")
    parser.add_argument("--piRNA", required=True, help="piRNA sequence")
    parser.add_argument("--targetRNA", default="", help="Optional target RNA sequence")
    parser.add_argument("--outdir", default="jobs", help="Output directory for job folders")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        for i, row in enumerate(reader, start=1):
            name = row["name"].strip()
            protein_seq = row["fasta_protein_sequence"].strip()

            job_dir = os.path.join(args.outdir, f"job{i}_{name}")
            data_inputs_dir = os.path.join(job_dir, "data_inputs")
            inference_inputs_dir = os.path.join(job_dir, "inference_inputs")

            os.makedirs(data_inputs_dir, exist_ok=True)
            os.makedirs(inference_inputs_dir, exist_ok=True)

            data = generate_fold_input(name, protein_seq, args.piRNA, args.targetRNA)

            json_path = os.path.join(data_inputs_dir, "fold_input.json")
            with open(json_path, "w") as f:
                json.dump(data, f, indent=2)

            print(f"✅ Created {json_path}")

    print(f"\nAll jobs created under: {os.path.abspath(args.outdir)}")

if __name__ == "__main__":
    main()