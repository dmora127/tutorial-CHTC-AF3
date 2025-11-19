# Predicting Protein Structures with AlphaFold3 on the CHTC GPU Capacity

A two-phase workflow: alignment generation → structure prediction

AlphaFold3 (AF3) is a next-generation biomolecular structure prediction system capable of modeling multichain protein complexes, DNA/RNA interactions, ligands, and modified residues (e.g., 2′-O-methylated piRNAs). This guide provides a step-by-step workflow for running AF3 on the CHTC HTC ecosystem, including how to structure your project, prepare input JSON files, run many inference jobs with HTCondor, and store your results using OSDF.

AlphaFold3 workloads in high-throughput environments are best organized into two separate job types:
* **Step 1: Generating the Alignments (Data-Only Pipeline)**
  * MMseqs2 search
  * HMMER search
  * Template retrieval (if enabled)
  * MSA + templates packaged into an AF3 “fold_input.json” 
* **Step 2: Predicting the Structure (Inference Pipeline)**
  * Load the precomputed fold_input.json 
  * Run the AF3 diffusion model 
  * Produce PDB models, ranking, trajectories, and metrics


You will learn how to:

* Basecall raw Nanopore reads using the latest GPU-accelerated Dorado basecaller  
* Use the OSPool's GPU capacity to accelerate basecalling with Dorado  
* Break down massive bioinformatics workflows into many independent smaller tasks  
* Submit hundreds to thousands of jobs with a few simple commands  
* Use the Open Science Data Federation (OSDF) to manage file transfers during job submission  

All of these steps run across hundreds (or thousands) of jobs using the HTCondor workload manager and Apptainer containers to execute your software reliably and reproducibly at scale. The tutorial uses realistic genomics data and emphasizes performance, reproducibility, and portability. You will work with real data and see how high-throughput computing (HTC) can accelerate your workflows.

![Workflow_Diagram.png](.images/Workflow_Diagram.png)

> [!NOTE]
> If you are new to running jobs on the OSPool, complete the HTCondor ["Hello World"](https://portal.osg-htc.org/documentation/htc_workloads/submitting_workloads/tutorial-quickstart/) tutorial and our ["Submit Jobs to the OSPool"](https://portal.osg-htc.org/documentation/htc_workloads/workload_planning/htcondor_job_submission/) guide before starting this tutorial.

**Let’s get started!**

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [tutorial-ONT-Basecalling](#tutorial-ont-basecalling)
- [Long-Read Genomics on the OSPool](#long-read-genomics-on-the-ospool)
   * [Tutorial Setup](#tutorial-setup)
      + [Assumptions](#assumptions)
      + [Materials](#materials)
   * [Understanding Basecalling in Oxford Nanopore Sequencing](#understanding-basecalling-in-oxford-nanopore-sequencing)
      + [From Signal to Sequence: The Role of Basecalling](#from-signal-to-sequence-the-role-of-basecalling)
   * [Basecalling on the OSPool by Sequencing Channel](#basecalling-on-the-ospool-by-sequencing-channel)
   * [Recommended Directory Structure](#recommended-directory-structure)
   * [Basecalling Oxford Nanopore Long Reads Using Dorado](#basecalling-oxford-nanopore-long-reads-using-dorado)
      + [Set Up Your Software Environment](#set-up-your-software-environment)
      + [Data Wrangling and Splitting Reads](#data-wrangling-and-splitting-reads)
         - [Downloading the Dorado Basecalling Models](#downloading-the-dorado-basecalling-models)
         - [Split Your Reads for Basecalling](#split-your-reads-for-basecalling)
      + [Submit Your Basecalling Jobs](#submit-your-basecalling-jobs)
      + [Post-Basecalling Steps](#post-basecalling-steps)
   * [Next Steps](#next-steps)
      + [Software](#software)
      + [Data](#data)
      + [GPUs](#gpus)
   * [Getting Help](#getting-help)

<!-- TOC end -->


## Tutorial Setup

### Assumptions

This tutorial assumes that you:

* Have basic command-line experience (e.g., navigating directories, using bash, editing text files)
* Have a working OSPool account and can log into an Access Point (e.g., ap40.uw.osg-htc.org)
* Are familiar with HTCondor job submission, including writing simple `.sub` files and tracking job status with `condor_q`
* Understand the general workflow of long-read sequencing analysis: basecalling → mapping → variant calling
* Have access to a machine with a GPU-enabled execution environment (provided automatically via the OSPool)
* Have sufficient disk quota and file permissions in your OSPool home and OSDF directories

> [!TIP]
> You do not need to be a genomics expert to follow this tutorial. The commands and scripts are beginner-friendly and self-contained while reflecting real-world research workflows.

### Prerequisites

1. [ ] A CHTC HTC account. If you do not have one, request access at the [CHTC Account Request Page](https://chtc.cs.wisc.edu/uw-research-computing/form.html).
2. [ ] Basic familiarity with HTCondor job submission. If you are new to HTCondor, complete the [HTCondor "Hello World" Tutorial](https://portal.osg-htc.org/documentation/htc_workloads/submitting_workloads/tutorial-quickstart/) and read the [Submit Jobs to the OSPool Guide](https://portal.osg-htc.org/documentation/htc_workloads/workload_planning/htcondor_job_submission/).
3. [ ] AlphaFold3 Model Weights. Request the AF3 model weights from the [DeepMind AlphaFold Team](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-model-parameters).

> [!WARNING]
> Requesting AlphaFold3 model weights requires agreeing to DeepMind's terms of service. Ensure you comply with all licensing and usage restrictions when using AF3 for research. This tutorial does not distribute AF3 model weights. Requesting the weights can take up to several weeks. Ensure you have them before starting the tutorial.

Log into your OSPool account:

    ```bash
    ssh user.name@ap##.uw.osg-htc.org
    ```

To obtain a copy of the tutorial files, you can:

* Clone the repository:

  ```bash
  git clone https://github.com/dmora127/tutorial-CHTC-AF3.git
  cd tutorial-CHTC-AF3/
  bash tutorial-setup.sh <username>
  ```
  _This script creates the directory structure in your home directory `/home/<user.name>/tutorial-ONT-Basecalling/` and OSPool directory `/ospool/ap40/<user.name>/tutorial-ONT-Basecalling/`, along with several subdirectories used in this tutorial._

* Or download the toy dataset using Pelican **CURRENT UNAVAILBLE - WORK IN PROGRESS**:

  ```bash
  pelican object get pelican://osg-htc.org/ospool/uw-shared/OSG-Staff/osg-training/tutorial-ospool-genomics/data/path/to/pod5/files ./
  ```


## Understanding the AlphaFold3 Workflow

[INSERT NARRATIVE ON ALPHAFOLD3 WORKFLOW HERE]

### The CPU-Only Pipeline: Generating Alignments (Step 1)

[EXPLAIN THE ALIGNMENT GENERATION STEP HERE]

## The GPU-Accelerated Pipeline: Structural Prediction (Step 2)

[EXPLAIN THE STRUCTURE PREDICTION STEP HERE]

## Recommended Directory Structure

Use a layout like the one below (directories may already exist if you ran the setup script):

```bash
./tutorial-CHTC-AF3/
├── scripts                            # scripts for preprocessing and running dorado jobs
│   ├── preprocessing_pod5s.sh
│   └── run_dorado.sh
├── list_of_pod5_files.txt             # one POD5 path per line
├── logs                               # condor .log files
├── README.md                          # Tutorial step-by-step documentation and instructions
├── data_pipeline.sub                  # HTCondor submit file (CPU-Intensive Alignment Generation)
├── inference_pipeline.sub             # HTCondor submit file (GPU-Intensive Structure Prediction)
├── software                           # containers recipes 
│   └── alphafold3.def
└── tutorial-setup.sh                  # optional helper to set up the tree
```

You also have a companion **OSDF** directory for storing large files, such as containers and Dorado models:

```bash
/staging/<netID>/tutorial-CHTC-AF3/
├── data/
│   ├── dna_r10.4.1_e8.2_400bps_fast@v4.2.0_5mCG_5hmCG@v2.tar.gz
│   ├── dna_r10.4.1_e8.2_400bps_fast@v4.2.0.tar.gz
│   ├── rna004_130bps_sup@v5.2.0.tar.gz
├── software/
│   └── dorado_build1.2.0_27OCT2025_v1.sif
```

Run the included `tutorial-setup.sh` script in the companion repository to create this structure.

## Basecalling Oxford Nanopore Long Reads Using Dorado

### Set Up Your Software Environment
Before basecalling, set up your software environment to run Dorado inside an Apptainer container.

1. Set your temporary Apptainer build directory to `/home/tmp/`. Once you've built your container, you can delete the contents of this directory to reduce quota usage on `/home`. 

   On the Access Point (AP), run:
    ```bash
    mkdir -p $HOME/tmp
    export TMPDIR=$HOME/tmp
    export APPTAINER_TMPDIR=$HOME/tmp
    export APPTAINER_CACHEDIR=$HOME/tmp
    ```
> [!CAUTION]
> Run these commands **every time you log in or build a new container**. Building Apptainer containers without setting these variables places excessive strain on shared storage resources and **violates OSPool usage policies**. Failure to follow these steps may result in restricted access.

2. Change to your `software/` directory in your tutorial folder:
    ```bash
    cd ~/tutorial-ONT-Basecalling/software/
    ```

3. Create a definition file for Apptainer to build your Dorado container. Open a text editor, such as `vim` or `nano`, and save the following as `software/dorado.def`:

    ```bash
    Bootstrap: docker
    From: nvidia/cuda:13.0.1-cudnn-runtime-ubuntu22.04
    
    %post
        DEBIAN_FRONTEND=noninteractive
    
        # system packages
        apt-get update -y
        apt-get install -y python3-minimal curl
        curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
        python3 get-pip.py --break-system-packages
        rm get-pip.py 
        apt-get install -y bedtools
    
        # install Dorado and POD5
        cd /opt/
        curl -L https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz -o ./dorado-1.2.0-linux-x64.tar.gz
        tar -zxvf dorado-1.2.0-linux-x64.tar.gz
        rm dorado-1.2.0-linux-x64.tar.gz
    
        # install POD5 using pip
        pip install pod5 --break-system-packages
    
    %environment
        # set up environment for when using the container
        # add Dorado to $PATH variable for ease of use
        export PATH="/opt/dorado-1.2.0-linux-x64/bin/:$PATH"
    ```

    This definition file uses the Nvidia CUDA 13.0.1 libraries on an Ubuntu 22.04 base image and installs necessary packages to run Dorado and POD5 in an Apptainer container.


4. Build your Apptainer container on the Access Point (AP):
    ```bash
    apptainer build dorado_build1.2.0_27OCT2025_v1.sif dorado.def
   ```
   > [!WARNING]
   > This will take a few minutes depending on system usage. If you encounter any errors during the build process, double-check that you have set your temporary directories correctly (see step 1). Contact your RCF team if issues persist.

5. Move your finalized container image, `dorado_build1.2.0_27OCT2025_v1.sif`, to your `OSDF` directory
    
    ```bash
   mv dorado_build1.2.0_27OCT2025_v1.sif /ospool/ap40/data/<user.name>/tutorial-ONT-Basecalling/software/
   ```
   
### Data Wrangling and Preparing AlphaFold3 Inputs

Oxford Nanopore sequencing runs generally yield POD5 files. Each POD5 file is generated about once an hour throughout the
duration of the sequencing run. This output format does not scale very well, as data output usually plateaus after 24-48hrs.
This would mean that POD5 files that are generated from earlier in the sequencing run, will be larger in file size compared to
files later in the run. Additionally, this division of data does not allow for _Duplex_ read basecalling. As a result prior
to running Dorado, we must first reorganize the data contained within all the POD5 files.

#### Setting Up AlphaFold3 Input JSONs and Job Directories
Alphafold3 requires input JSON files that contain the input query sequence(s) along with their corresponding metadata, such as chain IDs, sequence names, and molecule type. Additionally, if you are using precomputed MSAs and templates, these JSON files should also reference the paths to them. For this tutorial, we will create individual JSON files for each protein sequence. Each AlphaFold3 job will have its own directory containing the input JSON file and any associated supporting files, such as

```bash
./AF3_Jobs/Job3_ProteinZ/
├── data_inputs                            # data_inputs directory for the data pipeline
│   └── fold_input.json              # input JSON for data pipeline
└── inference_inputs                       # inference_inputs directory for the GPU inference pipeline
```

As your jobs progress through the two phases of the AlphaFold3 workflow, the output files will be organized as follows:

```bash
./AF3_Jobs/Job3_ProteinZ/
├── data_inputs                            # data_inputs directory for the data pipeline
│   └── fold_input.json              # input JSON for data pipeline
├── inference_inputs                       # inference_inputs directory for the GPU inference pipeline
│   └── node9.data_pipeline.tar.gz   # tar.gz output from data pipeline (input for inference pipeline)
├── data_pipeline_job.err                  # stderr from data pipeline runs
├── data_pipeline_job.out                  # stdout from data pipeline runs
├── inference_pipeline_job.err             # stderr from inference pipeline runs
├── inference_pipeline_job.out             # stdout from inference pipeline runs
└── node9.inference_pipeline.tar.gz        # tar.gz output from inference pipeline (final results)
```

We need to setup a job directory for each AlphaFold3 job we plan to run. Each job directory will contain the necessary input JSON files and will be used to store the output files generated during the two phases of the AlphaFold3 workflow. We can use a simple bash script to automate the creation of these job directories and input JSON files.

1. Setup a CSV manifest containing your protein sequences in FASTA format in the `data/protein_sequences/` directory. Each FASTA file should contain a single protein sequence. The CSV should follow this format: 

    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinA,protein,A,MKTAYIAKQRQIS
    ProteinB,protein,A,GAVLILALLAVF
    ```
    
    If you plan to model multiple molecules, simply add more columns to the CSV manifest. For example, if you have two proteins and one DNA molecule, your CSV might look like this:
    
    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinA,protein,A,MKTAYIAKQRQIS,protein,B,MKTAYIAKQRQIS
    ProteinB,protein,A,GAVLILALLAVF,dna,B,GCGTACGTAGCTAGC
    ```
    
    You can also model multimeric complexes by including multiple chain_ids in the CSV manifest. Each chain_id must be wrapped in double-quotes `(")` and the full chain_id set must be wrapped in single quotes `(')`. For example, if you have a trimeric protein complex, your CSV might look like this:
    
    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinComplex,protein,'"A", "B", "C"',MKTAYIAKQRQIS
    ```

2. Run the helper script in `scripts/setup_af3_jobs.sh` to read the CSV manifest and create the necessary job directories and input JSON files for each protein sequence.:

    ```bash
   ./scripts/setup_af3_jobs.sh --manifest data/protein_sequences/manifest.csv --output_dir AF3_Jobs/
    ```

3. Verify that the job directories and input JSON files have been created correctly:

    ```bash
   tree AF3_Jobs/
   ```

   You should see a list of job directories corresponding to each protein sequence in your CSV manifest.

#### Preparing Your _List of (AlphaFold) Jobs_

Often you will hear the term "list of jobs" when working with HTCondor. A "list of jobs" is a simple way to specify multiple jobs to be run using a single HTCondor submit file. If you are familar with HTCondor, you may have used the `queue` command in your submit files to specify multiple jobs. HTCondor uses this `queue` statement to process multiple job submissions in one single submit file. Each job is defined by a set of arguments passed to the executable. This constitutes a "list of jobs".

There are many ways to create a "list of jobs" directly in your HTCondor submit files. For example, you can use the `queue` command with a `from` clause to read job arguments from a file or use the `queue` command with a `matching files` clause to generate jobs based on files in a directory. However, when working with large numbers of jobs, it is often more convenient to use a separate file to list the jobs you want to run. This file is referred to as a "list of jobs".

 In this tutorial, we will create a "list of jobs" file that contains the names of each AlphaFold3 job directory we plan to run. This file will be used in our HTCondor submit files to specify which jobs to execute.

1. Generate the "list of jobs" file by listing the job directories in your `AF3_Jobs/` directory:

    ```bash
    ls AF3_Jobs/ > list_of_af3_jobs.txt
    ```

1. Examine the contents of your "list of jobs" file:

    ```bash
    cat list_of_af3_jobs.txt
    ```

> [!NOTE]  
> Each line in `list_of_af3_jobs.txt` corresponds to a single AlphaFold3 job that will be executed using HTCondor. You can modify this file to add or remove jobs as needed. This file, functionally, is a one-column comma-separated values (CSV) file where each line represents a job name. You could add additional columns to this file if you wanted to pass more variables to your HTCondor jobs.

6. Create a list of POD5 files to iterate through during job submission:

    ```bash
    ls split_pod5_subsets > ~/tutorial-ONT-Basecalling/list_of_pod5_files.txt
   ```
   
    If you `head` this new file you should see an output similar to this:

    ```
    [user.name@ap40 user.name]$ head ~/tutorial-ONT-Basecalling/list_of_pod5_files.txt
    channel-100.pod5
    channel-101.pod5
    channel-102.pod5
    channel-103.pod5
    channel-104.pod5
    channel-105.pod5
    channel-106.pod5
    channel-107.pod5
    channel-108.pod5
    channel-109.pod5
    [user.name@ap40 user.name]$ 
   ```

### Submit Your AlphaFold3 Jobs - CPU-Intensive Alignment Generation (Step 1)

[INSERT INSTRUCTIONS FOR SUBMITTING THE DATA PIPELINE JOBS HERE]

[INSERT GRAPHIC OF ALPHAFOLD3 WORKFLOW WITH DATA PIPELINE HIGHLIGHTED]

1. Change to your `tutorial-CHTC-AF3/` directory:
    ```bash
    cd ~/tutorial-CHTC-AF3/
   ```

2. Review your Data Pipeline executable script `scripts/data_pipleine.sh`. Generally, no changes will be necessary. You should, however, review the following information and options below as your AF3 jobs may require additional non-default options.

    #### Overview: AlphaFold3 Data Pipeline Executable (Data-Only Stage)

    This script implements the **data-generation portion** of the AlphaFold3 workflow — the part that creates alignments, MSA features, and template features required for structure prediction. It **does not** run structure inference step. It is designed for execution inside HTCondor jobs on CHTC systems and supports both containerized and pre-staged database workflows.
    
    The script handles:

   - Preparing a clean per-job working directory  
   - Detecting database availability (local `/alphafold3` vs. copying from staging)  
   - Copying and/or binding input JSON files  
   - Running the AlphaFold3 **data pipeline only**  
   - Tarballing the resulting `data_pipeline` outputs for transfer back to the Access Point  
   - Cleaning up the job sandbox to minimize disk use

    #### 1. Checking for Pre-Staged AF3 Databases

    CHTC hosts a copy of the AlphaFold3 databases locally on certain machines. Machines with these databases advertise this resource availability using the `HasAlphafold3` HTCondor MachineAd. The script is able to run on machines with/without these databases pre-loaded. The script inspects `.machine.ad` to see whether the matched machine advertises the availbility of AF3 databases:
    
    ```
    HasAlphafold3 = true
    ```
    
    The script also checks that the `/alphafold3/` path is not empty. If the MachineAd is set to `true` **and** `/alphafold3/` is populated:
    
    ✔ Uses the pre-extracted AF3 databases  
    ✘ Otherwise falls back to decompressing databases locally.

    The databases hosted by CHTC include the full AF3 default database suite:

    - PDB mmCIF  
    - MGnify Clusters 
    - UniRef90  
    - UniProt (Full UniProt Database) 
    - BFD First Non-Consensus Sequences
    - RNAcentral  
    - NT RNA clusters  
    - Rfam
   
    Databases may be updated periodically by CHTC staff. If you require a custom database set (e.g., reduced-size databases for testing), you can modify the script to use your own database tarballs. Contact the CHTC Research Computing Facilitation team for assistance.

    #### 2. Command-Line Options
    
    | Flag | Meaning |
    |------|---------|
    | `--work_dir_ext <name>` | Name appended to working directory (`work.<name>`) |
    | `--tmpdir <path>` | Override scratch location |
    | `--verbose` / `--silent` | Adjust verbosity |
    | `--no_copy` | Don’t copy containers/databases locally |
    | `--container <image>` | Path to Apptainer image |
    | `--random_sleep_minutes <N>` | Randomized delay to prevent thundering herd |
    | `--smalldb` | Use a reduced-size database set |
    | `--extracted_database_path <path>` | Use provided database directory |

    #### 3. Working Directory Setup
    
    Creates:
    
    ```
    work.<ext>/
        af_input/
        af_output/
        models/
        public_databases/
        tmp/
    ```
    
    Moves all `*.json` into `af_input/`.
    
    #### 4. Running the AlphaFold3 Data Pipeline
    
    Inside the container:
    
    ```
    python run_alphafold.py --run_data_pipeline=true --run_inference=false
    ```
    
    Generates all alignment and feature files under `af_output/<job_name>/`.
     
    #### 5. Packaging Results
    
    Each output directory is archived:
    
    ```
    <target>.data_pipeline.tar.gz
    ```
    
    These are returned to the submit host.

3. Create your submit file `data_pipeline.sub`. The submit file below works out-of-the-box (if you've setup your directories as specified in section [Setting Up AlphaFold3 Input JSONs and Job Directories]()). You can specify additional parameters for the exectuable in the `arguments` attribute as needed. Refer to the [AlphaFold3 Data Pipeline Executable - Command-Line Options]() section above for available options.

    ```bash
    # CHTC maintained container for AlphaFold3 as of January 2025
    #container_image = osdf:///ospool/uw-shared/OSG-Staff/public/alphafold3/alphafold3.minimal.22Jan2025_v1.sif
    container_image = file:///staging/groups/glbrc_alphafold/af3/alphafold3.minimal.22Jan2025.sif
    
    executable = scripts/data_pipeline.sh
    
    log = ../logs/data_pipeline.log
    output = data_pipeline_$(Cluster)_$(Process).out
    error  = data_pipeline_$(Cluster)_$(Process).err
    
    initialdir = $(directory)
    # transfer all files in the data_inputs directory
    transfer_input_files = data_inputs/
   
   # transfer output files back to the submit node
   transfer_output_files = $(directory).data_pipeline.tar.gz
   transfer_output_remaps = "$(directory).data_pipeline.tar.gz=inference_inputs/$(directory).data_pipeline.tar.gz"
    
    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    
    # We need this to transfer the databases to the execute node
    Requirements = (Target.HasCHTCStaging == true) && (TARGET.HasAlphafold3 == true)
    
    if defined USE_SMALL_DB
      # testing requirements
      request_memory = 8GB
      request_disk = 16GB
      request_cpus = 4
      arguments = --smalldb --work_dir_ext $(Cluster)_$(Process) --verbose
    else
      # full requirements
      request_memory = 8GB
      # Request less disk if matched machine already has AF3 DB preloaded (650GB savings)
      request_disk = 700000 - ( (TARGET.HasAlphaFold3?: 1) * 650000)
      request_cpus = 8
      arguments = --work_dir_ext $(Cluster)_$(Proc) 
    endif
    
    #queue directory matching job*
    queue directory from list_of_af3_jobs.txt
   ```

This submit file will read the contents of `list_of_af3_jobs.txt`, iterate through each line, and assign the value of each line to the variable `$(directory)`. This allows you to programmatically submit _N_ jobs, where _N_ equals the number of AlphaFold3 job directories you previously created. Each job processes one AlphaFold3 job directory and uses the CHTC-maintained AlphaFold3 container image, which is transferred to the Execution Point (EP) by HTCondor.

The submit files will attempt to match to machines that are advertising the `HasAlphafold3` resource, which ensures that the necessary AlphaFold3 databases are transferred to the EP for each job. You can remove the `&& (TARGET.HasAlphafold3 == true)` clause from the `Requirements` attribute if you want to run jobs on machines with and without pre-staged databases. Removing this clause will increase the number of available machines for your jobs, leading to faster _start_ times for your jobs, but may also increase job runtime due to long database transfer times.

The script will also check the matched machine's MachineAd, after the job has matched, to see if the `HasAlphafold3` attribute is set to `true`. If it is, the submit file will request significantly less disk space, as the databases are already present on the machine. If the attribute is not set to `true`, the submit file will request more disk space to accommodate the database transfer. Lower disk requests can lead to increased number of running jobs, as the scheduler has more flexibility in matching jobs to machines.

    > [!TIP]  
    > AlphaFold3 jobs can be resource-intensive, especially when using the very conserved query sequences. Conserved sequenced can generally **very deep alignments** which will require signifantly more memory. If you encounter out-of-memory errors during job execution, consider increasing the `request_memory` attribute in your submit file. You can also utilize the `retry_request_memory = <memory/expression>` command in your submit file to request a retry if the job holds for an out-of-memory error. For more information on how to use `retry_request_memory`, visit our [Request variable memory](https://chtc.cs.wisc.edu/uw-research-computing/variable-memory#use-retry_request_memory) documentation page.

```bash

4. Submit your data-pipeline jobs:

    ```
   condor_submit scripts/data_pipeline.sub
   ```

    > [!TIP]  
    > You can test the data pipeline using reduced-size databases by defining the `USE_SMALL_DB=1` variable when submitting your jobs. This is useful for debugging and testing purposes, as it reduces resource requirements and speeds up job execution. To use the small database set, submit your jobs with the following command: 
       > ```
       > condor_submit USE_SMALL_DB=1 data_pipeline.sub
       > ```

5. Track your job progress:

    ```
   condor_watch_q
   ```

### Submit Your AlphaFold3 Jobs - GPU-Accelerated Structural Prediction (Step 2)

[INSERT INSTRUCTIONS FOR SUBMITTING THE INFERENCE PIPELINE JOBS HERE]

[INSERT GRAPHIC OF ALPHAFOLD3 WORKFLOW WITH INFERENCE PIPELINE HIGHLIGHTED]

1. Change to your `tutorial-CHTC-AF3/` directory:
    ```bash
    cd ~/tutorial-CHTC-AF3/
   ```

2. Review your Data Pipeline executable script `scripts/inference_pipleine.sh`. Generally, no changes will be necessary. You should, however, review the following information and options below as your AF3 jobs may require additional non-default options.

    #### Overview: AlphaFold3 Inference Pipeline Executable

    This script implements the **structure inference portion** of the AlphaFold3 workflow. The inference pipeline loads model weights, runs the deep learning model on previously generated MSA/template features, and produces predicted 3D structures. It **does not** generate alignments or template features. This script is designed for execution inside HTCondor jobs and supports both containerized workflows and user-specified model weight files.

    The script handles:

   - Preparing a clean per-job working directory on the GPU execute point
   - Expanding `.data_pipeline.tar.gz` results from the previous data pipeline stage 
   - Detecting and loading AlphaFold3 model parameters (**User supplied**) 
   - Running AlphaFold3 **inference only**  
   - Packaging predicted structures and metadata  
   - Cleaning up the job sandbox to minimize disk use

    #### 1. Model Parameters and Containers

    Unlike the data pipeline, the inference stage does **not need the full AlphaFold3 databases**. It only requires:

    - A valid AlphaFold3 **model parameter file** (`.bin.zst`)
    - A valid AlphaFold3 **container image** (if running outside the container)

    The script allows users to specify a model file using:

    ```
    --model_param_file <name_of_model_weights_file.zst>
    ```

    If the model weights are compressed (`.zst`), the script automatically decompresses them into the working directory. 

    #### 2. Command-Line Options

    | Flag | Meaning |
    |------|---------|
    | `--work_dir_ext <name>` | Name appended to working directory (`work.<name>`) |
    | `--verbose` / `--silent` | Adjust verbosity |
    | `--no_copy` | Do not copy container or model parameters locally |
    | `--container <image>` | Path to Apptainer/Singularity image |
    | `--model_param_file <path>` | Path to AlphaFold3 model weights |
    | `--user_specified_alphafold_options` | Additional arguments passed directly to `run_alphafold.py` |
    | `--enable_unified_memory` | Enable unified GPU memory mode for large jobs |

    #### 3. Working Directory Setup

    Creates:

    ```
    work.<ext>/
        af_input/
        af_output/
        models/
        public_databases/
    ```

    The inference stage expects the output from the data pipeline in the form of:

    ```
    *.data_pipeline.tar.gz
    ```

    These archives are extracted into `af_input/`, reconstructing the full feature directory structure required by AlphaFold3.

    #### 4. Preparing Model Weights

    Valid model weights may be provided in either uncompressed or `.zst` form. The script handles both:

    - If compressed:  
      ```
      zstd --decompress > models/af3.bin
      ```
    - If uncompressed:  
      ```
      cp af3.bin models/
      ```

    #### 5. GPU Capability Handling

    The script inspects the GPU compute capability using:

    ```
    python -c "import jax; print(jax.local_devices(backend='gpu')[0].compute_capability)"
    ```

    For GPUs with compute capability **7.x**, it automatically:

    - Forces flash attention to `xla`
    - Sets required XLA flags to avoid runtime errors

    Users may also enable unified memory mode if working with especially large complexes.

    #### 6. Running the AlphaFold3 Inference Stage

    Inside the container, the script executes:

    ```
    python run_alphafold.py  --run_data_pipeline=false --run_inference=true --input_dir=/root/af_input --model_dir=/root/models --output_dir=/root/af_output
    ```

    Users may append custom AF3 flags through:

    ```
    --user_specified_alphafold_options "<options>"
    ```
   
    This allows for advanced customization of the inference run. It can be used to modify parameters such as:
    - Enabling different size model compilation bucket sizes (`--buckets 5132, 5280, 5342`)
    - Altering or Disabling JAX GPU Memory Preallocation (useful for larger complexes over 8k tokens, see our section on [Handling Large Complexes]()
    - Other AlphaFold3-specific options
    For more information on available AlphaFold3 options, refer to the [AlphaFold3 GitHub Repository](https://github.com/google-deepmind/alphafold3/)

    Output includes:

    - Predicted 3D structures  
    - Ranking scores  
    - Seed/sample-level predictions  

    #### 7. Packaging Results

    Each prediction directory under `af_output/` is archived:

    ```
    <target>.inference_pipeline.tar.gz
    ```

    These tarballs are returned to the submit host for downstream use.

    #### 8. Cleanup

    After packaging the results, the script removes:

    - `work.<ext>/`  
    - Shell history files (`.bash_history`, `.lesshst`, etc.)

    This ensures minimal disk usage on execute nodes and prepares the job sandbox for automatic cleanup by HTCondor.

3. Create your submit file `inference_pipeline.sub`. You will need to edit the `MODEL_WEIGHT_PATH` **and** `gpus_minimum_memory`. You can specify additional parameters for the exectuable in the `arguments` attribute as needed. Refer to the [AlphaFold3 Data Pipeline Executable - Command-Line Options]() section above for available options.

    ```bash
    container_image = file:///staging/groups/glbrc_alphafold/af3/alphafold3.minimal.22Jan2025.sif
    
    executable = inference_pipeline.sh
    
    environment = "myjobdir=$(directory)"
    
    MODEL_WEIGHTS_PATH = /staging/damorales4/af3/weights/af3.bin.zst
    
    log = ../logs/inference_pipeline.log
    output = inference_pipeline_$(Cluster)_$(Process).out
    error  = inference_pipeline_$(Cluster)_$(Process).err
    
    initialdir = $(directory)
    # transfer all files in the inference_inputs directory
    transfer_input_files = inference_inputs/, $(MODEL_WEIGHTS_PATH)
    
    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    
    request_memory = 8GB
    # need space for the container (3GB) as well
    request_disk = 10GB
    request_cpus = 4
    request_gpus = 1
    
    # we should be able to run short jobs on CUDA_CAPABILITY=7.x but need
    # other environment variables and options to be set. This is done automatically
    # in inference_pipeline.sh
    gpus_minimum_memory = 0
    
    # short jobs 4-6 hours so it is okay to use is_resumable
    +GPUJobLength = "short"
    +WantGPULab = true
    +is_resumable = true
    
    # Use --user-specified-alphafold-options to pass any extra options to AlphaFold3, such as
    # arguments = --model_param_file af3.bin.zst --work_dir_ext $(Cluster)_$(Process) --user-specified-alphafold-options "--buckets 5982"
    arguments = --model_param_file af3.bin.zst --work_dir_ext $(Cluster)_$(Process)
    
    queue directory from list_of_af3_jobs.txt   
   ```

This submit file will read the contents of `list_of_af3_jobs.txt`, iterate through each line, and assign the value of each line to the variable `$(directory)`. This allows you to programmatically submit _N_ jobs, where _N_ equals the number of AlphaFold3 job directories you previously created. Each job processes one AlphaFold3 job directory and uses the CHTC-maintained AlphaFold3 container image, which is transferred to the Execution Point (EP) by HTCondor.

The GPU-accelerated inference jobs do not require the full AlphaFold3 databases, so there is no need to check for pre-staged databases on the execute nodes. However, it does require the model weights to be transferred to the EP. The submit file above transfers the model weights specified in the `MODEL_WEIGHTS_PATH` variable to the EP for each job.

Beyond the model weights, your AlphaFold3 jobs may require different resource requests depending on the size and complexity of the proteins you are modeling. Typically, this mainly involves adjusting the `request_disk` and `gpus_minimum_memory` attributes in your submit file. Estimating the required GPU memory can be challenging, as it depends on factors such as sequence length, number of chains, and model parameters. As a starting point, you can refer to the following general guidelines:

| Tokens     | Estimated GPU Memory Requirement |
|------------|----------------------------------|
| Up to 1200 | 8-10 GB                          |
| 1200-1850  | 15-20 GB                         |
| 2000-3000  | 35-40 GB                         |
| Over 3000  | 70+ GB                           |

You can generally estimate number of tokens as approximately equivalent to 1.2x the sequence length. For example, a protein with 1000 amino acids would have roughly 1200 tokens. If you are modeling multimeric complexes, you will need to account for the combined sequence lengths of all chains. For structures that include nucleic acids (DNA/RNA), each nucleotide is considered one token. Post-translational modifications (PTMs) and additional ligands may also increase token counts, depending on the specific modifications involved. 

For very large complexes exceeding 10k tokens, you may need to enable unified memory mode using the `--enable_unified_memory` flag in the executable script. This allows AlphaFold3 to utilize system RAM in addition to GPU memory, which can help accommodate larger models. However, it may also lead to **significantly slower performance** due to increased data transfer times between system RAM and GPU memory. You will also need to ensure that the execute node has sufficient system RAM to support unified memory mode by increasing the `request_memory` attribute in your submit file. We recommend reaching out to the CHTC Research Computing Facilitation team for assistance with very large complexes.

    > [!TIP]  
    > CHTC has a number of smaller RTX series GPUs (e.g., RTX 4000, RTX 5000) with only 8-16GB of GPU memory. If your job requires less than 10GB of GPU memory, you can set `gpus_minimum_memory = 10000` in your submit file to allow your jobs to match to these smaller GPUs. This can help increase the number of available machines for your jobs, leading to faster start times.
    > 
    > If your jobs require more GPU memory than is available on these smaller GPUs, you can adjust the `gpus_minimum_memory` attribute in your submit file to request machines with larger GPUs (e.g., A100, V100). For example, to request machines with at least 32GB of GPU memory, you can set `gpus_minimum_memory = 32000` in your submit file. This will help ensure that your jobs are matched to machines with sufficient GPU memory to run successfully.

```bash

4. Submit your data-pipeline jobs:

    ```
   condor_submit scripts/inference_pipeline.sub
   ```

5. Track your job progress:

    ```
   condor_watch_q
   ```


### Post-Data Pipeline Data Wrangling

Running the data pipeline will generate a tarball for each job containing the necessary input files for the inference pipeline. These tarballs will be named `<job_name>.data_pipeline.tar.gz` and will be located in the respective job directories.

You should now have a directory of basecalled FASTQ and BAM files in your outputs folder. You'll likely want to perform additional steps after basecalling, such as checking read quality, mapping to a reference genome, or calling variants. We recommend merging your basecalled FASTQ files into a single file for downstream analysis. You can do this using the `cat` command:

```
for f in outputs/basecalledFASTQs/*.fastq; do
    cat "$f" >> outputs/merged_basecalled_reads.fastq
done
```

You can use this merged FASTQ file for running FastQC, mapping, or variant calling in the next sections. The recommended next step is to run FastQC to assess the quality of your basecalled reads. You can find a step-by-step tutorial on running [FastQC on the OSPool](https://portal.osg-htc.org/documentation/software_examples/bioinformatics/tutorial-fastqc/) on our documentation portal.

## Next Steps

Now that you've completed this long-read genomics tutorial on the OSPool, you're ready to adapt these workflows for your own data and research questions. Here are some suggestions for what you can do next:

🧬 Apply the Workflow to Your Own Data
* Replace the tutorial datasets with your own POD5 files and reference genome.
* Modify the basecalling, mapping, and variant calling submit files to fit your data size, read type (e.g., simplex vs. duplex), and resource needs.

🧰 Customize or Extend the Workflow
* Incorporate quality control steps (e.g., filtering or read statistics) using FastQC.
* Use other mappers or variant callers, such as ngmlr, pbsv, or cuteSV.
* Add downstream tools for annotation, comparison, or visualization (e.g., IGV, bedtools, SURVIVOR).

📦 Create Your Own Containers
* Extend the Apptainer containers used here with additional tools, reference data, or dependencies.
* For help with this, see our [Containers Guide](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers/).

🚀 Run Larger Analyses
* Submit thousands of basecalling or alignment jobs across the OSPool.
* Explore data staging best practices using the OSDF for large-scale genomics workflows.
* Consider using workflow managers (e.g., [DAGman](https://portal.osg-htc.org/documentation/htc_workloads/automated_workflows/dagman-workflows/) or [Pegasus](https://portal.osg-htc.org/documentation/htc_workloads/automated_workflows/tutorial-pegasus/)) with HTCondor.

🧑‍💻 Get Help or Collaborate
* Reach out to [support@osg-htc.org](mailto:support@osg-htc.org) for one-on-one help with scaling your research.
* Attend office hours or training sessions—see the [OSPool Help Page](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) for details.

### Software

In this tutorial, we created several *starter* Apptainer containers, including tools like: Dorado, SAMtools, Minimap, and Sniffles2. These containers can serve as a *jumping-off* for you if you need to install additional software for your workflows. 

Our recommendation for most users is to use "Apptainer" containers for deploying their software.
For instructions on how to build an Apptainer container, see our guide [Using Apptainer/Singularity Containers](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers-singularity/).
If you are familiar with Docker, or want to learn how to use Docker, see our guide [Using Docker Containers](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers-docker/).

This information can also be found in our guide [Using Software on the Open Science Pool](https://portal.osg-htc.org/documentation/htc_workloads/using_software/software-overview/).

### Data

The ecosystem for moving data to, from, and within the HTC system can be complex, especially if trying to work with large data (> gigabytes).
For guides on how data movement works on the HTC system, see our [Data Staging and Transfer to Jobs](https://portal.osg-htc.org/documentation/htc_workloads/managing_data/overview/) guides.

### GPUs

The OSPool has GPU nodes available for common use, like the ones used in this tutorial. If you would like to learn more about our GPU capacity, please visit our [GPU Guide on the OSPool Documentation Portal](https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/gpu-jobs/).

## Getting Help

The OSPool Research Computing Facilitators are here to help researchers using the OSPool for their research. We provide a broad swath of research facilitation services, including:

* **Web guides**: [OSPool Guides](https://portal.osg-htc.org/documentation/) - instructions and how-tos for using the OSPool and OSDF.
* **Email support**: get help within 1-2 business days by emailing [support@osg-htc.org](mailto:support@osg-htc.org).
* **Virtual office hours**: live discussions with facilitators - see the [Email, Office Hours, and 1-1 Meetings](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) page for current schedule.
* **One-on-one meetings**: dedicated meetings to help new users, groups get started on the system; email [support@osg-htc.org](mailto:support@osg-htc.org) to request a meeting.

This information, and more, is provided in our [Get Help](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) page.
