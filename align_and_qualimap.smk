import pandas as pd

run_ids_df = pd.read_csv("not_prealigned_metadata.tsv", sep="\t", usecols=["Run", "Experiment"])
run_ids_df['Experiment'] = run_ids_df['Experiment'].str.strip()  # Strip whitespace from Experiment identifiers
run_ids_df['Run'] = run_ids_df['Run'].str.strip()  # Strip whitespace from Run identifiers
run_ids = run_ids_df["Run"].tolist()  # Extract run IDs as a list


# Process experiments and runs
experiment_to_run = {}
experiment_to_runs = {}


# Assuming 'run_ids_df' is your DataFrame from which you populate 'experiment_to_run'
for experiment, runs in run_ids_df.groupby("Experiment")["Run"]:
    runs_list = runs.tolist()
    clean_experiment = experiment.strip()  # Remove whitespace from experiment ID
    clean_runs_list = [run.strip() for run in runs_list]  # Remove whitespace from each run ID
    if len(clean_runs_list) == 1:
        experiment_to_run[clean_experiment] = clean_runs_list[0]
    else:
        experiment_to_runs[clean_experiment] = clean_runs_list

experiment_to_runs["EXTEST"] = ["TEST2", "TEST2"]

# Count experiments based on the number of runs
single_run_experiments = sum(len(runs) == 1 for runs in experiment_to_runs.values())
two_run_experiments = sum(len(runs) == 2 for runs in experiment_to_runs.values())
three_run_experiments = sum(len(runs) == 3 for runs in experiment_to_runs.values())
four_run_experiments = sum(len(runs) == 4 for runs in experiment_to_runs.values())
more_than_four_run_experiments = sum(len(runs) > 4 for runs in experiment_to_runs.values())

# Print counts of experiments based on the number of associated runs
print(f"Number of experiments with a single run: {single_run_experiments}")
print(f"Number of experiments with 2 runs: {two_run_experiments}")
print(f"Number of experiments with 3 runs: {three_run_experiments}")
print(f"Number of experiments with 4 runs: {four_run_experiments}")
print(f"Number of experiments with more than 4 runs: {more_than_four_run_experiments}")


def get_bam_path_for_experiment(experiment):
    """
    Construct the file path for a BAM file based on the experiment name.
    Ensures that any whitespace in the experiment name or in the dictionary lookup is handled.
    """
    # Ensure we strip any whitespace from the experiment name before lookup
    clean_experiment = experiment.strip().replace(" ", "")
    
    # Lookup the run ID in the dictionary, ensuring to strip any potential whitespace
    run_id = experiment_to_run[clean_experiment].strip().replace(" ", "")
    
    # Construct and return the file path
    file_path = f"aligned/single/{run_id}.bam"

    # log all vars
    print(f"experiment: `{experiment}`")
    print(f"clean_experiment: `{clean_experiment}`")
    print(f"run_id: `{run_id}`")


    return file_path



# Reference genome information (not directly used in the Snakemake rules below but could be part of the workflow)
REF_GENOME = "canFam3"
print(f"processing {len(run_ids)} runs")

# Snakemake rules
rule all:
    input:
        "canFam3.fa",
        "canFam3.fa.pac",
        expand("fasta/{run_id}.fasta", run_id=run_ids),
        expand("aligned/single/{run_id}.bam", run_id=run_ids),
        expand("aligned/merged/{experiment}.bam", experiment=experiment_to_runs.keys()),
        expand("aligned/{experiment}.bam", experiment=experiment_to_run.keys())

rule download_genome:
    output: "canFam3.fa"
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz
        gunzip canFam3.fa.gz
        """

rule index_genome:
    input: "canFam3.fa"
    output: "canFam3.fa.pac", "canFam3.fa.ann", "canFam3.fa.amb"
    shell:
        """
        bwa index -a bwtsw {input}
        """

rule download_sra:
    output: "sra/{run_id}.sra"
    shell:
        r"""
        mkdir -p sra
        # Directly use Snakemake's wildcard in the shell commands

        # Attempt to download .sralite format first
        echo "Attempting to download .sralite format for {wildcards.run_id}..."
        if ! aws s3 cp --no-sign-request s3://sra-pub-zq-3/{wildcards.run_id}/{wildcards.run_id}.lite.1 sra/{wildcards.run_id}.sra 2>/dev/null; then
            echo ".sralite format download failed for {wildcards.run_id}. Attempting direct download in .sra format..."
            # If .sralite download fails, attempt to download .sra format
            if ! aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/{wildcards.run_id}/{wildcards.run_id} sra/{wildcards.run_id}.sra 2>/dev/null; then
                echo "Failed to download {wildcards.run_id} in any recognized format."
                exit 1
            fi
        fi
        echo "Successfully downloaded {wildcards.run_id} as .sra"
        """


rule dump_to_fasta:
    priority: 100
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 15 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        partition = "bmm"
    input:
        sra_file = "sra/{run_id}.sra"
    output:
        fasta = "fasta/{run_id}.fasta",
    params:
        outdir = "fasta",
    shell:
        """
        fasterq-dump  --skip-technical \
        --fasta-unsorted --threads {threads} --bufsize 1000MB --curcache 10000MB --mem {resources.mem_mb} \
        --outdir {params.outdir} \
        {input.sra_file}
        """

ruleorder: bwa_align > merge_bams_for_multi_run_experiment

rule bwa_align:
    priority: 200
    threads: 12  # Allows dynamic adjustment based on available resources
    resources:
        mem_mb=lambda wildcards, attempt: 32 * 1024 * attempt,  # Memory allocated dynamically based on attempt
        time=lambda wildcards, attempt: 8 * 60 * attempt,  # Time allocated dynamically based on attempt
        runtime=lambda wildcards, attempt: 8 * 60 * attempt,  # Runtime allocated dynamically based on attempt
        partition="bmm"  # Specify the partition for job scheduling if needed
    # conda: "env.yaml"  # Specifies the Conda environment file to use
    input:
        fasta="fasta/{run_id}.fasta",
        ref_index="canFam3.fa.pac",  # Reference index file for BWA
        ref_genome="canFam3.fa"  # Reference genome file
    output:
        aligned="aligned/single/{run_id}.bam"  # Output BAM file path
    shell:
        """
        bwa mem -t {threads} {input.ref_genome} {input.fasta} | \
        samtools sort -@ {threads} -o {output.aligned} - && \
        samtools index {output.aligned}
        """


rule link_bam_for_single_run_experiment:
    priority: 250
    input:
        bam=lambda wildcards: get_bam_path_for_experiment(wildcards.experiment)
    output:
        bam="aligned/{experiment}.bam"
    shell:
        """
        BASEDIR=$(pwd)
        ln -sf ${{BASEDIR}}/{input.bam} ${{BASEDIR}}/{output.bam}
        samtools index ${{BASEDIR}}/{output.bam}
        """


rule merge_bams_for_multi_run_experiment:
    priority: 250
    input:
        bams=lambda wildcards: ["aligned/single/" + run_id + ".bam" for run_id in experiment_to_runs.get(wildcards.experiment, [])]
    output:
        bam="aligned/merged/{experiment}.bam"
    params:
        num_bams=lambda wildcards: len(["aligned/single/" + run_id + ".bam" for run_id in experiment_to_runs.get(wildcards.experiment, [])]),
        bam_list=lambda wildcards: " ".join(["aligned/single/" + run_id + ".bam" for run_id in experiment_to_runs.get(wildcards.experiment, [])])
    shell:
        """
        echo "Merging BAMs for experiment {wildcards.experiment}..."
        echo "Number of BAMs to merge: '{params.num_bams}'"
        NUM_BAMS="{params.num_bams}"
        if [ $NUM_BAMS -gt 1 ]; then
            samtools merge -@ {{threads}} {output.bam} $(echo "{params.bam_list}")
        else
            echo "Error: Attempted to merge BAMs for an experiment with less than 2 runs."
            exit 1
        fi
        samtools index {output.bam}
        """