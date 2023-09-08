#!/usr/bin/env bash

# Ensure the samplesheet.tsv file is empty or create it
> samplesheet.tsv

# Loop through all fastq.gz and fq.gz files
for file in *.{fastq,fq}.gz; do
    # Check if the file actually exists
    if [[ ! -e $file ]]; then
        continue
    fi

    echo "Processing file: $file"

        # Skip if file is an R2 or _2 file
    if [[ $file =~ "_R2.fastq.gz"  || $file =~ "_R2.fq.gz" || $file =~ "_2_val_2.fastq.gz" || $file =~ "_2_val_2.fq.gz" || $file =~ "_2.fastq.gz" || $file =~ "_2.fq.gz" ]]; then
        continue
    fi


    paired_file=""
    # Determine if the file is R1, _1_val_1, or single-end and set paired_file accordingly
    if [[ $file =~ "_R1.fastq.gz" || $file =~ "_R1.fq.gz" ]]; then
         id=$(echo $file | sed 's/_R1.*//')
        paired_file="${id}_R2.${file#*.}"
    elif [[ $file =~ "_1_val_1.fastq.gz" || $file =~ "_1_val_1.fq.gz" ]]; then
        id=$(echo $file | sed 's/_1_val_1.*//')
        paired_file="${id}_2_val_2.${file#*.}"
    elif [[ $file =~ "_1.fastq.gz" || $file =~ "_1.fq.gz" ]]; then
        id=$(echo $file | sed 's/_1.*//')
        paired_file="${id}_2.${file#*.}"      
    else
        id=$(basename "$file" | rev | cut -d "." -f 3- | rev)
    fi

    # Get absolute path for the file
    abs_file=$(realpath "$file")

    # If it's a paired-end read, check if paired_file exists and get its absolute path
    if [[ -n $paired_file && -e $paired_file ]]; then
        abs_paired=$(realpath "$paired_file")
        echo -e "$id\t$id\t$abs_file\t$abs_paired" >> samplesheet.tsv
    elif [[ -z $paired_file ]]; then
        # For single-end reads, only output the single file path
        echo -e "$id\t$id\t$abs_file" >> samplesheet.tsv
    fi
done
  