process EXTRACT_KEY_VALUE {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.json"), emit: key_value_json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def extract_key_value(input_file):
        key_value_pairs = {}
        with open(input_file, 'r') as file:
            for line in file:
                parts = line.split('|')
                if len(parts) >= 2:
                    #value is the first part of each line separated by '|'.
                    value = parts[0].split()[-1]

                    #The key is the third to last field in each line separated by spaces, and then separated again by '|' to take the last part.
                    last_third_field = line.split()[-3]

                    key = last_third_field.split('|')[-1]

                    key_value_pairs[key] = value

        return key_value_pairs

    # Generate output JSON file name based on input file
    output_json_name = os.path.basename('${input_file}').replace('.tracking', '_keyvalue.json')

    extracted_data = extract_key_value('${input_file}')

    with open(output_json_name, 'w') as json_file:
        json.dump(extracted_data, json_file, indent=4)
    """
}

// workflow {
//     // Define the input file path
//     ch_input_file = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/null.tracking')

//     // Process the file
//     key_value_json = EXTRACT_KEY_VALUE(ch_input_file.map { [ [:], it ] })
// }
