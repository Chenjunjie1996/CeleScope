version 1.0

import "../tools/structs.wdl"

task count {
    input {
        String sample_name
        File in_bam
        String genomeDir
        File in_data
        Int mem_on_bam

        Runtime runtime_count
    }

    runtime {
        cpu: runtime_count.cpu
        memory:  if mem_on_bam > runtime_count.memory_gb then mem_on_bam + "GiB" else runtime_count.memory_gb + "GiB"
        docker: runtime_count.docker
        queue: runtime_count.queue
    }

    command {
        set -euo pipefail
        mv "~{in_data}" ".data.json"
        celescope rna count \
        --outdir "05.count" --sample "~{sample_name}" --assay rna \
        --bam "~{in_bam}" \
        --genomeDir "~{genomeDir}" \
        --thread "~{runtime_count.cpu}" 

        wc -l "05.count/~{sample_name}_matrix_10X/barcodes.tsv" | cut -f 1 -d ' ' > "cell_num.txt"
        wc -l "05.count/~{sample_name}_matrix_10X/genes.tsv" | cut -f 1 -d ' ' > "gene_num.txt"
        tar cf "~{sample_name}_matrix_10X.tar" --directory="05.count/~{sample_name}_matrix_10X" .
    }

    output {
        File out_data = ".data.json"
        File out_matrix = "~{sample_name}_matrix_10X.tar"
        Int cell_num = read_int("cell_num.txt")
        Int gene_num = read_int("gene_num.txt")
        Int mem_on_mtx = ceil(cell_num * gene_num * 0.00000003 + 2)
    }
}
