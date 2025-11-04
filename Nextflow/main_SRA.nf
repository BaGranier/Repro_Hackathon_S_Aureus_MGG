process TEST_SRA {
    container 'local/sra_toolkit_docker'

    output:
        path "output.txt"

    script:
    """
    echo "=== Test container SRA ===" > output.txt
    fastq-dump --version >> output.txt 2>&1
    echo "Hello from SRA Toolkit!" >> output.txt
    """
}

workflow {
    TEST_SRA()
}
