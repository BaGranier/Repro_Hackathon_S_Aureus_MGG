nextflow.enable.dsl=2

/*
 * MASTER PIPELINE
 * 1️⃣ Build Docker SRA Toolkit
 * 2️⃣ Lancer le pipeline SRA (main_SRA.nf)
 * 3️⃣ Lancer le script R
 */

// ===========================
// Définition des channels
// ===========================
Channel.value(true)
    .set { trigger_build }

// ===========================
// Étape 1 : Build Docker
// ===========================
process BUILD_DOCKER {
    container 'ubuntu:22.04'

    input:
        val x from trigger_build

    output:
        val true into docker_built

    script:
    """
    echo "=== [1] Build de l'image Docker SRA Toolkit ==="
    docker build -t local/sra_toolkit_docker -f Dockerfile_SRA .
    echo "Docker build terminé."
    """
}

// ===========================
// Étape 2 : Lancer le pipeline SRA
// ===========================
process RUN_SRA {
    container 'local/sra_toolkit_docker'

    input:
        val y from docker_built

    output:
        val true into sra_done

    script:
    """
    echo "=== [2] Lancement du pipeline main_SRA.nf ==="
    nextflow run main_SRA.nf -with-docker
    echo "Pipeline main_SRA terminé."
    """
}

// ===========================
// Étape 3 : Lancer le script R
// ===========================
process RUN_R {
    container 'rocker/r-base:3.4.1'

    input:
        val z from sra_done

    output:
        path "output_R.txt"

    script:
    """
    echo "=== [3] Étape R ===" > output_R.txt
    R --version >> output_R.txt
    echo "Hello world from R!" >> output_R.txt
    """
}

// ===========================
// Workflow global
// ===========================
workflow {
    BUILD_DOCKER()
    RUN_SRA(BUILD_DOCKER.out)
    RUN_R(RUN_SRA.out)
}
