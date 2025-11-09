sudo apt update
sudo apt install -y openjdk-21-jdk
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version
nextflow run main.nf -with-dockerdocke