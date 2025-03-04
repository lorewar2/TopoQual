# TopoQual

TopoQual polishes Pacific Biosciences (PacBio) Circular Consensus Sequencing (CCS) data and accurately predicts quality scores.

## Setting up / Requirements

Rust 1.65.0+ and Python 3.9+ should be installed.
The following programs/packages are required for topoqual to run properly,

*   [`samtools`](http://www.htslib.org/) - to read the subread bam file.
*   [`pytorch`](https://pytorch.org/) - to evalute the quality scores using the deeplearning model.
*   [`numpy`](https://numpy.org/) - for computation.
*   [`pysam`](https://github.com/PacificBiosciences/actc) - for writing the modified bam file in python.

you can install samtools, pytorch using conda and numpy, pysam using pip:

```bash
conda install bioconda::samtools
conda install pytorch::pytorch torchvision torchaudio -c pytorch
pip install pysam
pip install numpy
```

## Usage

Download the repository.

```bash
git clone https://github.com/lorewar2/TopoQual.git
```

Configure the thread count in script.sh (Decrease/Increase the thread count depending on the memory availability, 1 thread requires ~10GB of memory)

TEST DATA:

Run the test sample with Topoqual

```bash
bash script.sh
```

REAL DATA:

Modify input/ouput variables to point to your data in script.sh

Run the real sample with Topoqual

```bash
bash script.sh
```

## How to cite

If you are using TopoQual in your work, please cite:

[TopoQual polishes circular consensus sequencing data and accurately predicts quality scores](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-06020-0)
