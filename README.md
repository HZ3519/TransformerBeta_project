# TransformerBeta

![TransformerBeta](./media/Figure_strategy.png)

paper: [Computational design of target-specific linear peptide binders with TransformerBeta](https://arxiv.org/abs/2410.16302)

## Introduction

TransformerBeta is a generative model based on Transformer architecture, developed to generate complementary binder for linear peptide epitopes of length 8 in an antiparallel beta strand conformation. The model is trained on a curated dataset of length 8 antiparallel beta strand pairs from the AF2 Beta Strand Database.

This page provides instructions for local installation and usage of TransformerBeta.

For detailed insight into TransformerBeta, we suggest referring to our research paper (yet to be released). 

To use TransformerBeta effortlessly without installation, visit: [TransformerBeta on Google Colab](https://colab.research.google.com/github/HZ3519/TransformerBeta/blob/main/notebooks/peptide_design_colab.ipynb)

For individual usage of the AF2 Beta Strand Database, check out: [AF2 Beta Strand Database on Huggingface](https://huggingface.co/datasets/hz3519/AF2_Beta_Strand_Database/tree/main).

For accessing the model weights and corresponding training, validation, and test sets mentioned in the paper, refer to: [TransformerBeta on Huggingface](https://huggingface.co/hz3519/TransformerBeta).

## Setup Guidance of local installation

This section provides instructions to setup the TransformerBeta project on your local machine.

#### Create a virtual environment

Creating a virtual environment is recommended because it helps avoid conflicts between package dependencies. Use the command below to create a virtual environment in Python:

```bash
python3 -m venv [name_of_virtual_environment]
```

To activate your virtual environment, use the appropriate command for your operating system and shell:

- For Windows Command Prompt:
```bash
[name_of_virtual_environment]\Scripts\activate.bat
```
- For Windows PowerShell:
```bash
[name_of_virtual_environment]\Scripts\Activate.ps1
```

- For Unix or MacOS using Bash or Zsh:
```bash
source [name_of_virtual_environment]/bin/activate
```

If you're using a different shell, the command may vary. Note that on Windows, if you're using the Windows Subsystem for Linux or Git Bash, you can use the Unix-style command.

#### Install the dependencies and TransformerBeta package

With the virtual environment activated, we can install the necessary packages and dependencies. First, clone the TransformerBeta repository and install the required packages using the following commands:

```bash
git clone https://github.com/HZ3519/TransformerBeta_project.git
pip install d2l==0.17.5 --no-deps
pip install -r ./TransformerBeta_project/requirements.txt
pip install ./TransformerBeta_project
```

## Usage
To use the model, please refer to our guide, `TransformerBeta_project/notebooks/peptide_design_local.ipynb`. This notebook provides a detailed walkthrough of how to generate complementary peptide binders for a given peptide epitope. Various functionalities are made available to aid in the generation of complementary peptides for your epitope of interest.

## Contribution
Contributions to this project are always welcome! If you come across any bugs or issues, feel free to open an issue or submit a pull request. We appreciate your efforts in improving our project. For directly reporting bugs or issues, you can also reach out to:

Haowen Zhao at hz362@cam.ac.uk

## License

The TransformerBeta models and AlphaFold 2 Beta Strand Database is made available under the terms of the MIT License.

## Citing our research

If you use our work, please cite:

```bibtex
@article{zhao2024computational,
  title={Computational design of target-specific linear peptide binders with TransformerBeta},
  author={Zhao, Haowen and Aprile, Francesco A and Bravi, Barbara},
  journal={arXiv preprint arXiv:2410.16302},
  year={2024}
}
```
