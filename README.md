# TransformerBeta: Complementary Peptide Generation using Transformer Architecture

## Table of Contents
1. [Introduction](#introduction)
2. [Setup Guidance for Local Installation](#setup-guidance-for-local-installation)
    - [Create a Virtual Environment](#create-a-virtual-environment)
    - [Install the Required Packages and Dependencies](#install-the-required-packages-and-dependencies)
    - [Clone the Repository and Install the TransformerBeta Package](#clone-the-repository-and-install-the-transformerbeta-package)
3. [Usage](#usage)
4. [Model Architecture](#model-architecture)
5. [Results and Evaluation](#results-and-evaluation)
6. [Contribution](#contribution)
7. [License](#license)
8. [Acknowledgments](#acknowledgments)

## Introduction

TransformerBeta is a generative model based on Transformer architecture, developed to generate complementary binder for linear peptide epitopes of length 8 in an antiparallel beta strand conformation. The model is trained on a curated dataset of length 8 antiparallel beta strand pairs from the AF2 Beta Strand Database.

![TransformerBeta.png](data:image/png;base64,iVBO)

This page provides instructions for local installation and usage of TransformerBeta.

For detailed insight into TransformerBeta, we suggest referring to our research paper (yet to be released). 

To use TransformerBeta effortlessly without installation, visit: [TransformerBeta on Google Colab](https://colab.research.google.com/github/HZ3519/TransformerBeta/blob/main/notebooks/peptide_design_colab.ipynb)

For individual usage of the AF2 Beta Strand Database, check out: [AF2 Beta Strand Database on Huggingface](https://huggingface.co/datasets/hz3519/AF2_Beta_Strand_Database/tree/main).

For accessing the model weights and corresponding training, validation, and test sets mentioned in the paper, refer to: [TransformerBeta on Huggingface](https://huggingface.co/hz3519/TransformerBeta).

## Setup Guidance of local installation

This section provides instructions to setup the TransformerBeta project on your local machine.

#### Create a virtual environment

Creating a virtual environment is recommended as it helps avoid conflicts between package dependencies. The below command can be used for creating a virtual environment in Python:

python3 -m venv [name_of_virtual_environment]

To activate your virtual environment, you can use the following commands based on your operating system:

- For Windows: [name_of_virtual_environment]\Scripts\activate
- For Unix or MacOS: source [name_of_virtual_environment]/bin/activate

#### Install the dependencies and TransformerBeta package

With the virtual environment activated, we can install the necessary packages and dependencies. Use the following commands to do so:

git clone https://github.com/HZ3519/TransformerBeta_project.git
pip install d2l==0.17.5 --no-deps
pip install -r ./TransformerBeta_project/requirements.txt
pip install ./TransformerBeta_project

## Usage
Include a few examples of how to use your project. This would typically include code snippets and screen shots.


## Results and Evaluation
Discuss the results of the project and any evaluation metrics that you used to measure its success.

## Contribution
Outline the ways in which someone might contribute to your project, the process for submitting pull requests to you, and how they should report bugs or issues they find.

## License

The TransformerBeta models and AlphaFold 2 Beta Strand Database is made available under the terms of the MIT License.

## Citing our research

