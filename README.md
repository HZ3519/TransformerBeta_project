# TransformerBeta: Complementary Peptide Generation using Transformer Architecture

## Table of Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
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

#### Create a virtual environment

#### Install the required packages and dependencies

pip install d2l==0.17.5 --no-deps
pip install -r ./TransformerBeta_project/requirements.txt

#### Clone the repository and install the TransformerBeta package




## Usage
Include a few examples of how to use your project. This would typically include code snippets and screen shots.


## Results and Evaluation
Discuss the results of the project and any evaluation metrics that you used to measure its success.

## Contribution
Outline the ways in which someone might contribute to your project, the process for submitting pull requests to you, and how they should report bugs or issues they find.

## License

The TransformerBeta models and AlphaFold 2 Beta Strand Database is made available under the terms of the MIT License.

## Citing our research

