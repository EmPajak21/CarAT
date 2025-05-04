
# <img src="assets/carat_logo.png" alt="CarAT Logo" width="100" height="auto" style="vertical-align: middle;"> &mdash; Carbon Atom Tracker

[![Run Tests](https://github.com/EmPajak21/carat/actions/workflows/run-tests.yml/badge.svg)](https://github.com/EmPajak21/carat/actions/workflows/run-tests.yml/badge.svg)
[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


CarAT (Carbon Atom Tracker) is an automated tool that tracks biogenic carbon content (BCC) across complex industrial value chains. With the Together for Sustainability (TfS) consortium mandating BCC reporting by 2026 ([TfS Guidelines][reference1]), CarAT provides a scalable solution to this critical industry challenge, and also serves as a decision support tool for decarbonisation strategies by substituting fossil carbon with biogenic carbon.

## 🏭 Why the Chemicals Industry Needs CarAT
The chemicals industry is facing increasing pressure to reduce its carbon footprint and move toward sustainability. Accurately calculating BCC is challenging because:

- Products' BCC depends on all contributing feedstocks throughout the value chain
- Process efficiencies and variables affect carbon distribution
- Value chains are complex and interconnected
- Manual calculations are time-consuming and error-prone

CarAT addresses these challenges by leveraging existing Enterprise Resource Planning (ERP) data through a systematic three-stage approach.

## ⚙️ How CarAT Works

![CarAT Methodology](assets/methodology.png)
CarAT's methodology consists of three main stages:

1. **Data Preparation**: Processes value chain data from existing ERP systems.
2. **Atom Mapping**: Uses RXNMapper from [Schwaller et al.][reference2] to track atoms through chemical reactions.
3. **Optimisation**: Employs linear programming to calculate BCC across the entire value chain.

## 🚀 Getting Started  

To set up the environment for this project, we use an `environment.yml` file. This ensures that all dependencies are installed and configured correctly.  

### Setting Up the Environment  

1. **Create the Environment**  
    Run the following command in your terminal to create the environment:  
    ```bash  
    conda env create -f environment.yml  
    ```  

2. **Activate the Environment**  
    After the environment is created, activate it using:  
    ```bash  
    conda activate carat  
    ```


## 📈 Example Usage

For a complete walkthrough of how to use CarAT, please see the included Jupyter notebook:
[tdi_value_chain.ipynb](tdi_value_chain.ipynb)

This example demonstrates the full workflow including:
- Loading and visualising a value chain graph.
- Computing the bill of atoms with RXNMapper.
- Preprocessing data for LP formulation.
- Solving the LP optimisation model.
- Generating Sankey diagrams for visualisation of results.

*N.B. The linear program formulation created for CarAT can be found in `carat/core/linear_program.md`.*

[reference1]: https://www.tfs-initiative.com/app/uploads/2024/03/TfS_PCF_guidelines_2024_EN_pages-low.pdf "Together for Sustainability (TfS). (2024). Product Carbon Footprint (PCF) Guidelines. TfS Initiative."

[reference2]: https://www.science.org/doi/10.1126/sciadv.abe4166 "Schwaller, P., Hoover, B., Reymond, J.-L., Strobelt, H., & Laino, T. (2021). Extraction of organic chemistry grammar from unsupervised learning of chemical reactions. Science Advances, 7(15), eabe4166."
