# All_Supported_Nd_Vectors_BMCIF

## Overview
This repository provides a Python implementation of the algorithms for determining all Supported Nondominated Vectors in Bi-objective Minimum Cost Integer Flow (BMCIF) problems. The methods implemented are based on decision-space and/or outcome-space approaches presented in the following paper:

**Reference:**  
Könen, D., & Stiglmayr, M. (2023). *Output-sensitive Complexity of Multi-Objective Integer Network Flow Problems*. arXiv preprint arXiv:2312.01786.  
[View on arXiv](https://arxiv.org/abs/2312.01786)

## Repository Structure
The repository contains the following folders:
- **`src/`**: Contains Python source code implementing the algorithms.
- **`instances/`**: Includes test instances used for benchmarking.
- **`output/`**: Stores the results generated by the algorithm. 

## Features
- Implements decision-space and outcome-space methods for solving BMCIF problems.
- Supports various test instances for validation and performance evaluation.
- Outputs the number supported nondominated vectors (and extreme supported vectors) and needed time for each of the different methods.

## Installation & Usage
### Requirements
Ensure you have Python installed (>=3.8) and the necessary dependencies

### Running the Algorithm
To execute the algorithm, use the following command: (To do!)

python src/main.py --input instances/example_instance.txt --output output/results.txt (Verbessern!)

### Example Usage (To do!)

python src/main.py --


## Citation
If you use this implementation in your research, please cite:


@article{konen2023output,
  title={Output-sensitive Complexity of Multi-Objective Integer Network Flow Problems},
  author={K{\"o}nen, D. and Stiglmayr, M.},
  journal={arXiv preprint arXiv:2312.01786},
  year={2023}
}

## License

This repository is licensed under the MIT License. See LICENSE for details.

## Contact

For questions or issues, please contact koenen@uni-wuppertal.de or open an issue in this repository.




