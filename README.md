# Arnoldi with MATLAB
![matlab-badge](https://img.shields.io/badge/MATLAB-Solutions-blue
)

## Table of Contents
 - [About the Project](#about-the-project)
    - [Key Features](#key-features)
 - [Environment Setup](#environment-setup)
 - [Author](#author)

## About the project
This project introduces two MATLAB functions designed to extract a few eigenpairs from large-scale matrices, focusing on efficiency and accuracy. The functions are:

`myarnoldi.m`: Implements the Arnoldi method for eigenvalue computation.\
`myarnoldiro.m`: Enhances the Arnoldi method by incorporating reorthogonalization to improve numerical stability and accuracy.

### Key Features
- __Sparse Matrix Support:__ Utilizes MATLAB's sparse matrix capabilities to efficiently handle large-scale problems.
- __Demonstration with west0479 Matrix:__ The correctness of the implemented methods is verified using west0479, a MATLAB sparse demonstration matrix, showcasing the functions' ability to accurately compute eigenvalues and eigenvectors.
- __Residuals Computation:__ Calculates residuals to ensure the accuracy of the computed eigenpairs, with expectations of achieving machine precision for small iteration counts (e.g., 20, 40).
- __Ritz Values Computation:__ Involves computing and plotting Ritz values against the exact eigenvalues of the matrix to visually verify the effectiveness of the Arnoldi processes.
- __Error Analysis:__ Employs relative error analysis for the top three largest eigenvalues, illustrating the convergence behavior and the benefits of reorthogonalization through iteration-based plots.
- __Large Matrix Testing:__ Extends the validation to large matrices from the "SuiteSparse Matrix Collection," demonstrating the methods' scalability and efficiency on real-world datasets.
Shift-and-Invert Spectral Transformation: Incorporates a shift-and-invert strategy in the reorthogonalized Arnoldi method to target eigenvalues near a specified value, enhancing the method's utility for specific applications.

## Environment Setup
To successfully run the provided MATLAB functions and test the algorithms with your matrices, follow these setup instructions carefully. This guide assumes you have MATLAB installed. If not, please visit the official MATLAB website to download and install MATLAB before proceeding.

The project uses west0479, a sparse demonstration matrix included in MATLAB, as a primary example for testing. To load this matrix into your workspace, execute the following commands in the MATLAB Command Window:
```matlab
load west0479;
A = west0479;
```
## Author
Timothy Hwang