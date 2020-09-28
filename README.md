# Fixed-Hessian-Logistic-Regression

This repo implements three methods of Privacy Preserving Logistic Regression Training using Homomorphic Encryption. The first is closely based on [1] and uses Gradient Descent minimisation. The second follows the approach of [2] and uses Nesterov Accelerated Gradient Descent with a compact database encoding. The third is closest to [3], and uses a fixed Hessian approximation to Newton-Raphson minimisation. All methods are built using the CKKS encryption scheme implemented in SEAL (https://github.com/Microsoft/SEAL). To run this code, you will first need to install SEAL.

[1] Kim, M., Song, Y., Wang, S., Xia, Y., & Jiang, X. (2018). Secure logistic regression based on homomorphic encryption: Design and evaluation. JMIR medical informatics, 6(2), e19.
[2] Kim, A., Song, Y., Kim, M., Lee, K., & Cheon, J. H. (2018). Logistic regression model training based on the approximate homomorphic encryption. BMC medical genomics, 11(4), 83.
[3] Bonte, C., & Vercauteren, F. (2018). Privacy-preserving logistic regression training. BMC medical genomics, 11(4), 86.
