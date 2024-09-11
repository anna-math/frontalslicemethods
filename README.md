
# Frontal Slice Methods for T-product Systems

**Content**
This is the code repository for the research publication "Tensor Frontal Slice  Sketching Approaches to Large-Scale Linear Systems". 
The manuscript of this paper can be accessed at https://arxiv.org/abs/1908.08864. 

 - We provided a set of illustrative code of sketching methods for solving t-product tensor linear system that serves as a proof of concept, and also a set of robust code that can be executed for accompanying datasets.
 -  Reproducible MATLab Code
	 - fullfrontalDescent.m, implements the full frontal descent algorithm for solving tensor linear systems $\mathcal{A} * \mathcal{X} = \mathcal{B}$ defined by tensor t-product.
 - General Usage (fullfrontalDescent.m implements a handy wrapper)
	 `[error, time, X_star] = frontalDescentAll(A, X, B, 'cyclic', paras);`
	 - This solves the system $\mathcal{A} * \mathcal{X} = \mathcal{B}$ supplied as `A,X,B` in the parameters.
- Inputs
	- AA: 3rd order measurement tensor of size n1 x n2 x nn
	- XX: 3rd order solution tensor of size  n2 x n3 x nn
	- BB: 3rd order resulting tensor of size  n1 x n3 x nn
	- variations: 'cyclic', 'leverageSampling', 'randomSampling',  'full' are supported 
	- paras:
		- .maxiter: maximum number of iterations
		- .alpha: learning rate
		- .stepsize: step size for all variations (default 1)
		- .blocksize: number of frontal slices to use at each iteration, default 1
		- .bandsize: effective number of frontal slices, default nn
		- .tolerance: tolerance for early stopping (default 0)
- Outputs
	- err: (paras.maxiter)x1 vector containing approximation error at each iteration
	- timeVec: (paras.maxiter)x1 vector containing elapsed cputime time at each iteration
	- Xt: approximation of solution to system

**Abstract**
Inspired by the row and column action methods for solving large-scale linear systems, in this work, we explore the use of frontal slices for solving tensor linear systems. In particular, this paper presents a novel approach for using frontal slices of a tensor $\mathcal{A}$ to solve tensor linear systems $\mathcal{A} * \mathcal{X} = \mathcal{B}$ where $*$ denotes the t-product. In addition, we consider variations of this method, including cyclic, block, and randomized approaches, each designed to optimize performance in different operational contexts. Our primary contribution lies in the development and convergence analysis of these methods. Experimental results on synthetically generated and real-world data, including applications such as image and video deblurring, demonstrate the efficacy of our proposed approaches and validate our theoretical findings.

**Citation**
We provided both iPynb illustrative code, Python production code for reproducible and experimental purposes under [LICENSE](https://github.com/hrluo/TensorDecisionTreeRegressor/blob/master/LICENSE).
Please cite our paper using following BibTeX item:

    @article{2024tensorsketch,
        title={Tensor Frontal Slice Sketching Approaches to Large-Scale Linear Systems},
        author={Hengrui Luo, Anna Ma},
        year={2024},
        eprint={https://arxiv.org/abs/2408.13547},
        archivePrefix={arXiv},
        primaryClass={math.LA}
    }

Thank you again for the interest and please reach out if you have further questions.
