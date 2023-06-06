# SpecBM
Spectral Bundle Method - Primal and Dual formulations (SBMP and SBMD).

Note: the code requires the installation of [Sedumi](https://sedumi.ie.lehigh.edu/) and [Mosek](https://www.mosek.com/).

To get a quick start, try Experiment_1_1.m


To access large scale data, please visit [Google drive](https://drive.google.com/drive/folders/101KqJ56fwcZMuYuTTpwUASnevcnB2frt?usp=drive_link).

\newcommand{\Trace}{\mathop{\bf tr}}
\newcommand{\tr}{\mathsf{T}}

# Description
SBMP and SBMD consider the standard primal and dual vectorized semidefinite programs

		minimize 	<C,X>						maximize 	b'y
	(1)	subject to	A(X) = b,				  (2)	subject to	A'y + Z = c,	
				X \in PSD							Z \in PSD


# Description - SBMP
SBMP solves the penalized primal problem 
```math
\min_{X \in \mathcal{X}_0} \quad \langle C,X\rangle + \rho \max \{\lambda_{\max}(-X),0\},
```
where $` \mathcal{X}_0 =\{X \in \mathbb{S}^n_+ \mid \mathcal{A}(X) = b\} `$.
The parameter $`\rho`$ should be chosen as $` \rho > \sup_{Z^{\star} \in \mathcal{D}^\star} \mathop{\bf tr}(Z^{\star}) `$

# Description - SBMD
SBMD solves the penalized dual problem 
```math
\min_{y \in \mathbb{R}^m} \quad -b^{\mathsf{T}} y + \rho \max \{\lambda_{\max}(\mathcal{A}^{*}y-C),0\},
```
