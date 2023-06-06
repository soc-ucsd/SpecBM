# SpecBM
Spectral Bundle Method - Primal and Dual formulations (SBMP and SBMD).

Note: the code requires the installation of [Sedumi](https://sedumi.ie.lehigh.edu/) and [Mosek](https://www.mosek.com/).

To get a quick start, try Experiment_1_1.m


To access large scale data, please visit [Google drive](https://drive.google.com/drive/folders/101KqJ56fwcZMuYuTTpwUASnevcnB2frt?usp=drive_link).

\newcommand{\Trace}{\mathop{\bf tr}}
\newcommand{\tr}{\mathsf{T}}

# Description
SBMP and SBMD consider the standard primal and dual vectorized semidefinite programs
```math 
	\min_{X}\quad \langle C,X \rangle, \quad \mathrm{subject~to}\quad \mathcal{A}(X) = b,\; X \in \mathbb{S}^n_+. \quad [P]
```
```math
	\min_{y}\quad b^{\mathsf{T}}y, \quad \mathrm{subject~to}\quad C-\mathcal{A}^{*}y(X) = Z,\; Z \in \mathbb{S}^n_+. [D]
```

# Description - SBMP
SBMP solves the penalized primal problem 
```math
\min_{X \in \mathcal{X}_0} \quad \langle C,X\rangle + \rho \max \{\lambda_{\max}(-X),0\},
```
where $` \mathcal{X}_0 =\{X \in \mathbb{S}^n_+ \mid \mathcal{A}(X) = b\} `$.

The parameter $` \rho `$ should be chosen as $$\rho > \sup_{Z^{\star} \in \mathcal{D}^\star} \mathop{\bf tr}(Z^{\star}),$$ where $` \mathcal{D}^\star = \left\{(y,Z) \in \mathbb{R}^m \times \mathbb{S}^{n} \mid d^\star = b^{\mathsf{T}} y, Z+\mathcal{A}^* (y) = C, Z \in \mathbb{S}^n_+\right\}`$ is the optimal solution set of the dual problem (2).


# Description - SBMD
SBMD solves the penalized dual problem 
```math
\min_{y \in \mathbb{R}^m} \quad -b^{\mathsf{T}} y + \rho \max \{\lambda_{\max}(\mathcal{A}^{*}y-C),0\}.
```
The parameter $` \rho `$ should be chosen as $$\rho > \sup_{X^{\star} \in \mathcal{P}^\star} \mathop{\bf tr}(X^{\star}),$$
where $` \mathcal{P}^\star= \left\{X \in \mathbb{S}^{n} \mid p^\star = \langle C, X\rangle, \mathcal{A}(X) = b, X \in \mathbb{S}^n_+\right\}`$ is the optimal solution set of the primal problem (1).

