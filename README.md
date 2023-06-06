# SpecBM
Spectral Bundle Method - Primal and Dual formulations (SBMP and SBMD).

Note: the code requires the installation of [Sedumi](https://sedumi.ie.lehigh.edu/) and [Mosek](https://www.mosek.com/).

To get a quick start, try Experiment_1_1.m


To access large scale data, please visit [Google drive](https://drive.google.com/drive/folders/101KqJ56fwcZMuYuTTpwUASnevcnB2frt?usp=drive_link).


# Description
SBMP and SBMD consider the standard primal and dual vectorized semidefinite programs
```
		minimize 	$\langle C,X \rangle$						maximize 	b'y
	(1)	subject to	A(X) = b,				         (2)	subject to	A'y + z = c,	
				x \in PSD$$							z \in PSD
```

# Description - SBMP
SBMP solve the penalized primal problem
```
	minimize c'x 
```


# Description - SBMD
