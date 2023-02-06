phiquadmv paper scripts
========================

Matlab scripts to reproduce the results in [1].

**Note:** this repo contains the scripts the way they were used to produce the results in [1].
While the scripts in this repo are self-sufficient, the software developed in [1] is mantained in a separate repo, see [2].

Dependencies
------------

* Matlab >= r2021b
* Chebfun (available for free at https://www.chebfun.org/)

Installation
-------------

Just make sure that chebfun is in the Matlab path.

Usage
-----

In order to produce the numerical experiments from [1], run the following scripts:
* Problem 1 - Heat equation in 3D --> Heat_3D.m
* Problem 2 - Advection-di usion problem with a Sishkin mesh --> AdvDiff_2D.m
* Problem 3 - Hochbruck-Ostermann equation --> HochOster_2D_gauss.m (phiquadmv and Gauss quadrature), 
                                               HochOster_2D_cheb.m (phiquadmv and Clenshaw-Curtis quadrature) and HochOster_2D_expmv.m (expmv routine)

References
----------

[1] M. Croci & J. Muñoz-Matute, Exploiting Kronecker structure in exponential integrators: fast approximation of the action of $φ$-functions of matrices via quadrature, ArXiV preprint (https://arxiv.org/abs/2211.00696), 2022.
[2] https://github.com/jmunoz022/phiquadmv

Acknowledgements
----------------

This material is based upon work supported by the Department of Energy, NNSA under Award Number DE-NA0003969 also by the European Union’s Horizon 2020 research and innovation program under the Marie Sklodowska-Curie individual fellowship No. 101017984 (GEODPG).
