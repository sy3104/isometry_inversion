# isometry_inversion

This repository contains companion code for the following paper:

- Satoshi Yoshida, Akihito Soeda, and Mio Murao, Universal construction of decoders from encoding black boxes, [Quantum 7, 957 (2023)][8], [arXiv:2110.00258][1].

The MATLAB codes are modified from Marco Túlio Quintino's [unitary_inverse][2], accompanied the following paper:

- Marco Túlio Quintino, Qingxiuxiong Dong, Atsushi Shimbo, Akihito Soeda, and Mio Murao, Reversing unknown quantum transformations: A universal protocol for inverting general unitary operations, [Phys. Rev. Lett. 123, 210502 (2019)][3], [arXiv:1810.06944][4].

## Requirement

These codes are written in MATLAB and require one of the following interpreters to conduct SDP:

- [CVX][5]: a Matlab-based convex modeling framework
- [YALMIP][6]: A Toolbox for Modeling and Optimization in MATLAB

These codes also use functions of QETLAB ([QETLAB][7]: A MATLAB Toolbox for Quantum Entanglement), but all used functions are contained in the subfolder [QETLAB_used_functions](https://github.com/sy3104/isometry_inversion/tree/main/QETLAB_used_functions).

It has been tested on MATLAB R2021b, CVX 2.2, and YALMIP R20210331.

## Description

The main scripts are contained in the subfolders [isometry_inversion_cvx](https://github.com/sy3104/isometry_inversion/tree/main/isometry_inversion_cvx) (The interpreter CVX is used to conduct SDP) and [isometry_inversion_yalmip](https://github.com/sy3104/isometry_inversion/tree/main/isometry_inversion_yalmip) (The interpreter YALMIP is used to conduct SDP).

- [run_isometry_inversion_cvx.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_inversion_cvx.m) / [run_isometry_inversion_yalmip.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_inversion_yalmip.m): Code that obtains the maximal success probability of transforming 'k' uses of a set of 'n' random isometry operations from 'd'-dimensional space to 'D'-dimensional space into its inverse map. The interpreter CVX/YALMIP is used to conduct SDP.
- [run_isometry_inversion_sod_cvx.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_inversion_sod_cvx.m) / [run_isometry_inversion_sod_yalmip.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_inversion_sod_yalmip.m): Code that obtains the maximal success probability of transforming 'k' uses of a set of 'n' random isometry operations from 'd'-dimensional space to 'D'-dimensional space into its inverse map in "success-or-draw" manner. The interpreter CVX/YALMIP is used to conduct SDP.
- [run_isometry_transposition_cvx.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_transposition_cvx.m) / [run_isometry_transposition_yalmip.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_transposition_yalmip.m): Code that obtains the maximal success probability of transforming 'k' uses of a set of 'n' random isometry operations from 'd'-dimensional space to 'D'-dimensional space into its transpose map. The interpreter CVX/YALMIP is used to conduct SDP.
- [run_isometry_complex_conjugation_cvx.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_complex_conjugation_cvx.m) / [run_isometry_complex_conjugation_yalmip.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_complex_conjugation_yalmip.m): Code that obtains the maximal success probability of transforming 'k' uses of a set of 'n' random isometry operations from 'd'-dimensional space to 'D'-dimensional space into its complex conjugate map. The interpreter CVX/YALMIP is used to conduct SDP.
- [run_isometry_pseudo_complex_conjugation_cvx.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/run_isometry_pseudo_complex_conjugation_cvx.m) / [run_isometry_pseudo_complex_conjugation_yalmip.m](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_pseudo_complex_conjugation_yalmip.m): Code that obtains the maximal success probability of transforming 'k' uses of a set of 'n' random isometry operations from 'd'-dimensional space to 'D'-dimensional space into its pseudo complex conjugate map. The interpreter CVX/YALMIP is used to conduct SDP.

These scripts make use of the following subfunctions:

- [find_dimension_isometry](https://github.com/sy3104/isometry_inversion/blob/main/find_dimension/find_dimension_isometry.m): Function that evaluate the dimention of the linear space spanned by 'k' copies of isometry operations from 'd'-dimensional space to 'D'-dimensional space
- [is_parallel_protocol_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/is_parallel_protocol_cvx.m) / [is_parallel_protocol_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/is_parallel_protocol_yalmip.m): Code that imposes SDP constrains on a variable C to enforce that C is a valid parallel superchannel
- [is_sequential_protocol_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/is_sequential_protocol_cvx.m) / [is_sequential_protocol_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/is_sequential_protocol_yalmip.m): Code that imposes SDP constrains on a variable C to enforce that C is a valid sequential superchannel
- [is_general_protocol_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/is_general_protocol_cvx.m) / [is_general_protocol_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/is_general_protocol_yalmip.m): Code that imposes SDP constrains on a variable C to enforce that C is a valid general superchannel
- [maxp_isometry_inversion_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/maxp_isometry_inversion_cvx.m) / [maxp_isometry_inversion_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_inversion_yalmip.m): Code that evaluates the maximal probability of transforming the set of k isometries encoded in Choi operators Jin and obtain the inverse maps
- [maxp_isometry_inversion_sod_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/maxp_isometry_inversion_sod_cvx.m) / [maxp_isometry_inversion_sod_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_inversion_sod_yalmip.m): Code that evaluates the maximal probability of transforming the set of k isometries encoded in Choi operators Jin and obtain the inverse maps in "success-or-draw" manner
- [maxp_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/maxp_cvx.m) / [maxp_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_yalmip.m): Code that evaluates the maximal probability of transforming the set of k isometries encoded in Choi operators Jin and obtain the isometry operators Jout
- [maxp_isometry_pseudo_complex_conjugation_cvx](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_cvx/maxp_isometry_pseudo_complex_conjugation_cvx.m) / [maxp_isometry_pseudo_complex_conjugation_yalmip](https://github.com/sy3104/isometry_inversion/blob/main/isometry_inversion_yalmip/maxp_isometry_pseudo_complex_conjugation_yalmip.m): Code that evaluates the maximal probability of transforming the set of k isometries encoded in Choi operators Jin and obtain the pseudo complex conjugate maps


## License

This code is under the [MIT license](https://opensource.org/licenses/MIT).


[1]:https://arxiv.org/abs/2110.00258
[2]:https://github.com/mtcq/unitary_inverse
[3]:https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.210502
[4]:https://arxiv.org/abs/1810.06944
[5]:http://cvxr.com
[6]:https://yalmip.github.io
[7]:https://qetlab.com
[8]:https://doi.org/10.22331/q-2023-03-20-957
