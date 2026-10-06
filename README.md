# Dirichlet-L-functions
Fortran code for computing some Dirichlet L-functions (and their derivatives) in the complex plane, as well as quadrature coefficients for computing these Dirichlet L-functions in high-precision. 

# Dirichlet-L-functions

Supplementary materials and Fortran code for the paper **“Simple and accurate approximations to Dirichlet L-functions”** by A. Kuznetsov (2026); see [Reference 1](#references).

This repository provides:

- [Fortran modules](Fortran-modules) for computing Dirichlet $L$-functions $L(s,\chi)$ and their derivatives $L'(s,\chi)$ for complex $s$, using quadruple precision.
- [Quadrature coefficients](omega_lambda_coefficients) for the approximations $L_p(s,\chi)$ described in the paper, supplied at several precisions.

## Supported characters

The code and coefficient tables cover all 16 primitive Dirichlet characters with modulus $3\leq q\leq10$. Character names $\chi_{q,a}$ use the Conrey labels employed by the [LMFDB](https://www.lmfdb.org/).

| Modulus $q$ | Character labels $a$ |
|---|---|
| 3 | 2 |
| 4 | 3 |
| 5 | 2, 3, 4 |
| 7 | 2, 3, 4, 5, 6 |
| 8 | 3, 5 |
| 9 | 2, 4, 5, 7 |

For example, $\chi_{3,2}$ has values $(1,-1,0)$ on the residue classes $1,2,3$, and $\chi_{4,3}$ has values $(1,0,-1,0)$ on $1,2,3,4$. The latter gives the Dirichlet beta function.

## Quadrature coefficients

The directory [omega_lambda_coefficients](omega_lambda_coefficients) contains one ZIP archive for each character. For example, `chi_3_2.zip` contains the coefficients for $\chi_{3,2}$.

Each archive contains subdirectories for the following values of $p$:

| Number of quadrature nodes $p$ | Significant decimal digits stored per real or imaginary component |
|---|---|
| 20 | 20 |
| 30 | 30 |
| 40 | 40 |
| 50 | 50 |
| 70 | 80 |
| 100 | 100 |

Here $p$ is the total number of quadrature nodes. The number of stored digits describes the data format; the accuracy of the resulting approximation also depends on $s$, $\chi$, and $p$.

### Filenames and indexing

For each $p$ and each residue class $r=0,\ldots,q-1$, the archive contains two files:

```text
lambda_chi_{q,a}_p=<p>_r=<r>.txt
omega_chi_{q,a}_p=<p>_r=<r>.txt
```

For example, after extracting `chi_3_2.zip`, the files for $p=20$ and $r=0$ are:

```text
chi_3_2/p=20/lambda_chi_{3,2}_p=20_r=0.txt
chi_3_2/p=20/omega_chi_{3,2}_p=20_r=0.txt
```

The lambda file contains the nodes $\lambda_{\chi,p,j}^{(r)}$, and the omega file contains the weights $\omega_{\chi,p,j}^{(r)}$, in increasing index order $j=1,\ldots,p$. **Each file contains $p$ complex coefficients, stored as $2p$ real numbers.** In each consecutive pair, the first number is the real part and the second is the imaginary part. Nodes and weights with the same index $j$ belong together.


### Example

The first two numbers in `lambda_chi_{3,2}_p=20_r=0.txt` are

```text
2.6331058366723794812e0
-2.4826329891642160635e0
```

They represent the complex coefficient

$$
\lambda_{\chi_{3,2},20,1}^{(0)}
=2.6331058366723794812
-i 2.4826329891642160635.
$$

Similarly, the first two numbers in `omega_chi_{3,2}_p=20_r=0.txt` are

```text
-1.8668713937563482763e-8
-1.7367471766222593860e-8
```

They represent

$$
\omega_{\chi_{3,2},20,1}^{(0)}
=-1.8668713937563482763\times10^{-8}
-i 1.7367471766222593860\times10^{-8}.
$$

Here $i$ denotes the imaginary unit.

### Continuation lines

The files for $p=70$ and $p=100$ use the MPFUN2020 output format, in which a backslash `\` at the end of a line continues the same number on the next line.

For example, the first number in `lambda_chi_{3,2}_p=70_r=0.txt` is written as

```text
5.401150991772538617196592479085455275359490671052862625718202331524952280740539\
9e0
```

and should be interpreted as

```text
5.4011509917725386171965924790854552753594906710528626257182023315249522807405399e0
```

Remove each continuation backslash and the following newline, joining the adjacent digit strings without inserting a space, before splitting the file into numbers and grouping consecutive values into real/imaginary pairs. Use an appropriate precision when reading the coefficients to retain the stored digits. See the [MPFUN2020 documentation](https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf) for details.

## Fortran modules

The directory [Fortran-modules](Fortran-modules) contains the following source files:

| Module file | Functions computing $L(s,\chi)$ |
|---|---|
| [L_chi_3_2_module.f90](Fortran-modules/L_chi_3_2_module.f90) | `L_3_2` |
| [L_chi_4_3_module.f90](Fortran-modules/L_chi_4_3_module.f90) | `L_4_3` |
| [L_chi_5_2_and_chi_5_3_module.f90](Fortran-modules/L_chi_5_2_and_chi_5_3_module.f90) | `L_5_2`, `L_5_3` |
| [L_chi_5_4_module.f90](Fortran-modules/L_chi_5_4_module.f90) | `L_5_4` |
| [L_chi_7_2_and_chi_7_4_module.f90](Fortran-modules/L_chi_7_2_and_chi_7_4_module.f90) | `L_7_2`, `L_7_4` |
| [L_chi_7_3_and_chi_7_5_module.f90](Fortran-modules/L_chi_7_3_and_chi_7_5_module.f90) | `L_7_3`, `L_7_5` |
| [L_chi_7_6_module.f90](Fortran-modules/L_chi_7_6_module.f90) | `L_7_6` |
| [L_chi_8_3_module.f90](Fortran-modules/L_chi_8_3_module.f90) | `L_8_3` |
| [L_chi_8_5_module.f90](Fortran-modules/L_chi_8_5_module.f90) | `L_8_5` |
| [L_chi_9_2_and_chi_9_5_module.f90](Fortran-modules/L_chi_9_2_and_chi_9_5_module.f90) | `L_9_2`, `L_9_5` |
| [L_chi_9_4_and_chi_9_7_module.f90](Fortran-modules/L_chi_9_4_and_chi_9_7_module.f90) | `L_9_4`, `L_9_7` |

Each function has a corresponding function with `_prime` appended to its name. For example:

- `L_3_2(s)` returns the scalar value $L(s,\chi_{3,2})$.
- `L_3_2_prime(s)` returns an array of length two: element 1 is $L(s,\chi_{3,2})$, and element 2 is $L'(s,\chi_{3,2})$.

All inputs are scalars of type `complex(kind=qp)`, with

```fortran
integer, parameter :: qp = selected_real_kind(33, 4931)
```

The functions are not vectorized. The parameter `qp` is private to each module, so the calling program should declare the same kind parameter.

Each module is self-contained and includes its required coefficients and helper routines. A Fortran compiler supporting the requested real kind is sufficient; MPFUN2020 and the coefficient ZIP archives are not required to compile or run these modules.

### Example program

Save the following as `example.f90` in the repository root:

```fortran
program example
    use L_chi_3_2_module, only: L_3_2, L_3_2_prime
    implicit none

    integer, parameter :: qp = selected_real_kind(33, 4931)
    complex(kind=qp) :: s, value, values(2)

    s = cmplx(0.5_qp, 100.0_qp, kind=qp)

    value = L_3_2(s)
    print *, 'L(s) = ', value

    values = L_3_2_prime(s)
    print *, 'L(s) = ', values(1)
    print *, "L'(s) = ", values(2)
end program example
```

With GNU Fortran, compile and run from the repository root:

```bash
gfortran -O3 -ffree-line-length-none -fno-fast-math -ffp-contract=off \
    Fortran-modules/L_chi_3_2_module.f90 example.f90 -o example
./example
```

The option `-ffree-line-length-none` accommodates source lines longer than 132 columns.

## Method and accuracy

The modules combine direct summation, Euler–Maclaurin summation, and the quadrature approximation $L_{40}(s,\chi)$. They select the method according to estimated computational cost, with the switch between the Euler–Maclaurin and quadrature regions at $|\operatorname{Im}(s)|=400$. The functional equation and complex conjugation extend the computation to the other regions of the complex plane; see Section 3 of the paper.

For fixed $\chi$ and $p$, evaluating the quadrature approximation in a fixed vertical strip requires $O(t^{1/2})$ numerical operations as $t\to+\infty$, excluding coefficient precomputation.


## References

1. A. Kuznetsov. *Simple and accurate approximations to Dirichlet L-functions*. Manuscript, October 6, 2026.
2. D. H. Bailey. *MPFUN2020: A thread-safe arbitrary precision package with special functions*. [Software](https://www.davidhbailey.com/dhbsoftware/) and [documentation](https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf).

## Author and license

Alexey Kuznetsov, York University, Toronto, Canada.  
[Website](https://kuznetsovmath.ca/) · [Email](mailto:akuznets@yorku.ca)

Released under the [BSD 3-Clause License](LICENSE).
