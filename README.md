
# Variational Bayes for NMR

A framework for approximate Bayesian inference from off-grid, non-uniformly,
non-simultaneously sampled multidimensional NMR data.

The **vbnmr** library extends the
[variational feature learning](http://github.com/geekysuavo/vfl) library
with a **quad** factor type and a **vfgp** model type. The **quad**
factor encodes multidimensional multicomplex quadrature cosines, which
are a fundamental component of the NMR signal model. The **vfgp**
model type was originally used for active learning, but this has
been superseded by the GPU-assisted search structure in VFL.

While the quadrature factor is suited for both on-grid and off-grid
learning _via_ gradient-based optimization, it also supports fast
mean-field inference when the data lie on an integer grid.

More details are in preparation for submission to:

> Worley, B., Malliavin, T. E., Nilges, M., _Active nonuniform sampling_,
> Journal of Magnetic Resonance, 2017.

## Introduction

FIXME.

### Modeling NMR data

FIXME.

### Installation

FIXME.

## Licensing

The **vbnmr** library is released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

And as always, enjoy!

*~ Brad.*

