
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

In short, **vbnmr** extends [VFL](http://github.com/geekysuavo/vfl) to
support parameteric models of the kind introduced by G. L. Bretthorst,
R. Chylla, and J. L. Markley. Thus, we directly obtain an engine for
parameter inference (_via_ fixed-form variational Bayes), in the
NMR signal model.

### Modeling NMR data

In the simplest possible case, we may treat the NMR signal as a sum of
[decaying](https://en.wikipedia.org/wiki/Exponential_decay)
[multicomplex](http://en.wikipedia.org/wiki/Multicomplex_number)
[quadrature](https://en.wikipedia.org/wiki/Quadrature_phase)
[sinusoids](https://en.wikipedia.org/wiki/Sine_wave). In most
cases, we have a good prior estimates of:

 * Number of signals (specified as a count, _M_)
 * Measurement noise (specified a precision, _tau_)
 * Signal-to-noise ratio (specified as a power ratio, _nu_)

In addition, we should have more-or-less decent prior estimates of:

 * Mean and variance of decay rates (specified via _A_ and _B_)
 * Mean and variance of signal frequencies (specified via _U_ and _T_)

This state of knowledge may be represented in **vflang** as follows:

```
dat = data();
# ... fill the dataset with measurements ...
mdl = tauvfr(
  tau: tau, nu: nu, data: dat,
  factors: M * [decay(alpha: A, beta: B, fixed: true),
                quad(mu: U, tau: T, ftsize: 65536)]
);
```

The keen observer will notice that we have fixed the decay factor and
specified the **ftsize** property of each quadrature factor. Doing so
enables the use of fast mean-field inference when the data lie on an
integer grid:

```
opt = mf(model: mdl, maxIters: 1);
opt.execute();
```

If required, we can then perform full inference using natural gradient
ascent over all the factor parameters:

```
for j in std.range(n: mdl.M) {
  mdl.factors[j][0].fixed = false; # un-fix each decay factor.
}
opt = fg(model: mdl); # run natural gradient optimization.
opt.execute();
```

The above **vflang** is a small taste of the possible methods of inference
using the **vbnmr** extension to **vfl**. More complete documentation will
be made available in the near future.

## Installation

The **vbnmr** library is written in C99-compliant source code (with GNU
extensions). The [FFTW](http://www.fftw.org) library is required for
compilation.

You can compile and install **vbnmr** as follows:

```bash
git clone git://github.com/geekysuavo/vbnmr.git
cd vbnmr
make
sudo make install
```

Following installation, **vbnmr** can be used by inclusion as a shared
library in C programs, or by importing it as a module into **vflang**
scripts.

## Licensing

The **vbnmr** library is released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

And as always, enjoy!

*~ Brad.*

