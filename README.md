# lazy-langevin [![Build Status](https://secure.travis-ci.org/jtobin/lazy-langevin.png)](http://travis-ci.org/jtobin/lazy-langevin)

Sampling via gradient-based diffusion.

See the *Examples* folder for example usage.

## TODO

- In general the observed information will not be positive definite; this is causing a problem in the matrix square root calculation in the RiemannianLangevin module.  Replace the observed information with a pd approximating matrix.
- Note that the inclusion of hmatrix requires this to be GPL, so that needs to be changed.

