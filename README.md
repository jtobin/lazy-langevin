Haskell implementation of some Langevin diffusion samplers.

- In general the observed information will not be positive definite; this is causing a problem in your matrix square root calculation in the RiemannianLangevin module.  Better approximation available here?
- Replacing the observed information with a different metric isn't really a general-purpose solution.
- According to Lange's NAFS, replace the observed information with a pd approximating matrix.

