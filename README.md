# Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources

### Authors: T. Caussade, E. Paduro, M. Courdurier, E. Cerpa, W.M. Grill, L.E. Medina

This repository contains all data and codes to reproduce the article: *Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources*

## Download

From a command line:
```
git clone https://github.com/tcaussade/quasistatic-khz-sources
```

To compile the MRG model fiber, run
```
nrnivmodl nrnmechanism
```

## Contents
The scripts are organized in two main folders:

* `Database`: Contains all relevant generated data to produce the figures. See `readme_database.txt` for more details.
* `Codes`: Contains all relevant codes to reproduce the simulations. See `readme_codes.txt` for more details.

Other secondary folders are necessary to compile the experiments correctly.
* `inh_src`: Contains all the functionalities employed to compute the exact solution of the inhomogeneous Helmholtz equation, and then applying the generated electric potential to the MRG model fiber
* `nrnfiberlib`, `nrnlib`: Contains .hoc files (MRG model)
* `nrnmechanism` : Contains .mod files.

For a quick use of the improved quasi-static approximation see `Corrected Conductivity.ipynb`.


## References

[1] Caussade T., Paduro E., Courdurier M., Cerpa E., Grill W. M., Medina L. E. (2023). "Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources". Journal of Neural Engineering (submission/review pending)

## License

[MIT](LICENSE)
