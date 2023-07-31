# Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources
### Authors: Thomas Caussade, Esteban Paduro, Matías Courdurier, Eduardo Cerpa, Warren M. Grill, Leonel E. Medina

This repository contains all data and codes to reproduce the article: "Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources"

## Download

From a command line:
```
git clone https://github.com/tcaussade/quasistatic-khz-sources
```

## Contents
The scripts are organized in two folders:

* `Database`: Contains all relevant generated data to produce the figures. See `readme_database.txt` for more details.
* `Codes`: Contains all relevant codes to reproduce the simulations. The subfolder `inh_src` contains all the functionalities employed to compute the exact solution of the inhomogeneous Helmholtz equation, and then applying the generated electric potential to the MRG model fiber. See `readme_codes.txt` for more details.

For a quick use of the improved quasi-static approximation see `Codes/improved_sqs.py`.

## References

[1] (in review process)

## License

[MIT](LICENSE)
