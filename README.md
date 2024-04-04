E3SM chemistry diagnostic package (ChemDyg) is designed for E3SM chemistry development. There are total 12 types of plot and 3 types of table for model-to-model and model-to-observation comparison.

Index page and example figures:
https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.lee1061/20220914.PAN.MZThet.v2.LR.bi-grid.amip.chemUCI_Linozv3/e3sm_chem_diags/plots/

ChemDyg is executed by zppy, which is a post-processing toolchain for E3SM written in Python (https://github.com/E3SM-Project/zppy).

*How to run ChemDyg via zppy script on chrysalis?
```
source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
wget https://raw.githubusercontent.com/E3SM-Project/ChemDyg/main/ChemDyg_example_script.cfg
zppy -c ChemDyg_example_script.cfg
```
**Note that the users need to change the'e3sm_unified' path on different zppy-supported machines (https://e3sm-project.github.io/zppy/_build/html/main/getting_started.html#activate-e3sm-unified-environment). 

For people who cannot access zppy-supported machines, the observational data used in ChemDyg are available at https://doi.org/10.5281/zenodo.10905766 and https://web.lcrc.anl.gov/public/e3sm/diagnostics/observations/Atm/ChemDyg_inputs/.  

## License

Copyright (c) 2023, Energy Exascale Earth System Model Project
All rights reserved

SPDX-License-Identifier: (BSD-3-Clause)

See [LICENSE](./LICENSE) for details

Unlimited Open Source - BSD 3-clause Distribution
'LLNL-CODE-819717`
