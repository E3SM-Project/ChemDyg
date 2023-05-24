2023/03/01
Hsiang-He Lee (lee1061@llnl.gov)
2023/05/24 updated

E3SM chemistry diagnostic package version 0.1.0 (ChemDyg v0.1.0) is designed for E3SM chemistry development. There are total 11 types of plot and 3 types of table for model-to-model and model-to-observation comparison.

Index page and example figures:
https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.lee1061/20220914.PAN.MZThet.v2.LR.bi-grid.amip.chemUCI_Linozv3/e3sm_chem_diags/plots/

ChemDyg is executed by zppy, which is a post-processing toolchain for E3SM written in Python (https://github.com/E3SM-Project/zppy).

# How to run ChemDyg via zppy script on chrysalis?

*Clone repository

```git clone git@github.com:E3SM-Project/zppy.git```

```cd zppy```

*Checkout branch
```git checkout golaz/plugins```

*Create environment

```module load anaconda3/2020.11```

```conda env create --name zppy_dev_chem -f conda/dev.yml```

*Activate

```source /gpfs/fs1/soft/chrysalis/manual/anaconda3/2020.11/etc/profile.d/conda.sh```

```conda activate zppy_dev_chem```

*Install chemdyg package

```install -y -c e3sm/label/chemdyg_dev chemdyg```

*Install zppy in environment

```python -m pip install .```

*Once this is done, you call is as

```zppy -c ChemDyg_example_script.cfg```


** Additional python script "DEBUG_chem_diags_timestep.py" is designed for debug purpose. The user gives the path of E3SM .h0 and .h1 outputs and other necessary information. The script will generate chem tendency table for each timestep.

See [documentation](https://e3sm-project.github.io/zppy) for more details.

## License

Copyright (c) 2023, Energy Exascale Earth System Model Project
All rights reserved

SPDX-License-Identifier: (BSD-3-Clause)

See [LICENSE](./LICENSE) for details

Unlimited Open Source - BSD 3-clause Distribution
'LLNL-CODE-819717`
