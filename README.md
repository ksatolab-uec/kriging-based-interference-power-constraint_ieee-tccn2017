# Demonstration of Kriging-Based Interference Power Constraint
This code demonstrates the Kriging-based transmission power optimization for spectrum sharing, proposed in the following paper:
* K. Sato and T. Fujii, IEEE Trans. Cogn. Commun. Netw., vol.3, no.1, pp.13-25, March 2017.
* https://ieeexplore.ieee.org/abstract/document/7817747 (Open Access)

# Note
This code simplifies the system model and method from the original article.
(We implemented the original code with C)
The performances of this demonstration may be different in part from results in the article.
For example...
* Single SU transmitter (TCCN: multiple SUs)
* Binning-based semivariogram modeling (TCCN: residual maximum likelihood)
* ```kriging.py``` is based on https://github.com/ksatolab-uec/radiomap-construction

# How To Use
This code works with the following command:
```
$python main.py
```

# License

The MIT License (MIT)

Copyright (c) 2022 Koya SATO.
