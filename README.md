# PySCES Polariton Add-on

This package serves as an add-on for [PySCES (Python Suite for Computational Electronic Structure)](https://github.com/AnanthGroup/PySCES), extending its functionality to handle polariton simulations.

## Overview

PySCES Polariton Add-on integrates with the main PySCES package to provide specialized tools for modeling and simulating polaritons - quasiparticles that result from strong coupling between photons and excited states of matter.

## Requirements

- Python 3.9+
- NumPy
- SciPy
- PySCES

## Installation

1. First, install the core PySCES package:
   ```bash
   git clone https://github.com/AnanthGroup/PySCES.git
   cd PySCES
   pip install -e .
   ```

2. Then, install this polariton add-on:
   ```bash
   git clone https://github.com/yourusername/pysces_polariton_add_on.git
   cd pysces_polariton_add_on
   pip install -e .
   ```

## Usage

```python
import pysces as ps
# Import the polariton add-on
import pysces_polariton as psp

# Example code for polariton simulation
# ...
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the same license as PySCES - please see the main PySCES repository for details.

```