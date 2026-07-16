# Backbone Rigid Invariants

## Introduction

This program is enables:

1. Cleaning the PDB dataset (from [RCSB](https://www.rcsb.org/)) using the regulations we defined to obtain the cleaned protein chains.
2. Computing Backbone Rigid Invariants (BRIs) with triangular invariants (trins) and Average Invariant (Brain) for a given cleaned chain.
3. Quick comparison of BRIs, generating distances between chains to find near-duplicates of proteins.

## Getting started

### Prerequisites

Python versions >=3.10.

### Installation

Use pip to install backbone-rigid-invariant:

> Currently, you may have to install with released distributions.

Install with a `wheel` file:
```shell
pip install backbone_rigid_invariant-1.1.0-py3-none-any.whl
```

Or install with a `.tar.gz` compression file:
```shell
pip install backbone_rigid_invariant-1.1.0.tar.gz
```

Then import the package by ```import bri```.

### Read chains & compute invariants

This example reads one single protein chain from PDB and computes the Backbone Rigid Invariant:

```python
from bri import Chain


# create a Chain instance using pdb_id, entity_id, model_id, chain_id, start of residue_id, and chain length
reference_1hho = Chain('1hho', 1, 1, 'A', 1, 141)

# get Backbone Rigid Invariant
invariant = reference_1hho.invariant
```

The variable `invariant` is a `Pandas.DataFrame` object with *n* rows, where *n* represents the length of the chain.

### Implement perturbation & validate Lipschitz Constant

To perturb a protein chain and validate Lipschitz Constants:

```python
from bri import Chain, get_e_lipschitz_constant, get_theo_lipschitz_constant


# create chain instances
reference_1hho = Chain('1hho', 1, 1, 'A', 1, 141)	
perturb_1hho = Chain('1hho', 1, 1, 'A', 1, 141)

# set perturbation
perturb_1hho.perturb_radius = 0.01

# compute lipschitz_constant
experimental_constant = get_e_lipschitz_constant(reference_1hho, perturb_1hho)
theoretical_constant = get_theo_lipschitz_constant(reference_1hho, perturb_1hho)
```

Specifying the attribute `Chain.perturb_radius = epsilon` induces a perturbation within a sphere with `epsilon` radius, and updates `Chain.invariant` at the same time.

`get_e_lipschitz_constant()` and `get_theo_lipschitz_constant()` use the `Chain` objects they get to compute experimental and theoretical Lipschitz Constant namely.

### Plot Backbone Invariant Diagram (BID) & Backbone Invariant Barcode (BIB)

BID and BIB are plots that can straight show some features of a chain. Class `Chain` integrates methods to generate these plots powered by library `matplotlib`:

```python
from bri import Chain

# initial chain instances
reference_1hho = Chain('1hho', 1, 1, 'A', 1, 141)

# show and save Backbone Invariants Diagram
reference_1hho.generate_BID().savefig('1HHO_BID.png')
reference_1hho.generate_BIB().savefig('1HHO_BIB.png')
```

`.generate_BID()` and `.generate_BIB()` show generated plots and return `matplotlib.pyplot.1` objects. So it is easy to save with its `.savefig()` method, with arguments it supports.

<!-- **For use details， please see the [example notebook](examples/BRI_examples.ipynb)** -->

