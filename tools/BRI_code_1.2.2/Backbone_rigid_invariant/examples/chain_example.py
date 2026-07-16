# -*- coding = utf-8 -*-
# @Time: 2024/7/9 14:14
# @File: chain_example.PY
from util.pdbx2df import Chain
from util.invariant import get_e_lipschitz_constant, get_theo_lipschitz_constant


# initial chain instances
reference_1hho = Chain('1hho', 1, 1, 'A', 1, 141)
perturb_1hho = Chain('1hho', 1, 1, 'A', 1, 141)

# set perturbation
perturb_1hho.set_perturbation(0.01)

# show perturbed coordinates and invariants
print(perturb_1hho.coordinates)
print(perturb_1hho.invariant)

# check lipschitz_constant
print(get_e_lipschitz_constant(reference_1hho, perturb_1hho))
print(get_theo_lipschitz_constant(reference_1hho, perturb_1hho))

# show and save Backbone Invariants Diagram
reference_1hho.generate_BID().savefig('1HHO_BID.png')

