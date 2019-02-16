# PyBaMM

[![Build Status](https://travis-ci.org/tinosulzer/PyBaMM.svg?branch=master)](https://travis-ci.org/tinosulzer/PyBaMM)
[![Documentation Status](https://readthedocs.org/projects/pybamm/badge/?version=latest)](https://pybamm.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/tinosulzer/PyBaMM/branch/master/graph/badge.svg)](https://codecov.io/gh/tinosulzer/PyBaMM)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tinosulzer/PyBaMM/master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Python Battery Mathematical Modelling solves continuum models for batteries, using both numerical methods and asymptotic analysis.

## How do I use PyBaMM?

The easiest way to use PyBaMM is to run a 1A constant-current discharge with a model of your choice with the default settings:
```python3
model = pybamm.lead_acid.LOQS()
sim = pybamm.Simulation(model)
sim.run()
sim.plot()
```
However, much greater customisation is available. It is possible to change the parameter values, geometry, submesh type, number of submesh points, methods for spatial discretisation and solver for time-stepping:
```python3
model = pybamm.lead_acid.LOQS()
parameter_values = pybamm.ParameterValues(
    "input/parameters/lithium-ion/parameters/LCO.csv",
    {"I_typ": 1, "current function": "constant_current.py"},
)

geometry = pybamm.Geometry1DMacro()
submesh_types = {
  "negative electrode": pybamm.Uniform1DSubMesh,
  "separator": pybamm.Uniform1DSubMesh,
  "positive electrode": pybamm.Uniform1DSubMesh,
}
submesh_pts = {
    "negative electrode": {"x": 40},
    "separator": {"x": 25},
    "positive electrode": {"x": 35},
}
spatial_methods = {
    "negative electrode": pybamm.FiniteVolume,
    "separator": pybamm.FiniteVolume,
    "positive electrode": pybamm.FiniteVolume,
}
solver = pybamm.ScipySolver(method="RK45")
sim = pybamm.Simulation(
    model,
    parameter_values=parameter_values,
    geometry=geometry,
    submesh_types=submesh_types,
    submesh_pts=submesh_pts,
    spatial_methods=spatial_methods,
    solver=solver,
)
sim.run()
sim.plot()
```
These scripts are available as [Jupyter notebooks](examples/getting-started.ipynb), where you can also try out different settings.
Further details on how to customise simulations and add your own models and settings are provided in the [examples and documentation](#examples-and-documentation).

### Examples and documentation

PyBaMM comes with a number of [detailed examples](examples/README.md), hosted here on github. In addition, there is a [full API documentation](http://pybamm.readthedocs.io/), hosted on [Read The Docs](readthedocs.io).

## How can I install PyBaMM?

You'll need the following requirements:

- Python 2.7 or Python 3.4+
- Python libraries: `numpy` `scipy` `pandas` `matplotlib`

These can easily be installed using `pip`. To do this, first make sure you have the latest version of pip installed:

```
$ pip install --upgrade pip
```

Then navigate to the path where you downloaded PyBaMM to, and install both PyBaMM and its dependencies by typing:

```
$ pip install .
```

Or, if you want to install PyBaMM as a [developer](CONTRIBUTING.md), use

```
$ pip install -e .[dev,docs]
```

To uninstall again, type

```
$ pip uninstall pybamm
```

## How can I contribute to PyBaMM?

If you'd like to help us develop PyBaMM by adding new methods, writing documentation, or fixing embarrassing bugs, please have a look at these [guidelines](CONTRIBUTING.md) first.

## Licensing

PyBaMM is fully open source. For more information about its license, see [LICENSE](./LICENSE.txt).
