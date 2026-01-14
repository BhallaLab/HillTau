![alt text](./Images/HillTau_Logo4_360px.png?raw=true "HillTau logo")

# HillTau
Sparse, efficient, non-ODE approximation for biochemical signaling

Copyright (C) Upinder S. Bhalla and NCBS-TIFR, 2020


## Contents of this repository:

README.md: this file

LICENSE: License terms. GPL 3.0.

PythonCode: The code to run HillTau.

CppCode: The code to run HillTau using python bindings via Pybind11 to C++.

MASH: Model Abstraction from SBML to HillTau. Utility program for model 
reduction, see 

[Documentation for MASH](Mashdoc.md)

ht2sbml: Utility program to convert HillTau (JSON format) model to a near-
	approximation in SBML, suitable for running on several ODE simulators
	like COPASI.

htgraph: Utility program to generate png or svg graphs to display reaction 
structure of HillTau model.

### HillTau

[Documentation for HillTau](Documentation.md)


The following paper discusses design and examples in detail:

HillTau: A fast, compact abstraction for model reduction in biochemical
signaling networks.<br>
**Bhalla US.**<br>
*PLoS Comput Biol.* 2021

[DOI](https://doi.org/10.1371/journal.pcbi.1009621)


## Installation

**Simple and slow Python version:**
Copy the two files hillTau.py and hillTauSchema.json from PythonCode to your
target directory.<br>
If you want to run the utilities, also copy *mash.py*, *ht2sbml.py* and 
*htgraph.py* to the target directory.

**Complex and fast C++ version:**
pip install hillTau<br>
This will help install the _much_ faster C++ version of HillTau, with the same
friendly Python interface.

*Limitations:* 
1. 	The pip install version has only been configured for Linux 
	systems. Mac will come soon, and Windows is a work in progress.
2. 	The *htgraph.py* script requires a separate installation of GraphViz.
	This may be done on Linux systems using *sudo apt-get install GraphViz*
3. 	The *mash.py* script requires a separate installation of the MOOSE
	development branch.
	This may be done on Linux systems by accessing the GitHub website
	for MOOSE, checking it out, and compiling as per the instructions in 
	INSTALL.md

	```git clone https://github.com/BhallaLab/moose-core.git```

	

Once the pip install is done, you can use *import hillTau* in any python script
where you need it. <br>
You can also run the standalone hillTau code from the command line like this:

```
hillTau model_file <arguments>
```

*Tested on:*:
-	Ubuntu 20.x
-	CentOS el7
-	More to come soon.

## Versions
The Python version of HillTau has been tested with Python 2.7.17 and Python 3.6.9<br>
The C++ version of HillTau uses Python 3.6 or higher. It is the recommended
version.

## Examples
Examples: Directory with examples

Examples/HT_MODELS: Examples of HillTau models

Examples/KKIT_MODELS: Examples of KKIT models which map to the HillTau ones.
	KKIT models are an old ODE biochemical model format compatible with
	GENESIS and MOOSE.

Examples/SBML_MODELS: Examples of SBML models which map to the HillTau ones.
	SBML is the Systems Biology Markup Language and is a standard for defining
	many kinds of systems biology models, including the current mass-action
	ones.

Examples/PaperFigures: Using the HillTau form to generate the figures for the
	reference paper for HillTau. Most of these require MOOSE to be 
	installed, but fig1.py just requires HillTau.

Other projects and papers that relate to HillTau: [Resources.md](Resources.md)
