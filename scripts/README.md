Validation Replication Package
==============================

This folder contains the scripts and data used for independent validation
of Bruehl and Villarroel (2025) findings.

Author: Brian Doherty
Contact: briandohertyresearch@gmail.com
Date: January 2026


Contents
--------

Scripts:
  nuclear_transient_correlation.py - Chi-square, negative binomial, permutation tests
  earth_shadow_validation.py - Earth shadow geometry and transient classification

Data (in data/ subfolder):
  Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx - Original dataset from Dr. Bruehl
  SUPERVIKTIG_HELAVASCO.csv - VASCO transient catalog
  SUPERVIKTIG_HELAVASCO_validated_v4.csv - Validated transient catalog


Requirements
------------

Python 3.10 or higher

Dependencies:
  pandas
  numpy
  scipy
  statsmodels
  astropy
  openpyxl (for Excel file reading)


Usage
-----

From this directory:

  python nuclear_transient_correlation.py
  python earth_shadow_validation.py

Results are saved to the results/ subfolder within this directory.


Methodology
-----------

See VALIDATION_METHODOLOGY.md for complete methodology documentation.


Data Sources
------------

Nuclear test dates: DOE/NV-209 Rev 16, Johnston's Archive
Transient data: Provided by Dr. Stephen Bruehl
VASCO catalog: VASCO project (Villarroel et al.)

