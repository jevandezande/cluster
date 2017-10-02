# Cluster Calculation Code

Generates clusters for X-ray computations following the model of Maganas et al.
(DOI: 10.1039/c3cp50709b).


## Description

Generates clusters that are separated into three regions: qc (where quantum
chemical methods are applied), br (which includes capped ECPs and point
charges), and pc (point charges). Currently only supports ORCA (PRs with
support for other codes would be greatly appreciated).


## Requirements

* sq - submission to cluster
* cclib - reading of output files
