---
sort: 2
---

# Subfind

## Groups files

Groups files are generated with subfind code. You can find a general description of output at [this website](). Dataset is saved in such a **'subhalo_tab_xxx.ppp'** file.
```
subhalo_tab_xxx.ppp:
xxx is the snapshot number from 0 to 127. Here 127 corresponds to redshift z=0. 
ppp is the number of divided regions. 
```
We show the blocks of **'subhalo_tab_xxx.ppp'** file in the following.

## blocks for file information

| name          | description        |
| ------------- | -------------------- | 
|`ngroups` | | 
|`nids` | |
|`nfiles` | | 
|`nsubs` | | 
|`totngrp`| |
|`totnsub`| |



## blocks for hosthalos

| name          | description        |
| ------------- | -------------------- | 
|`group_len` | | 
|`group_offset` | |
|`group_nr` | host halo ID | 
|`group_cm` | | 
|`group_vel `| | 
|`group_pos` | | 
|`group_m_mean200` | | 
|`group_m_crit200`| | 
|`group_m_tophat200` | | 
|`group_veldisp` | |
|`group_veldisp_mean200` | | 
|`group_veldisp_crit200` | | 
|`group_veldisp_tophat200` | | 
|`group_nsubs `| | 
|`group_firstsub` | | 

## blocks for subhalos

| name          | description        |
| ------------- | -------------------- | 
|`sub_len ` | | 
|`sub_offset` | |
|`sub_grnr` | host halo ID| 
|`sub_nr` | subhalo ID | 
|`sub_pos `| | 
|`sub_vel` | | 
|`sub_cm`| | 
|`sub_spin`| | 
|`sub_veldisp` | | 
|`sub_vmax` | |
|`sub_vmaxrad` | | 
|`sub_halfmassrad`| | 
|`sub_ebind` | | 
|`sub_pot` | | 
|`sub_idbm`| | 

## supplement
| name          | description        |   unit  |
| ------------- | -------------------- |  ------ |
|`group_mass`   | total host halo mass in $log_{10}$ | $M_{\odot}/h$ |
|`sub_mass`   | total subhalo mass in $log_{10}$ | $M_{\odot}/h$ |


