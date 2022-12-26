---
sort: 1
---

# HBT

## Subhalo files

Subhalo files are generated with HBT+ code ([Han et al. 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..604H/abstract)). You can find a general description of output at [this website](https://github.com/Kambrian/HBTplus/wiki/Outputs). Dataset is saved in such a **'Subhalo_xxx.ppp.hdf5'** file.
```
Subhalo_xxx.ppp.hdf5:
xxx is the snapshot number from 0 to 127. Here 127 corresponds to redshift z=0. 
ppp is the number of divided regions. 
```
We show the adopted parameters during HBT+ running in the following.
## Code parameters

| name          | value | description        |
| ------------- | -------------------- |  ------------------|
|`SnapshotPath`  | /home/cossim/Jiutian/M1000/simu | snapshot path|
|`HaloPath`  |/home/cossim/Jiutian/M1000/groups| host halo path|
|`SubhaloPath`  |/home/cossim/Jiutian/M1000/hbt| subhalo path|
|`SnapshotFileBase` | snapshot||
|`MaxSnapshotIndex` | 127| |
|`BoxSize`  | 1000$Mpc/h$| |
|`SofteningHalo`  | 0.004||
|`SnapshotFormat` | gadget||
|`GroupFileFormat` | gadget3_MXXL||
|`MaxConcurrentIO` | 512||
|`MinSnapshotIndex`  | 0||
|`MinNumPartOfSub` | 20 | minimum number of particles in subhalos|
|`MassInMsunh`  | 1e+10 $M_{\odot}/h$| mass unit|
|`LengthInMpch` | 1 Mpc/h| length unit |
|`VelInKmS`  | 1 km/s| velocity unit |
|`PeriodicBoundaryOn`  | 1||
|`SnapshotHasIdBlock`  | 1||
|`ParticleIdRankStyle`   | 0||
|`ParticleIdNeedHash` | 1||
|`SnapshotIdUnsigned`   | 0||
|`SaveSubParticleProperties` | 0||
|`MergeTrappedSubhalos` | 1||
|`MajorProgenitorMassRatio`  | 0.8||
|`BoundMassPrecision`  | 0.995||
|`SourceSubRelaxFactor`  | 3||
|`SubCoreSizeFactor`  | 0.25||
|`SubCoreSizeMin`  | 20||
|`TreeAllocFactor ` | 0.8||
|`TreeNodeOpenAngle`  | 0.45||
|`TreeMinNumOfCells`  | 10||
|`MaxSampleSizeOfPotentialEstimate` | 1000||
|`RefineMostboundParticle` | 1||



We show the information of **'Subhalo_xxx.ppp.hdf5'** file in the following.

## keys value

| name          | description        |
| ------------- | -------------------- | 
| `Cosmology`       |     cosmology parameters       |                   
| `NestedSubhalos` |    |                   
| `NumberOfFakeHalos`         |  |                   
| `NumberOfFiles`     |     number of files in present snapshot       |                   
| `NumberOfNewSubhalos`        |                 |                  
| `NumberOfSubhalosInAllFiles`   | number of subhalos in all files               |    
| `SnapshotId` |       snapshot ID        |  
| `SubhaloParticles` |      particles in subhalos        | 
| `Subhalos` |      details of found subhalos         |  

```note
The total number of subhalos in each file is not included!
```

## fields in `Cosmology`

| name          |  value        | description       |
| ------------- | -------------------- | ----------------- |
| `HubbleParam`  | 100 | H |
| `OmegaLambda0` | 0.6889 | $O_{\Lambda0}$ |
| `OmegaM0`    | 0.3111 | $O_{m0}$  |
| `ScaleFactor` |     1 |   |


## fields in `Subhalos`

| name          | description       | unit      |  data type  |
| ------------- | -------------------- | ----------------- | ----------|
| `TrackId` | unique subhalo ID | | int64 |
| `Nbound`   |number of bound particles | | int64 |
| `Mbound`  | total mass of bound particles | | float32 |
| `HostHaloId` | host halo ID| | int64 |
| `Rank` |  the order of subhaloes inside the group | | int64 |
|`Depth` | the level of the subhalo in the merging hierarchy | | int32 |
|`LastMaxMass`| | | float32 |
|`SnapshotIndexOfLastMaxMass`| | | int32 |
|`SnapshotIndexOfLastIsolation`| | | int32 |
|`SnapshotIndexOfBirth`| | | int32 |
|`SnapshotIndexOfDeath`| | | int32 |
|`SnapshotIndexOfSink`| | | int32 |
|`RmaxComoving`| | | float32 |
|`VmaxPhysical`| | | float32 |
|`LastMaxVmaxPhysical`| | | float32 |
|`SnapshotIndexOfLastMaxVmax`| | | int32 |
|`R2SigmaComoving`| | |float32 |
|`RHalfComoving`| | | float32|
|`BoundR200CritComoving` | R200 under critical density | Mpc/h | float32|
|`BoundM200Crit` | M200 under critica density  | $10^{10}$ $M_{\odot}/h$ | float32|
|`SpecificSelfPotentialEnergy`| | | float32|
|`SpecificSelfKineticEnergy`| | | float32|
|`SpecificAngularMomentum`| | |float32 |
|`InertialEigenVector`| | |float32 |
|`InertialEigenVectorWeighted`| | | float32|
|`InertialTensor`| | | float32|
|`InertialTensorWeighted`| | | float32|
|`ComovingAveragePosition` | average position from bound particles | Mpc/h | float32|
|`PhysicalAverageVelocity` | average velocity from bound particles  | km/s | float32|
|`ComovingMostBoundPosition` | average position from the most bound particles | Mpc/h | float32|
|`PhysicalMostBoundVelocity`| average velocity from the most bound particles | km/s | float32|
|`MostBoundParticleId`| | | int64|
|`SinkTrackId`| | | int64|

```tip
`Nbound`: Subhalos with `Nbound=1` are also listed. 

`Rank`: `Rank=0` indicating the most-massive subhalo inside each group (i.e., the main/central subhalo).

`Depth`: A central subhalo has `Depth=0`; those directly merged to the host halo of the central have `Depth=1` (i.e., sub-subhalos); those directly merged to depth=1 subhalos have `Depth=2` (i.e., sub-sub-subhalos). 
```


