# dynamo

dynamo - Dynamic programming for Adaptive Modeling and Optimization

### Install

* Install MATLAB (R2017a or latter preferred)
* Clone this repository
* Open the Home>Set Path dialog and click on `Add Folder` to add the following folders to the `PATH`:
  * `$DYNAMO_Root/src`
  * `$DYNAMO_Root/extern` (Add all subfolders for this one)

### Getting Started

Explore the `example` directory. The `multi_inventory` and `Storage_Size_for_PV` examples are best well tested.

### Tests

1. Add: `$DYNAMO_Root/example/multi_inventory` and `$DYNAMO_Root/example/Storage_Size_for_PV` to your MATLAB path
2. Change to the `$DYNAMO_Root/test` folder
3. Run the following command to run all unit and doctest-based tests

    dynamo_test_all
