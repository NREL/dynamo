**Problem:** Identifying Flexible Storage Sizing Options in the face of uncertain PV adoption.

This is a bit more complex, domain relevant example. In particular it includes inter-period coupling and hence requires more sophisticated uncertainty definitions. The default run for the example is a dynamo port of the analysis from: B. Palmintier, D. Krishnamurthy, and H. Wu, "Design Flexibility for  Uncertain Distributed Generation from Photovoltaics," in *Proceedings of the Innovative  Smart Grid Technologies Conference 2016*, Minneapolis, MN, 2016.

The original analysis for this work was coducted using a fairly complex Excel Spreadsheet with pre-computed distribution simulations using OpenDSS. Both the Excel spreadsheet and paper are included for reference.

Start with 

```[simple_storage_pv_problem , isgt2016_results] = SimpleStoragePv_demo('isgt_lookup', 'dp')```

**NOTE:** The OpenDSS link and examples of using the scalable problem definition are not yet fully documented or included