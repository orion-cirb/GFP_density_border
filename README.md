# VeCell

* **Developed for:** Nadia
* **Team:** Espeli
* **Date:** August 2026
* **Software:** Fiji


### Images description

3D images taken with a x60 objective

3 channels:
  1. *DAPI:* DNA
  2. *TL phase:* bacteria
  3. *GFP:* EGPF

### Plugin description

* Detect bacteria on the average intensity Z-projection of channel 2 with Omnipose
* Detect DNA on the average intensity Z-projection of channel 1 with Omnipose

* In each bacterium, return distances between bacterium centroid and DNA centroid
* Create X fragments with a thickness of Y at the border of detected bacteria
* Compute DAPI and GFP Integrated Density and ratio of fragments and interior


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Omnipose** conda environment + *bact_phase_omnitorch_0* and *bact_fluor_omnitorch_0* models


