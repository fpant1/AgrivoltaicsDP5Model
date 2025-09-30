# Simulation and Design Optimisation Model for Agrivoltaics 🌱☀️

This project develops a **digital twin of an agrivoltaic scenario**, simulating both **solar energy generation** and **agricultural crop yields** under shared land use. The model serves as a **simulation and design optimisation tool**, enabling exploration of design parameters to improve both photovoltaic (PV) and crop performance.  

The overarching goal is to demonstrate how agrivoltaic systems can be optimised to maximise land productivity, providing a balance between **renewable energy production** and **sustainable agriculture**.  

---

## 🔧 Model Setup

The model recreates a modular agrivoltaic design in MATLAB, defining a range of design parameters such as:

- **Module geometry**: panel size, tilt, and tracking configuration  
- **Array structure**: row spacing, ground coverage ratio (GCR), height  
- **Environmental inputs**: solar irradiance, weather data, soil properties  
- **Crop parameters**: light use efficiency, phenology, water balance  

The system is built as a **modular framework**, allowing independent validation and analysis of its core components:  

- **Solar energy modelling**  
- **Shading and irradiance reduction**  
- **Microclimate representation**  
- **Crop growth modelling**

### Visualisation of the Setup

<div align="center">
  <img src="figures/matlab_structure.png" width="45%" alt="MATLAB Model of Agrivoltaic Structure"/>
  <img src="figures/cad_design.png" width="45%" alt="CAD Design of Agrivoltaic Structure"/>
</div>

*Left: Simplified MATLAB recreation of the agrivoltaic module. Right: CAD design of the physical module.*  

---

## ☀️ Solar Modelling

The solar model incorporates **tracking algorithms** and computes **direct and diffuse irradiance** incident on both the panels and the ground beneath the structure.  

- **Shading Factor**: quantifies the reduction in irradiance at the ground caused by the array.  
- **Sensitivity Analysis**: demonstrates the relationship between GCR and shading, validating against reported values in literature.  

<div align="center">
  <img src="figures/shading_animation.gif" width="70%" alt="GIF of shading simulation over a day"/>
</div>

*Figure: Example shading simulation visualised as irradiance distribution beneath the array.*  

---

## 🌡️ Microclimate Modelling

The microclimate model captures the effects of shading on **temperature, humidity, and soil water balance** beneath the array.  

- **Evapotranspiration**: calculated dynamically to determine plant water use.  
- **Water Stress**: integrated into the crop growth model via rainfed vs irrigated scenarios.  

This allows analysis of how **reduced solar input** interacts with **water availability**, directly influencing crop performance.

---

## 🌾 Crop Growth Modelling

The crop model uses crop-specific physiological inputs (e.g. **tomato** and **wheat**) to simulate yield under agrivoltaic conditions.  

- **Inputs**: sourced from experimental studies and growth guides.  
- **Outputs**: yield reduction under shade, validated against reported crop trials.  
- **Scenarios**: open-field baseline, agrivoltaic with shading, and irrigated vs rainfed systems.  

<div align="center">
  <img src="figures/crop_yield_comparison.png" width="70%" alt="Crop yield comparison"/>
</div>

*Figure: Comparison of modelled vs expected yields for tomato (Spain) and wheat (UK).*  

---

## 📊 Key Findings

- **Shading Validation**:  
  Average reductions of ~30% in irradiance reported in literature align with model outputs (36–59% depending on GCR).  

- **Crop Validation**:  
  - Tomatoes (Spain): 28% reduction under agrivoltaics.  
  - Wheat (UK): 33% reduction under agrivoltaics.  
  - Irrigation boosts yields by 57% (Spain) and 24% (UK).  

- **Design Insight**:  
  Tomatoes show slightly greater tolerance to shading than wheat, suggesting crop choice significantly influences system optimisation.  

---

## 📈 Optimisation Perspective

By coupling **solar yield** and **crop growth** in one framework, the tool can be used to:  

- Explore the trade-off between energy and food production.  
- Optimise design parameters such as **row spacing, height, and tracking**.  
- Support **site-specific agrivoltaic deployment strategies**.  

---

## 📚 References

- Weselek, A., et al. (2021). *Agrivoltaic systems: applications, challenges, and opportunities.*  
- AHDB (2021). *Wheat Growth Guide.*  
- Dupraz, C., et al. (2011). *Combining solar energy capture with food production: the concept of agrivoltaics.*  

---
