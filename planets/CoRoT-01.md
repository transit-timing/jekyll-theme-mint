---
layout: post
style: planet
---
<script src="../js/planets.js"></script> 

## CoRoT-01

<!-- Tab links -->
<div class="tab">
  <button class="tablinks" onclick="openCity(event, 'Ephemeris')">Ephemeris</button>
  <button class="tablinks" onclick="openCity(event, 'Data')">Data</button>
  <button class="tablinks" onclick="openCity(event, 'Figures')">Figures</button>
</div>

<!-- Tab content -->
<div id="Ephemeris" class="tabcontent" markdown="1">
  <br/><br/>
  $$P=1.091419108 \pm 1.508968772 \cdot 8.3e-08 day$$
  $$T_0 = 2456268.99119 \pm 7 \cdot 0.00012 BJD TDB$$
  <br/><br/>
  <br/><br/>
  ![alt text](/images/WASP-012_Sector_20_a_TimeSeries.png)
</div>

<div id="Data" class="tabcontent" markdown="1">
  {% include data/data_CoRoT-01.md %}
</div> 
 
<div id="Figures" class="tabcontent" markdown="1">
  {% include figures/figures_wasp-12.md %}
</div>


<script src="../js/tabs.js"></script> 
