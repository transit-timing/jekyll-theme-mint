---
    layout: post
    style: planet
    ---
    <script src="../js/planets.js"></script> 

    ## CoRoT-09

    <!-- Tab links -->
    <div class="tab">
      <button class="tablinks" onclick="openCity(event, 'Ephemeris')">Ephemeris</button>
      <button class="tablinks" onclick="openCity(event, 'Data')">Data</button>
      <button class="tablinks" onclick="openCity(event, 'Figures')">Figures</button>
    </div>

    <!-- Tab content -->
    <div id="Ephemeris" class="tabcontent" markdown="1">
    <br/><br/>
      $$P=1.091419108 \pm 95.272725 \cdot 6e-05 day$$
      $$T_0 = 2454603.34543 \pm 7 \cdot 0.0006 BJD TDB$$
      <br/><br/>
      <br/><br/>
      ![alt text](/images/WASP-012_Sector_20_a_TimeSeries.png)
    </div>

    <div id="Data" class="tabcontent" markdown="1">|    |   Mid-point |   Uncertainty | Time System   | #   | Ref                 |
|---:|------------:|--------------:|:--------------|:----|:--------------------|
|  0 | 2.4546e+06  |       0.0001  | BJD_TDB       | 1   |                     |
|  1 | 2.45537e+06 |       0.00037 | BJD_TDB       | >1  | 2017A&A...603A..43B |
    </div> 
     
    <div id="Figures" class="tabcontent" markdown="1">
    {% include figures/figures_wasp-12.md %}
    </div>


    <script src="../js/tabs.js"></script> 

     