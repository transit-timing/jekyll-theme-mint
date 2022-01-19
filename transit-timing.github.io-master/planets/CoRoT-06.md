---
    layout: post
    style: planet
    ---
    <script src="../js/planets.js"></script> 

    ## CoRoT-06

    <!-- Tab links -->
    <div class="tab">
      <button class="tablinks" onclick="openCity(event, 'Ephemeris')">Ephemeris</button>
      <button class="tablinks" onclick="openCity(event, 'Data')">Data</button>
      <button class="tablinks" onclick="openCity(event, 'Figures')">Figures</button>
    </div>

    <!-- Tab content -->
    <div id="Ephemeris" class="tabcontent" markdown="1">
    <br/><br/>
      $$P=1.091419108 \pm 8.8866287 \cdot 6e-06 day$$
      $$T_0 = 2457190.5078 \pm 7 \cdot 0.0013 BJD TDB$$
      <br/><br/>
      <br/><br/>
      ![alt text](/images/WASP-012_Sector_20_a_TimeSeries.png)
    </div>

    <div id="Data" class="tabcontent" markdown="1">|    |   Mid-point |   Uncertainty | Time System   | #   | Ref                 |
|---:|------------:|--------------:|:--------------|:----|:--------------------|
|  0 | 2.4546e+06  |       0.0002  | BJD           | >1  | 2010A&A...512A..14F |
|  1 | 2.45539e+06 |       0.00085 | BJD_TDB       | 1   | 2019EPSC...13..595E |
|  2 | 2.4554e+06  |       0.00074 | BJD_TDB       | 1   | 2019EPSC...13..595E |
|  3 | 2.45719e+06 |       0.00698 | BJD_TDB       | 1   | 2019EPSC...13..595E |
|  4 | 2.45832e+06 |       0.00102 | BJD_TDB       | 1   | 2019EPSC...13..595E |
|  5 | 2.45866e+06 |       0.00242 | BJD_TDB       | 1   | 2019EPSC...13..595E |
    </div> 
     
    <div id="Figures" class="tabcontent" markdown="1">
    {% include figures/figures_wasp-12.md %}
    </div>


    <script src="../js/tabs.js"></script> 

     