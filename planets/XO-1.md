---
layout: post
style: planet
---
<script src="../js/planets.js"></script>

## XO-1

<!-- Tab links -->
<div class="tab">
<button class="tablinks" onclick="openCity(event, 'Ephemeris')">Ephemeris</button>
<button class="tablinks" onclick="openCity(event, 'Data')">Data</button>
<button class="tablinks" onclick="openCity(event, 'Figures')">Figures</button>
</div>

<!-- Tab content -->
<div id="Ephemeris" class="tabcontent" markdown="1">
<br/><br/>
$$P = 3.94150465(41) $$ day <br/>
$$T_0 = 2456055.57566(16) $$ BJD TDB
<br/><br/>
<br/><br/>
![alt text](/images/XO-1_o_c.png)
</div>


<div id="Data" class="tabcontent" markdown="1">

{% include data/data_XO-1.md %}

</div>

<div id="Figures" class="tabcontent" markdown="1">
{% include figures/figures_XO-1.md %}
</div>


<script src="../js/tabs.js"></script>


