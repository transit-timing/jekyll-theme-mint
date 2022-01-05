---
layout: post
---

<script src="../js/planets.js"></script> 

## Earth
Earth's columns are here:
<div id="planetdata"></div>

<script>
  var myData = planetData['Earth'].join(' | ');
  console.log(myData);
  document.getElementById("planetdata").innerHTML = myData;
</script>