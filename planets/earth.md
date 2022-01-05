---
layout: post
planetdata: true
---

## Earth
Earth's columns are here:
<div id="planetdata"></div>

<script>
  var myData = planetData['Earth'].join(' | ');
  console.log(myData);
  document.getElementById("planetdata").innerHTML = myData;
</script>