// pcplot.js
// Input: mds
//          PC_<N>, CC_<N>, CC<K>_<N>
//        distLegend
//        pageName
// Usage: data -> cgi (invokes mds/pca) -> pcplot.js -> mds.cgi (invokes mds/pca) -> pcplot.js ...

// Constant
var is_pca = mds.distType === undefined;
  // false <=> MDS
var is_psd =  mds.eigens.psd;
var Class = mds.class || ""; 



var qualName = null;


var explainedTab = function ()
// Output: mdsCorrection, qualName
{
  var dec = 2;  // PAR
  qualName = "&lambda;" + ifS (! is_psd, "<sup>2</sup>");
  var tab = tableStart ("Eigenvalues") + "\
<tr>\
<th><u><i>i</i></u></th>\
<th>" + qualName + ",%</th>\
<th>&Sigma; ,%</th>\
</tr>";
  var total = 0;
  for (var i = 0; i < mds.eigens.explained.length; i++)
  {
    var qual = mds.eigens.explained[i];
    total += Math.abs(qual);
    var color = qual >= 0 ? "black" : "red";
    tab += "\
<tr>\
<td align=center>" + (i + 1) + "</td>\
<td align=center><font color=" + color + ">" + ifS(qual < 0, "(-)") + (Math.abs(qual) * 100).toFixed(dec) + "</font></td>\
<td align=center>" + (total * 100).toFixed(dec) + "</td>\
</tr>";
  }
  if (mds.eigens.explainedFrac_next !== undefined)
    tab += "<tr><td align=center>Next</td><td align=center>" + (mds.eigens.explainedFrac_next * 100).toFixed(dec) + "</td><td><br></td></tr>";
  tab += "</table>";
  if (! is_psd)
    tab += "&Sigma; &lambda; = " + (mds.eigens.explainedVarianceFrac * 100).toFixed(dec) + " %<br>";
  if (mds.eigens.orthogonal == false)
    tab += "<font color=red>Not orthogonal</font><br>";
      
  document.writeln (tab);
}



var explainedPlot = function (x_log /*Boolean*/)
{
  var lambdaPlot = createPlot ("lambdaPlot" + (x_log ? 1 : 0), 0.25, true);
  
  lambdaPlot. clear ();
  lambdaPlot. x_log = x_log;  
  lambdaPlot. y_log = true;
  lambdaPlot. drawAxes ( 1
                       , mds.eigens.explained.length
                       , Math.abs(mds.eigens.explained [mds.eigens.explained.length - 1] * 100)
                       , Math.abs(mds.eigens.explained [0] * 100)
                       , "\uD835\uDC56"
                       , "\u03BB\u00B2,%"
                       );
  lambdaPlot.ctx.beginPath ();
  lambdaPlot.ctx.strokeStyle = 'green';
  for (var i = 0; i < mds.eigens.explained.length; i++)
  {
    var qual = Math.abs(mds.eigens.explained[i]) * 100;
    var x = lambdaPlot.x_world2canvas (i + 1);
    var y = lambdaPlot.y_world2canvas (qual);
    if (i)
      lambdaPlot.ctx.lineTo (x, y)
    else
      lambdaPlot.ctx.moveTo (x, y)
  }
  lambdaPlot.ctx.stroke ();

  var lr = mds.eigenValueLR;
  if (x_log && lr)
  {
    var tab = qualName + "(<i>i</i>) = " + (lr.intercept * 100).toFixed(1) + " <i>i</i> <sup>" + lr.coeff + " </sup> %<br>";
    tab += "R<sup>2</sup> = " + lr.log_scale_R2 + " (in log-log scale)<br>";
    document.writeln (tab);
  }
}



//----------------------------------- DOM -----------------------------------
if (! is_pca)
  if (mds.distAttr)
  {
    var distTypeS = "";
    switch (mds.distType)
    {
      case 0: distTypeS = "similarity"; break;
      case 1: distTypeS = "distance"; break;
      case 2: distTypeS = "squared distance"; break;
    }
    document.writeln ("MDS of attribute '" + mds.distAttr + "' as " + distTypeS + "<br>");
  }
  else
  {
    if (mds.distType)
      alert ("Wrong type: " + mds.distType);
    document.writeln ("PCA via MDS<br>");
  }
  
if (mds.comments)
{
//document.writeln ("<b>Problems:</b><br>");
  for (var i in mds.comments)
    document.writeln (mds.comments[i] + "<br>");    
}

var plot = createPlot ("Plot", 0.6, true);
plot.canvas.addEventListener('mousedown', function (evt) {pos2node (evt);}, false);  

document.writeln ("<pre>");
if (distLegend)
  document.writeln ("<b>Scale:</b> " + distLegend);
document.writeln ("# Points: " + mds.objs.length);
if (mds.outlier_evalue)
  document.writeln ("# Outlier max. E-value: " + mds.outlier_evalue);
document.writeln ("</pre>");

document.writeln ("<table><tr><td><a href=\"javascript:canonical=false;separationClusterNum=0;func()\"><form>Reset</form></a></td>");
document.writeln ("<td><br></td><td><form onchange='parent.showNames=showNames.checked;parent.func();return false'>Show names: <input type=\"checkbox\" id=\"showNames\"></form></td>");
if (mds.objs[0].attr["CC_1"] !== undefined)
  document.writeln ("<td><br></td><td><form onchange='parent.canonical=canonical.checked;parent.func();return false'>Make clusters separate: <input type=\"checkbox\" id=\"canonical\"></form></td>");
document.writeln ("</tr></table>");

if (mds.attr_mds)
  document.writeln (" <form action='attr_mds.cgi' method=post target=_blank>"
                       + "<input type=\"hidden\" name=\"attr_mds\" value='" + q2h(window.JSON.stringify(mds.attr_mds)) + "'>"
                       + "<input type=\"hidden\" name=\"name\" value='" + q2h(pageName) + "'>"
                       + "<input type=submit value='MDS of Attributes'>"
                       + "</form>");
// TODO ??
// button "Convert to PNG" 
// select: dimension 1, dimension 2

document.writeln ("<table><tr><td valign=top>");
explainedTab ();
document.writeln ("</td><td><br></td><td valign=top><br><br>");
explainedPlot (true);
document.writeln ("<br>");
explainedPlot (false);
document.writeln ("</td><td><br></td><td valign=top><object id=\"clustersTab\"><br></object></td></tr></table>");

document.write ("<p><object id=\"pointsTab\"><br></object>");
//---------------------------------------------------------------------------

  


var getMin = function (attr)
{
  var x = Infinity;
  for (var i in mds.objs)
  {
    var v = mds.objs[i].attr[attr];
    if (x > v)
      x = v;
  }
  return x;
}



var getMax = function (attr)
{
  var x = -Infinity;
  for (var i in mds.objs)
  {
    var v = mds.objs[i].attr[attr];
    if (x < v)
      x = v;
  }
  return x;
}



// Global
var points = [];  // {objName, class, comment, x, y}, objName is unique
var classes = {};  // <class>:1
var classNum = 0;  // = size of classes

var showNames = false;
var canonical = false;
var separationClusterNum = 0; 
  // 0 - all


var pointName = function (point)
{
  return point.objName + sPref (": ", (classNum >= 2 ? "[" + (point.class || "?") + "] " : "") + (point.comment || ""));
}



var pos2node = function (evt)
{
  if (evt.button)  // Mouse non-left button (IE: >=9)
    return;
  
  var pos = getMouseCanvasPos (evt, plot.canvas);

  var innerHTML = "";
  for (var i in points)
  {
    var point = points [i];
    if (   Math.abs (pos.x - point.x) <= 2 * plot.dX
        && Math.abs (pos.y - point.y) <= 2 * plot.dX
       )
      innerHTML += pointName (point) + "\n";
  }
  
  plot.showMenu (innerHTML, "Click on a point, please!", evt);
}



var markPoint = function (objName)
{
  for (var i in points)
  {
    var point = points [i];
    if (point. objName == objName)
      plot.markXY (point. x, point. y);
  }
}



var func = function ()
{
  var meanSd = function (mean, sd, decimals)
  {
    var sdS = sd.toFixed (decimals + 1);
    if (sdS.match(/^0.0*$/))
      return mean.toFixed(decimals);
    return mean.toFixed (decimals + 1) + ' &plusmn; ' + sdS;
  }
  

  var dimPref1 = "PC_";
  var dimPref2 = "CC_";
  var dimPref = canonical 
                  ? separationClusterNum 
                    ? "CC" + separationClusterNum + "_"
                    : dimPref2
                  : separationClusterNum 
                    ? "BC" + separationClusterNum + "_"
                    : dimPref1;


  var dim = 0;
  while (mds.objs[0].attr[dimPref + (dim + 1)] !== undefined)
    dim++;
  if (! canonical && ! separationClusterNum && dim != mds.eigens.explained.length)
    alert ("Bad dim = " + dim);
    

  // Init
  plot.clear ();
  plot.drawAxes ( getMin (dimPref + "1")
                , getMax (dimPref + "1")
                , dim < 2 ? 0 : (separationClusterNum ? 0 : getMin (dimPref + "2"))
                , dim < 2 ? 0 : getMax (dimPref + "2")
                , null
                , null
                );

  points = []; 
  classes = {}; 
  classNum = 0;  
  
  // attrs
  var attrs = [];
  if (mds.objs.length)
    for (var s in mds.objs[0].attr)
      if (   s != "Cluster"
        //&& s != "Class"
          && s.search(dimPref1) != 0
          && s.search(dimPref2) != 0
          && s.search(/CC\d+_/) != 0
          && s.search(/BC\d+_/) != 0
         )
        attrs.push (s);
  
  var z_min = dim < 3 ? 0 : getMin (dimPref + "3");
  var z_max = dim < 3 ? 0 : getMax (dimPref + "3");
  
  var pointsTab = tableStart ("Points") + "\
<tr>\
<th><u>Name</u></th>";
  for (var i in attrs)
    pointsTab += "<th>" + attrs[i].replace(/_/g, " ") + "</th>";
pointsTab += "\
<th>Comment</th>\
<th>X</th>\
<th>Y</th>\
<th>Z</th>";
  if (mds.clusters)
    pointsTab += "<th>Cluster</th>";
  pointsTab += "</tr>";
  var clusters = [];  // {name, x, y, count, names:{<name>:<count>}}
  if (mds.clusters)
    for (var i = 0; i < mds.clusters.length; i++)
      clusters.push ({name:mds.clusters[i].name, x:0, y:0, count:0, names:{}});
  for (var i in mds.objs)
  {
    var obj = mds.objs[i];
    var attr = obj.attr;
    var X_1 = dim >= 1 ? attr[dimPref + "1"] : 0;
    var X_2 = dim >= 2 ? attr[dimPref + "2"] : 0;
    var X_3 = dim >= 3 ? attr[dimPref + "3"] : 0;
  /*if (separationClusterNum)
      X_2 *= Math.random(); */
    var x = plot.x_world2canvas (X_1);
    var y = plot.y_world2canvas (X_2);
    var z = dim < 3 ? null : (X_3 - z_min) / (z_max - z_min);
    var clusterNum = 0;
    for (var j = 0; j < clusters.length; j++)
      if (clusters [j]. name == attr.Cluster)
      {
        clusterNum = j;
        break;
      }
    plot.drawPoint (x, y, clusterNum, z || 0);
    var point = { objName: obj.objName
                , class:   Class ? attr[Class] : null
                , comment: obj.comment
                , x:x
                , y:y
                };
    points.push (point);
    if (showNames)
      plot.ctx.fillText (pointName (point), x + plot.dX, y);
    // cluster
    var cluster = null;
    if (clusters.length)
    {
      cluster = clusters[clusterNum];
      cluster.x += x;
      cluster.y += y;
      cluster.count++;
    }
    var dataClass = point.class;  
    if (dataClass)
    {
      if (clusters.length)
        if (cluster.names[dataClass])        
          cluster.names[dataClass] ++;
        else
          cluster.names[dataClass] = 1;
      // classes, classNum
      if (! classes[dataClass])
      {
        classes[dataClass] = 1;
        classNum++;
      }
    }
    // pointTab
    pointsTab += "<tr><td><a href=\"javascript:markPoint('" + q2s(point.objName) + "')\">" + point.objName + "</a></td>";
    for (var j in attrs)
    {
      var val = attr[attrs[j]];
      pointsTab += "<td align=center>" + (val == null ? "<br>" : val) + "</td>";
    }
    pointsTab +=   "<td>" + (point.comment || " ") 
                 + "</td><td align=right>" + X_1 
                 + "</td><td align=right>" + X_2 
                 + "</td><td align=right>" + X_3 
                 + "</td>";
    if (mds.clusters)
      pointsTab += "<td align=center><font color=" + int2fillStyle (clusterNum) + "><b>" + attr.Cluster + "</font></td>";
    pointsTab += "</tr>";
  }
  document.getElementById("pointsTab").innerHTML = pointsTab;

  // clusters
  if (mds.clusters && mds.clusters.length >= 2)
  {
    var clustersTab = tableStart ("Clusters") + "\
<tr>\
<th><u>Cluster</u></th>\
<th># Points</th>\
<th>" + (Class || "Class") + "</th>\
<th>Characteristic attributes:<br>avg. in cluster &plusmn; SD vs.<br>avg. outside cluster &plusmn; SD</th>\
<th>Space</th>\
<th>Separation</th>\
</tr>";
    for (var i = 0; i < mds.clusters.length; i++)
    {
      var cluster = clusters[i];
      cluster.x /= cluster.count;
      cluster.y /= cluster.count;
      var clusterNamesShort = "";  
      var clusterNamesLong = "";  
      if (cluster.names)
      {   
        var mainClass = null;
        var mainClassCount = 0;
        for (var j in cluster.names)
          if (mainClassCount < cluster.names[j])
          {
            mainClassCount = cluster.names[j];
            mainClass = j;
          }
        for (var j in cluster.names)
        {
          var goodClass = cluster.names[j] >= 0.5 * mainClassCount;  // PAR
          if (goodClass)  
            clusterNamesShort += j + " "; 
          clusterNamesLong += (goodClass ? "<b>" : "") + j + (goodClass ? "</b>" : "") + (arraySize (cluster.names) > 1 ? "(" + cluster.names[j] + ")" : "") + " ";
        }
        plot.fillTextCenter (clusterNamesShort, cluster.x, cluster.y);  
        if (separationClusterNum)
        {
          plot.ctx.beginPath();
          plot.ctx.moveTo(cluster.x, plot.y_world2canvas(plot.y_plot_min));
          plot.ctx.lineTo(cluster.x, plot.y_world2canvas(plot.y_plot_max));
          var strokeStyle = plot.ctx.strokeStyle;
          plot.ctx.strokeStyle = int2fillStyle (i);
          plot.ctx.stroke ();
          plot.ctx.strokeStyle = strokeStyle;
        }
      }
      var dm = mds.clusters[i].dm;
      var charAttrs = "";
      var dependencies = mds.clusters[i].dependencies;
      if (dependencies)
        for (var j = 0; j < dependencies.length; j++)
        {
          var dep = dependencies[j];
          charAttrs += dep.attr + ': ' + meanSd (dep.cluster.mean, dep.cluster.sd, dep.decimals)
                                       + ' vs. ' 
                                       + meanSd (dep.other.mean, dep.other.sd, dep.decimals)
                                       + "<br>";
        }
      var separationVal = canonical ? mds.clusters[i].canonical : mds.clusters[i].between_center_d2;
      clustersTab +=   "<tr>" 
                     + "<td align=center><font color=" + int2fillStyle (i) + "><b>" + cluster.name  + "</td>"
                     + "<td align=center>" + cluster.count + "</td>"
                     + "<td>" + (clusterNamesLong || "<br>") + "</td>"
                     + "<td>" + (charAttrs || "<br>") + "</td>"
                     + "<td>" + (dm ? 
                        ("<form action='mds.cgi' method=post target=_blank style='display:inline;'>"
                       + "<input type=\"hidden\" name=\"dm\" value='" + q2h(dm) + "'>"
                       + "<input type=\"hidden\" name=\"class\" value='" + q2h(Class) + "'>"
                       + "<input type=\"hidden\" name=\"name\" value='" + q2h(pageName) + "/Cluster " + (i + 1) + "'>"
                       + ifS (! is_pca, "<input type=\"hidden\" name=\"distAttr\" value='" + mds.distAttr + "'>")
                       + ifS (! is_pca, "<input type=\"hidden\" name=\"distType\" value='" + mds.distType + "'>")
                       + "<input type=\"hidden\" name=\"distLegend\" value='" + q2h(distLegend) + "'>"
                       + "<input type=submit value=Space>"
                       + "</form>") : "")
                       + "</td>"
                     + "<td align=center>" + (separationVal !== undefined ? "<a href=\"javascript:separationClusterNum=" + (i + 1) + ";func()\">" + separationVal + "</a>" : "<br>") + "</td>"
                     + "</tr>";
    }
    clustersTab += "</table>";
    document.getElementById("clustersTab").innerHTML = clustersTab; 
  }
}
func ();



