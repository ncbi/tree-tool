// common.js


// Global constants
var debug = false;
// PAR
var fontHeight = 14;  
var fontName = "Lucida Console";
var fillStyle = 'black';



onerror = function (msg,url,l)
{
  alert ("\
There was an error on this page.\n\n\
Error: " + msg + "\n\
URL: " + url + "\n\
Line: " + l + "\n\n\
Click OK to continue.\n\n\
");
  return true;
}



var printProperties = function (x)
// Reflection
{
  let offset = 0;

  var printProperties_ = function (name/*string*/, x)
  {
    document.writeln ();
    for (let j = 0; j < offset; j++)
      document.write ("  ");
    if (name)
      document.write (name + ": ");
    const type = typeof x;
    document.write (type + ": ");
    if (type == "object")
    {
      offset++;
      for (const i in x)
      //if (x.hasOwnProperty(i))
          printProperties_ (i, x[i]);
      offset--;
    }
    else
      document.write (x);
  }
  
  document.write ("<pre>");
  printProperties_ (null, x);
  document.writeln ("</pre>");
}



var inputString = function (promptS, val_init)
{
  let s = prompt (promptS, val_init || '');
  if (s == null)
    return null;
  s = s.trim ();
  return s;
}



var getMax = function (a, b)
{
  if (a == null)
    return b;
  if (b == null)
    return a;
  if (a > b)
    return a;
  return b;
}



var getMin = function (a, b)
{
  if (a == null)
    return b;
  if (b == null)
    return a;
  if (a > b)
    return b;
  return a;
}



var minimize = function (x, field, y)
{
  if (y == null)
    return false;
  if (x [field] == null || x [field] > y)
  {
    x [field] = y;
    return true;
  }
  return false;
}



var maximize = function (x, field, y)
{
  if (y == null)
    return false;
  if (x [field] == null || x [field] < y)
  {
    x [field] = y;
    return true;
  }
  return false;
}



var isDigit = function (x)
{
  return    x.length == 1
         && x >= '0'
         && x <= '9';
}



var ifS = function (cond, s)
{
  return cond ? s : "";
}



var _2s = function (s)
{
  return ifS (s, s);
}



var nvl = function (s, nullS)
{
  return s ? s : nullS;
}



var sPref = function (pref, s)
{
  return s ? pref + s : "";
}


var addS = function (s, delimiter, addition)
{
  if (! addition)
    return s;
  if (s)
    return s + delimiter + addition;
  return addition;
}



var yesNo = function (flag, yes = "Yes", no = "No")  
{
  if (flag == null)
    return "?";
  return flag ? yes : no;
}


var yesNo_std = function (flag)  
{
  return yesNo (flag, "Yes", "No");
}


var str_color = function (str, cond, bold)
{
  if (! cond)
    return str;
  return "<font color=red>" + (bold ? "<b>" : "") + str + "</font>" + (bold ? "</b>" : "");
}



var strSpaces = function (s)
{
  let n = 0;
  for (let i = 0; i < s.length; i++)
    if (s [i] == ' ')
      n++;
  return n;
}



var q2s = function (s)
{
  if (! s)
    return s;
  return s.replace(/\'/g, "\\\'")
          .replace(/\"/g, '\\\"')
          ;
}



var q2h = function (s)
{
  if (! s)
    return s;
  return s.replace(/\'/g, "&#39;")
          .replace(/\"/g, "&#34;")
          .replace(/ /g,  "&#32;")
          .replace(/\t/g, "&#09;")
          .replace(/\n/g, "&#10;")
          .replace(/\r/g, "")   // For UNIX
          .replace(/</g,  "&lt;") 
          .replace(/>/g,  "&gt;") 
          ;
}



var q2cgi = function (s)
{
  if (! s)
    return s;
  return s.replace(/\'/g, "%27")
          .replace(/\"/g, "%22")
          .replace(/ /g,  "%20")
          .replace(/\t/g, "%09")
          .replace(/\n/g, "%0A")
          .replace(/\r/g, "")   // For UNIX
          .replace(/=/g,  "%3D")
          .replace(/&/g,  "%26")
          ;
}



var json = function (x)
{
  var toText = function (e)
  {
    if (e == null)
      return "null";
    if (typeof e == "string")
      return '"' + q2s(e) + '"';
    if (typeof e == "object")  // array
    {
      let s = "";
      for (const i in e)
      {
        if (s)
          s += ",";
        s += i;
      }
      return "[" + s + "]";
    }
    return e;
  }
  
  let s = "";
  for (const i in x)
  //if (typeof x[i] != "object")
    {
      if (s)
        s += ",";
      s += i + ":" + toText(x[i]);
    }
  return "{" + s + "}"; 
}
  
  

var hasMask = function (word)
{
  return    typeof word == "string" 
         && word 
         && (   word.indexOf('%') != -1
             || word.indexOf('_') != -1
             || word.indexOf('[') != -1
            );         
}



var numberWithCommas = function (x) 
{
  let parts = x.toString().split(".");
  parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  return parts.join(".");
}



var strContains = function (hay, needle, lowercaseP)
{
  if (! needle)
    return true;
  if (! hay)
    return false;
  if (lowercaseP)
  {
    hay    = hay.   toLowerCase();
    needle = needle.toLowerCase();
  }
  return hay.indexOf (needle) != -1;
}



var arrayEmpty = function (arr)
{
  for (const i in arr)
    return false
  return true;
}



var arraySize = function (arr)
{
  let n = 0;
  for (const i in arr)
    n++;
  return n;
}



var arrayContains = function (arr, obj) 
{
  for (const i in arr) 
    if (arr[i] == obj) 
      return true;
  return false;
}



var structArrayContains = function (arr, field, obj) 
{
  for (const i in arr) 
    if (arr[i] [field] == obj) 
      return true;
  return false;
}



var array2string = function (arr)
{
  let s = "";
  for (const i in arr) 
  {
    if (s)
      s += ", ";
    s += arr[i];
  }
  return s;  
}



var table_column2string = function (table, column)
// Return: concatenated strings prefixed with " "
{
  let s = "";
  for (const i in table)
    s += " " + table[i][column];
  return s;
}



var table_column2string_distinct = function (table, column)
// Return: concatenated strings prefixed with " "
{
  let map = {};
  for (const i in table)
  {
    const v = table[i][column];
    if (v)
      map [v] = 1;
  }

  let s = "";
  for (const i in map)
    s += " " + i;
  
  return s;
}



var int2hex = function (i)
// Input: i: between 0 and 255
{
  var hex = Math.round(i).toString(16);
  if (hex.length == 1)
    hex = '0' + hex;
  return hex;
}



var int2fillStyle = function (n)
{
  let n_ = (n % 6) + 1;  // n_ == 7 <=> white
    // 1..6
    
  let color_min = 0;
  let color_max = 255;
  switch ((Math.floor ((n + 0.1) / 6) % 3))
  {
    case 0: break;
    case 1: color_min = 128; break;
    case 2: color_max = 127; break;
  }
  
  const range = color_max - color_min;
  
  var lastDigitShift = function (prod)
  // Update: n_
  {
    const res = n_ % 2;
    n_ = Math.round (n_ / 2 - 0.1);
    return color_min + res * range;
  }

  const b = lastDigitShift ();
  const g = lastDigitShift ();
  const r = lastDigitShift ();

  return   '#' 
         + int2hex (r) 
         + int2hex (g) 
         + int2hex (b);  
}



const stnd_max = 3;  // PAR



var stnd2color = function (z /*standardized number*/)
{
  // z = 0 => "black"
  z = Math.max (Math.min (z, stnd_max), -stnd_max);
  const red   = Math.max(-z,0) / stnd_max * 255;
  const green = Math.max( z,0) / stnd_max * 255;
  const red_hex   = int2hex (red);
  const green_hex = int2hex (green);
  const blue_hex  = int2hex (Math.min (red, green));  // PAR  
    // int2hex (Math.max (Math.min(255,red / 0.5), green / 4));  // PAR  
  return   '#' 
         + red_hex
         + green_hex  
         + blue_hex;  
}



var frac2color = function (frac /*0..1*/)
{
  const i = (1 - frac) * 255;
  return   '#' 
         + int2hex (i)  // R
         + int2hex (i)  // G
         + int2hex (i); // B
}



var tableStart = function (caption, widthPercent, fontSize)
{
  const wp = widthPercent ? " style='table-layout:fixed;width:" + widthPercent + "%;'" : "";
  const fs = fontSize ? " style='font-size:" + (fontHeight * fontSize) + "px'" : "";
  return   "<table border=1 cellpadding=5 style='border-collapse:collapse;'" + wp + fs + ">"
         + ifS (caption, "<caption><b><big>" + caption + "</caption>");
}



var column2table = function (row, column, isKey)
{
  const v = row[column];
  if (v == '----')
    return "<tr><td colspan=2 align=center>" + column + "</td></tr>";
  return "<tr><td><b>" + column + "</b></td><td>" + ifS(isKey,"<b><font color=blue>") + (v == null ? "" : v) + ifS(isKey,"</b></color>") + "</td></tr>";
}



var queryRow2table = function (title, row, keyNum)
{
  let tab = tableStart (title) + "<tbody>";
  let i = 0;
  for (const col in row)
  {
    tab += column2table (row, col, i < keyNum);
    i++;
  }
  return tab + "</tbody></table>";
}



var setObject = function (objName, html)
{
  const obj = document.getElementById (objName);
  if (! obj)
  {
  	alert ("DOM object " + objName + " is not found");
  	return null;
  }
  obj.innerHTML = html;
  return obj;
}



var setFocusObject = function (objName, html)
{
  const obj = setObject (objName, html);
  if (obj)
    obj.scrollIntoView();
}



var table2DOM = function (tableName, tableHTML)
{
  setFocusObject (tableName, 
    "<p>" + tableHTML + "<form onsubmit='parent.remove_" + tableName + "();return false'><input type='submit' value='Remove'></form><p>");
}



var getMouseCanvasPos = function (event, canvas)
// http://miloq.blogspot.com/2011/05/coordinates-mouse-click-canvas.html
{
  let x, y;

  if (event.x != undefined && event.y != undefined)
  {
    x = event.x;
    y = event.y;
    x += window.pageXOffset;
    y += window.pageYOffset;
  }
  else // "Firefox"
  {
    x = event.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
    y = event.clientY + document.body.scrollTop  + document.documentElement.scrollTop;
  }

  x -= canvas.offsetLeft;
  y -= canvas.offsetTop;

  return {x:x, y:y};
}



var getParams = function ()
// Return: {<name>:<val>}
{
  let params = {};
  const search = window.location.search;
  if (search)
    if (search.substring(0,1) == "?")
    {
      const pars = search.substring(1).split("&");
      for (const i in pars)
      {
        const pair = pars[i].split("=");
        if (pair.length == 2)
          params [pair[0]] = unescape(pair[1]);
        else
          alert ("Bad parameter " + pars[i]);
      }
    }
    else
      alert ("Bad parameters: " + search);

  return params;
}



var getNow = function ()
{
  const now = new Date ();
  return (now.getMonth() + 1) + '/' + now.getDate() + '/' + now.getFullYear() + "  " + String(now.getHours()).padStart(2,'0') + ":" + String(now.getMinutes()).padStart(2,'0');
}



var createPlot = function (plotName, innerHeightFrac_arg, square_arg/*Boolean*/)
// Return: object with usage: {clear(), drawAxes(); ...}*
{
  document.writeln("<p><canvas id='" + plotName + "' style=\"border:1px solid #000000;\"></canvas>");
  
  
  var plot =
  {
    // init
    canvas: null,
    ctx: null,
    canvasMenu: null,  
    dX: null,
    innerHeightFrac: innerHeightFrac_arg,
    square: square_arg,
    
    // set manually
    x_unit_min: 0,  
    y_unit_min: 0,
    x_log: false,
    y_log: false,
    
    // automatic
    canvasStartX: null,
    canvasEndX: null,
    canvasStartY: null,
    canvasEndY: null,
    
    x_plot_min: null,  
    x_plot_max: null,  
    x_range: 0,
    x_origin: null,  
  
    y_plot_min: null,  
    y_plot_max: null,  
    y_range: 0,
    y_origin: null,
  
  
    clear: function ()
    {
      this.ctx.clearRect (0, 0, this.canvas.width, this.canvas.height);  
  
      this.canvasStartX = 2 * this.dX;
      this.canvasStartY = 2 * this.dX;
      this.canvasEndX = this.canvas.width  - 2 * this.dX;
      this.canvasEndY = this.canvas.height - 2 * this.dX;    
  
      this.x_unit_min = 0;
      this.y_unit_min = 0;
      
      this.x_plot_min = null;
      this.x_plot_max = null;
      this.x_range = 0;
      this.x_origin = null;
      
      this.y_plot_min = null;
      this.y_plot_max = null;
      this.y_range = 0;
      this.y_origin = null;
  
      this.canvasMenu.style.visibility = "hidden";    
    },
    
    
    resize: function ()
    {
      this.canvas.height = window.innerHeight * this.innerHeightFrac;  
      this.canvas.width  = this.square ? this.canvas.height : window.innerWidth * 0.97;  // PAR
    },
  
  
    x_world2canvas: function (x_world)  
    { 
      const rel = this.x_range > 0 
                    ? (this.x_log 
                        ? Math.log(x_world) - Math.log(this.x_plot_min)
                        :          x_world  -          this.x_plot_min
                      ) / this.x_range 
                    : 0.5;
      return this.canvasStartX + rel * (this.canvasEndX  - this.canvasStartX); 
    },
    
    
    y_world2canvas: function (y_world)  
    { 
      const rel = this.y_range > 0 
                    ? (this.y_log
                         ? Math.log(y_world) - Math.log(this.y_plot_min)
                         :          y_world  -          this.y_plot_min
                      ) / this.y_range 
                    : 0.5;
      return this.canvasStartY + (1 - rel) * (this.canvasEndY - this.canvasStartY); 
    },
    
    
    x_canvas2world: function (x)
    {
      const rel = (x - this.canvasStartX) / (this.canvasEndX  - this.canvasStartX);
      if (this.x_range <= 0)
        return null;
      x = (this.x_log ? Math.log(this.x_plot_min) : this.x_plot_min) + rel * this.x_range;
      if (this.x_log)
        return Math.exp(x);
      return x;
    },
    
  
    y_canvas2world: function (y)
    {
      const rel = (y - this.canvasStartY) / (this.canvasEndY - this.canvasStartY);
      if (this.y_range <= 0)
        return null;
      y = (this.y_log ? Math.log(this.y_plot_min) : this.y_plot_min) + rel * this.y_range;
      if (this.y_log)
        return Math.exp(y);
      return y;
    },
    
    
    drawVertLine: function (x_world, arrow/*Boolean*/, legend, color)
    {
      if (x_world == null)
        return;
      const x = this.x_world2canvas(x_world);
      const y = this.y_world2canvas(this.y_plot_max);
      this.ctx.beginPath ();
      this.ctx.moveTo (x, this.y_world2canvas(this.y_plot_min));
      this.ctx.lineTo (x, y);
      if (arrow)
      {
        this.ctx.lineTo (x - this.dX, y + 2 * this.dX);
        this.ctx.moveTo (x, y);
        this.ctx.lineTo (x + this.dX, y + 2 * this.dX);
      }
      const strokeStyle = this.ctx.strokeStyle;
      this.ctx.strokeStyle = color;
      this.ctx.stroke ();
      this.ctx.strokeStyle = strokeStyle;
      
      if (legend)
      {
        const font_old = this.ctx.font;
        this.ctx.font = (fontHeight * 1.1) + "px " + fontName;  // PAR
        const x0 = x - 3 * this.dX;
        const y0 = this.canvas.height / 2;
        this.ctx.translate (x0, y0);
        this.ctx.rotate (-Math.PI/2);
        this.ctx.fillStyle = color;
        this.fillTextCenter (legend, 0, 2 * this.dX);
        this.ctx.rotate (Math.PI/2);
        this.ctx.translate (- x0, - y0);
        this.ctx.font = font_old;
        this.ctx.fillStyle = fillStyle;
      }
    },
    
  
    drawAxes: function (x_min, x_max, y_min, y_max, xLegend, yLegend)
    // Output: x_origin, y_origin
    {
      let txt
        , x, x_world
        , y, y_world
        , w;
      
      var range2unit = function (range)
      // Return: {unit, dec}
      {
        if (range <= 0)
          return {unit:null, dec:null};
        const x = [2, 5, 10];
        let dec = Math.floor (Math.log(range) * Math.LOG10E) - 1;
        const p = Math.pow (10, dec); 
        let unit_best;
        let k_max = -1;
        for (const i in x)
        {
          const unit = x[i] * p;
          const k = Math.round (range / unit);
          if (k <= 10)  // PAR
            if (k_max < k)
            {
              k_max = k;
              unit_best = unit;
              if (x[i] == 10)
                dec++;
            }
        }
        return {unit:unit_best, dec: Math.max(0,-dec)};
      }  
      
      var range_unit2origin = function (min, max, unit, log)
      {
        if (min <= 0 && max >= 0)
          return 0;
        if (min > 0)
        {
          if (log && min < unit)
            return unit;
          return Math.floor (min / unit) * unit;
        }
        return -Math.floor (-max / unit) * unit;
      }
      
      var unit2min = function (min, max, unit)
      {
        const range = max - min;
        if (range < unit)
          return min + range / 2 - unit;
        return min;
      }
      
      var unit2max = function (min, max, unit)
      {
        const range = max - min;
        if (range < unit)
          return min + range / 2 + unit;
        return max;
      }
      
      if (x_min == x_max)
        return;
      
      this.x_plot_min = unit2min (x_min, x_max, this.x_unit_min);
      this.x_plot_max = unit2max (x_min, x_max, this.x_unit_min);
      this.y_plot_min = unit2min (y_min, y_max, this.y_unit_min);
      this.y_plot_max = unit2max (y_min, y_max, this.y_unit_min);
  
      const x_unit = range2unit (this.x_plot_max - this.x_plot_min);
      const y_unit = range2unit (this.y_plot_max - this.y_plot_min);
      
      x_unit.unit = Math.max (this.x_unit_min, x_unit.unit);
      y_unit.unit = Math.max (this.y_unit_min, y_unit.unit);
      
      this.x_origin = range_unit2origin (x_min - this.x_unit_min, x_max + this.x_unit_min, x_unit.unit, this.x_log);
      this.y_origin = range_unit2origin (y_min - this.y_unit_min, y_max + this.y_unit_min, y_unit.unit, this.y_log);
    
      // To make axes visible
      if (this.x_plot_min > this.x_origin)
          this.x_plot_min = this.x_origin;
      if (this.x_plot_max < this.x_origin)
          this.x_plot_max = this.x_origin;
      //
      if (this.y_plot_min > this.y_origin)
          this.y_plot_min = this.y_origin;
      if (this.y_plot_max < this.y_origin)
          this.y_plot_max = this.y_origin;
  
      this.x_range = this.x_log ? (Math.log(this.x_plot_max) - Math.log(this.x_plot_min)) : (this.x_plot_max - this.x_plot_min);
      this.y_range = this.y_log ? (Math.log(this.y_plot_max) - Math.log(this.y_plot_min)) : (this.y_plot_max - this.y_plot_min);
    
      
      if (xLegend)
        this.canvasEndY -= 4 * this.dX;
      if (yLegend)
        this.canvasStartX += 3 * this.dX;
      
      // x
      this.ctx.beginPath ();
      x = this.x_world2canvas(this.x_plot_max);
      y = this.y_world2canvas(this.y_origin);
      this.ctx.moveTo (this.x_world2canvas(this.x_plot_min), y);
      this.ctx.lineTo (x, y);
      // arrow
      this.ctx.lineTo (x - 2 * this.dX, y - this.dX);
      this.ctx.moveTo (x, y);
      this.ctx.lineTo (x - 2 * this.dX, y + this.dX);
      // ticks
      for (x_world = this.x_origin - x_unit.unit; x_world >= this.x_plot_min; x_world -= x_unit.unit)
      {
        x = this.x_world2canvas(x_world);
        this.ctx.moveTo(x,y-this.dX);
        this.ctx.lineTo(x,y+this.dX);
        txt = x_world.toFixed(x_unit.dec);
        w = this.ctx.measureText(txt).width / 2
        this.ctx.fillText (txt, x - w, y-2*this.dX); 
      }
      for (x_world = this.x_origin + x_unit.unit; x_world <= this.x_plot_max; x_world += x_unit.unit)
      {
        x = this.x_world2canvas(x_world);
        this.ctx.moveTo(x,y-this.dX);
        this.ctx.lineTo(x,y+this.dX);
        txt = x_world.toFixed(x_unit.dec);
        w = this.ctx.measureText(txt).width / 2
        this.ctx.fillText (txt, x - w, y-2*this.dX); 
      }
      this.ctx.stroke();
      if (xLegend)
      {
        const font_old = this.ctx.font;
        this.ctx.font = (fontHeight * 1.1) + "px " + fontName;  // PAR
        this.fillTextCenter (xLegend, this.canvas.width / 2, this.canvas.height - 2 * this.dX);
        this.ctx.font = font_old;
      }
    
      // y
      if (this.y_range > 0)
      {
        this.drawVertLine (this.x_origin, true, yLegend, "black");
        // ticks
        this.ctx.beginPath ();
        x = this.x_world2canvas(this.x_origin);
        for (y_world = this.y_origin - y_unit.unit; y_world >= this.y_plot_min; y_world -= y_unit.unit)
        {
          y = this.y_world2canvas(y_world);
          this.ctx.moveTo (x-this.dX, y);
          this.ctx.lineTo (x+this.dX, y);
          this.ctx.fillText (y_world.toFixed(y_unit.dec), x + 2 * this.dX, y + this.dX / 2); 
        }
        for (y_world = this.y_origin + y_unit.unit; y_world <= this.y_plot_max; y_world += y_unit.unit)
        {
          y = this.y_world2canvas(y_world);
          this.ctx.moveTo (x-this.dX, y);
          this.ctx.lineTo (x+this.dX, y);
          this.ctx.fillText (y_world.toFixed(y_unit.dec), x + 2 * this.dX, y + this.dX / 2); 
        }
        this.ctx.stroke();
      }
    },
  
  
    markXY: function (x, y)
    {
      const d = 1.5 * this.dX;  // PAR
      this.ctx.beginPath ();
      this.ctx.moveTo (x - d, y - d)
      this.ctx.lineTo (x + d, y + d);
      this.ctx.moveTo (x - d, y + d)
      this.ctx.lineTo (x + d, y - d);
      this.ctx.stroke ();
    },
  
  
    drawPoint: function (x, y, fillStyleNum, z/*[0,1]*/)
    {
      const delta = 0.3;  // PAR
      z = z == null ? 1 : (delta + (1 - z) * (1 - delta));  // z in [delta,1]
      this.ctx.beginPath ();
      this.ctx.arc (x, y, this.dX * z, 0, /* z * */ 2 * Math.PI);  
      this.ctx.fillStyle = int2fillStyle (fillStyleNum);
      this.ctx.fill ();
      this.ctx.fillStyle = fillStyle;
    },
  
  
    fillTextCenter: function (txt, x, y)
    {
      const w = this.ctx.measureText(txt).width;  // != ""
      x = Math.min(this.canvas.width - w, Math.max(0, x - w / 2));      
      this.ctx.fillText (txt, x, y);  
    },
    
    
    closeMenu: function ()
    {
      this.canvasMenu.style.visibility = "hidden";    
    },
    
    
    showMenu: function (innerHTML, errorMsg, evt)
    {
      const visibilityOld = this.canvasMenu.style.visibility;
      this.closeMenu();
    
      if (innerHTML)
        innerHTML = "<pre>" + innerHTML + "</pre>";
      else
        if (visibilityOld == "hidden")
          innerHTML = errorMsg;
        else 
          return;
    
      this.canvasMenu.style.left = evt.pageX + "px";
      this.canvasMenu.style.top  = evt.pageY + "px";
      this.canvasMenu.innerHTML = innerHTML;
      this.canvasMenu.style.visibility = "visible";  
    }
  };
  
  
  plot.canvas = document.getElementById(plotName);  
  plot.resize ();
  
  plot.ctx = plot.canvas.getContext("2d");
  plot.ctx.fillStyle = fillStyle;  // stnd2color (genomeSize_stnd);
  plot.ctx.strokeStyle = "black";

  
  // PAR
  plot.dX = plot.ctx.measureText("M").width / 2; 

  // iframe ??
  const menuName = plotName + "Menu";
  document.write("<div id='" + menuName + "' style='border:1px solid black;background-color:#EEEEEE;visibility:hidden;position:absolute;padding-left:10px;padding-right:10px'></div>");
  plot.canvasMenu = document.getElementById(menuName);  
  
  return plot;
}



var readonly = function ()  { return true; }



var text2file = function (text, filename)
{
  const file = new Blob([text], {type: 'text/plain'});
  if (window.navigator.msSaveOrOpenBlob) // IE10+
    window.navigator.msSaveOrOpenBlob (file, filename);
  else 
  { // Others
    const a = document.createElement("a");
    const url = URL.createObjectURL(file);
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(function() {
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);  
    }, 0); 
  }
}



var drawImage = function (ctx, image, x, y, width, height, name)
// Return: success
{
  image.onload = function (e) { ctx.drawImage (image, x, y, width, height); }
  try {
                                ctx.drawImage (image, x, y, width, height);                  
      }
    catch (e) 
    { 
      alert ("Cannot draw " + name); 
      return false;
    }
  return true;
}



// https://pixabay.com/sound-effects/search/sound/
/*
var beep = function ()
{
  const beep = new Audio('sound/interface-124464.mp3'); 
  beep.play();
}
*/



var audioCtx = new AudioContext();


var beep = function (freq = 400 /*520*/, duration = 200, vol = 3) 
{
  const oscillator = audioCtx.createOscillator();
  const gain = audioCtx.createGain();
  oscillator.connect(gain);
  oscillator.frequency.value = freq;
  oscillator.type = "square";
  gain.connect(audioCtx.destination);
  gain.gain.value = vol * 0.01;
  oscillator.start(audioCtx.currentTime);
  oscillator.stop(audioCtx.currentTime + duration * 0.001);
}



var getMinutes = function ()
{
  const d = new Date ();
  return d. getHours () * 60 + d. getMinutes ();
}



/* 
  Safari restrictions:
     default parameters do not work
     variable:
       bit
*/


