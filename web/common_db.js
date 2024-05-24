// common_db.js


// Global constants
var offline = false;  // PAR
var lite = offline ? true : false /*PAR*/;
// System
var server = "PROTEUS";
var dbuser_readonly = "anyone";
var dbuser = dbuser_readonly;
var dbpassword_readonly = "allowed";
var dbpassword = dbpassword_readonly;
var database = 'uniColl';
var userName = '';




var str2sql = function (s)
{
  var s1 = '';
  if (s && s.length)
    for (var i = 0; i < s.length; i++)
      s1 += s[i] == "'" ? "''" : s[i];
  return s1;
}



var str2sqlNull = function (s)
{
  if (s == null)
    return 'null';
  return "'" + str2sql (s) + "'";
}



var sqlInjected = function (sql)
{
  if (! sql || sql.indexOf("'") == -1)
    return false;
  alert ("A string cannot contain quotes")
  return true;
}



/*
var str2html = function (s)
{
  var s1 = '';
  if (s && s.length)
    for (var i = 0; i < s.length; i++)
    {
      c = s[i];
      switch (c)
      {
        case '=': c = "%3D"; break;
        case '&': c = "%26"; break;
        case '%': c = "%25"; break;
        case '+': c = "%2B"; break;
      }
      s1 += c;
    }
  return s1;
}
*/



var select2columns = function (sql)  
// Return: [columns], [] - fail
{
  const selectS = "SELECT ";
  if (sql. substr (0, selectS. length). toUpperCase() != selectS)
    return [];

  let i = selectS. length;

  const fromS = "FROM ";

  var skipExpression = function ()
  // Update: i
  // Return: success
  {
  //alert (sql);  
    while (i < sql. length && sql. charAt (i) == ' ')
      i++;
    if (i == sql. length)
      return false;
    let parens = 0;
    let comment = false;
    // "...", [...], '...''...' ??
    while (i < sql. length)
    {
    //alert (i);  
      switch (sql. charAt (i))
      {
        case '(': parens++; break;
        case ')': if (parens) parens--; else return false; break;
        case ',': if (! parens) return true; break;
        case '*': if (sql. charAt (i - 1) == '/') comment = true; break;
        case '/': if (sql. charAt (i - 1) == '*') comment = false; break;
        default: if (! parens && ! comment && sql. substr (i, fromS. length). toUpperCase() == fromS)
                   return true;
      }      
      i++;
    }
    return ! parens && ! comment;
  }
  
  var getCol = function (exprStart, stop)
  // Return: null - fail
  // Input: sql
  {
    for (;;)
    {
      while (stop >= exprStart && sql. charAt (stop) == ' ')
        stop--;
      // Remove comment
      if (sql. substr (stop - 1, 2) != "*/")
        break;
      const pos = sql. lastIndexOf ("/*", stop - 2);
      if (pos == -1)
        return null;
      stop = pos - 1;
    }
    let start = stop - 1;
    while (   start >= exprStart 
           && ! (   sql. charAt (start) == ' '
                 || sql. charAt (start) == '.'
                )
          )
      start--;
  //alert (start);  
  //alert (stop);   
    let col = sql. substr (start + 1, stop - start);
    if (! col)
      return null;
      
    const first = col. charAt (0);
    const last  = col. charAt (col. length - 1);
    if (   first == '[' && last == ']'
        || first == '"' && last == '"'
       )
      return col. substr (1, col. length - 2); 
    return col;
  }
  

  arr = [];
  for (;;)
  {
    const exprStart = i;
    if (! skipExpression (sql, i))
      return [];
  //alert (i);  
    const col = getCol (exprStart, i - 1);
  //alert (col);
    if (! col)
      return [];
    if (i == sql. length || sql. substr (i, fromS. length). toUpperCase() == fromS)
    {
      arr. push (col);
      return arr;
    }
    if (sql. charAt (i) == ',') 
    {
      arr. push (col);
      i++;
    }
  }
}



var readonly = function ()
{
  return dbuser == dbuser_readonly;
}



var query_ = function (server, dbuser, dbpassword, database, sql, timeout_sec)
// Requires: sql: columns have distinct names
//           result set must have <= 10000/*??*/ rows 
// Returns: map: 0-based numbers as strings -> map (columns -> values)
//            numbers with very large precision are converted to strings
{
  if (offline) alert (sql);  
  
//prompt ("Enter", sql);  
  
  var rs = [];  
  $.ajax
    (
      {  
        url: 'query.cgi'  // was: query.php
      , data: { server: server
              , user: dbuser // dbuser_readonly
              , password: dbpassword  // dbpassword_readonly
              , db: database
              , sql: sql  
              , dml: 0
              , timeout: timeout_sec  // Was used in PHP
              , log: 0  
              }
      , type: "GET"  
      , dataType: 'json'
      , async: false
      , cache: false
      , timeout: timeout_sec * 1000
    
      , success: function (res)        
          {
            if (res.msg)
            {
              alert (res.msg);
              rs = res;
            }
            else
            {
            //alert (res.data);  
            //alert (encodeURIComponent(sql));
              rs = res.data;
            }
          }
        
      , error: function (jqXHR, textStatus, errorThrown) 
          { 
            alert ("State: " + jqXHR.readyState + "\nStatus: " + textStatus + "\nError: " + errorThrown + "\nCannot query the database " + database + " on server " + server + "\nSQL:\n" + sql); 
          }
      }
    );
    
  return rs;
}



var query_sybase = function (server, dbuser, dbpassword, sql)
// PHP is not implemented ??
{
  if (offline) alert (sql);  
  
//prompt ("Enter", sql);  
  
  var rs = [];  
  $.ajax
    (
      {  
        url: 'query_sybase.php'
      , data: { server: server
              , user: dbuser // dbuser_readonly
              , password: dbpassword  // dbpassword_readonly
              , sql: sql
              }
      , type: "GET"
      , dataType: 'json'
      , async: false
      , cache: false
    
      , success: function (res)        
          {
            if (res.msg)
            {
              alert (res.msg);
              rs = res;
            }
            else
              rs = res.data;
          }
        
      , error: function (jqXHR, textStatus, errorThrown) 
          { 
            alert ("Status: " + textStatus + "\nError: " + errorThrown + "\nCannot query the Sybase server " + server); 
          }
      }
    );
    
  return rs;
}



var query = function (sql)
{
  return query_ (server, dbuser, dbpassword, database, sql);
}

  
    
  
var check_credentials_ = function (server, dbuser, dbpassword, database)
{
  var rs = query_ (server, dbuser, dbpassword, database, "select 1", 3);
  if (rs.length == 1)
    return true;    
  alert ("Cannot connect to the database " + database + " on the server " + server);
  return false;
}



var check_credentials = function ()
{
  return check_credentials_ (server, dbuser, dbpassword, database);
}



var dml_ = function (server, dbuser, dbpassword, database, sql)
// Return: number of rows affected or -1 if there is an error
// Input: sql: no transaction
{
  if (offline) alert (sql);  

  var count = -1;
  $.ajax
    (
      {  
        url: 'query.cgi'  // was: query.php
      , data: { server: server
              , user: dbuser
              , password: dbpassword
              , db: database
              , sql: sql
              , dml: 1
              , log: 0
              }
      , type: "GET"
      , dataType: 'json'
      , async: false
      , cache: false
    
      , success: function (res)        
          {
            if (res.msg)
              alert (res.msg + "\n\n" + sql);
            else
              count = res.count;
          }
        
      , error: function (jqXHR, textStatus, errorThrown) 
          { 
            alert ("Status: " + textStatus + "\nError: " + errorThrown + "\nCannot change data in the database " + database + " on server " + server + "\n" + sql); 
          }
      }
    );
    
  return count;
}



var dml = function (sql)
{
  return dml_ (server, dbuser, dbpassword, database, sql);
}
  


var editArea = function (text, updateFunction)
{
  if (readonly ()) 
    return text;
//return "<form onsubmit='" + updateFunction + ";return false'><textarea rows=5 cols=130 id='f' onfocusout='alert(\"Edited text will be lost!\")'>" + q2h(text) + "</textarea><br><input type=\"submit\" value=\"Save\"></form>";
  return "<form onsubmit='" + updateFunction + ";return false'><textarea rows=5 cols=130 id='f'>" + q2h(text) + "</textarea><br><input type=\"submit\" value=\"Save\"></form>";
}



var initSystem = function ()
// Output: debug, userName
{
  // debug, userName
  if (! lite)
  {
    $.ajax
      (
        {  
          url: 'user.php',        
        //data: {}, 
          type: "GET",
          dataType: 'json',               
          async: false,
          timeout: 1000,
      
          success: function (res)        
            {
              userName = res.client;
              debug = (res.client.toLowerCase() == "ncbipc8167.be-md.ncbi.nlm.nih.gov");
            },

          error: function () 
            { alert ("Cannot get user"); }
        }
      );
    if (sqlInjected (userName))
      userName = '';
  }

//if (userName && ! debug && ! lite /*&& server == "PROTEUS"*/)
  //dml_ ("PROTEUS", "anyone", "allowed", "uniColl", "insert into WEBLOG (username, page) values ('" + userName + "', '" + location.pathname + "')");  
}



