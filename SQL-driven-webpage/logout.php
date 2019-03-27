<?php

/***********************************************************************
* logout.php
*
* Renata Cummins
* Final Project
*
* Has a box that allows previous users to log in again to access the 
* site, and provides some useful links for navigation away from the site. 
**********************************************************************/

    // require common code
    require_once("includes/common.php");
        
    // log current user out, if any
    $_SESSION = array();
    if (isset($_COOKIE[session_name()]))
    setcookie(session_name(), "", time() - 42000, "/");
    session_destroy();
                                
?>

<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                  
<html xmlns="http://www.w3.org/1999/xhtml">
                  
  <head>
    <title>Harvard XC Tracker: Login</title>
    <link rel="stylesheet" type="text/css" href="css/mystyle.css"  />
  </head>
                              
  <body>
    <div align="center">
      <!-- Banner image at top of page -->
      <a href="index.php"><img alt="Harvard XC Tracker" border="0" src="images/banner.jpg" /></a>
    </div>
    
    <br />
    
    <table width="80%" align="center">
      <tr>

        <!-- Crimson login box for members -->
        <td width="300px">
          <div class="login">
              <br />
              <h2>Log in again:</h2>
              <br />
            <form action="login2.php" method="post">
              <table align="center" border="0">
                <tr>
                  <td class="field">Name: </td>
                  <td><input name="name" type="text" size="30" /></td>
                </tr>
                <tr><td> </td></tr>
                <tr><td> </td></tr>
                <tr><td> </td></tr>
                <tr>
                  <td class="field">Password: </td>
                  <td><input name="password" type="password" size="30" /></td>
                </tr>
              </table>
              <p>
                <br />
                  <input type="submit" value="Log In" /><br />
                <br />
              </p>
            </form>
          </div>  
        </td>
        
        <td>
          <div  class="info">
            You are logged out. <br />
            <br />
            Some useful links: <br />
            <a href="http://www.gocrimson.com" target="_blank">gocrimson.com</a><br /> 
            <a href="http://www.flotrack.org/auth/login" target="_blank">flotrackr</a><br /> 
            <a href="http://my.harvard.edu" target="_blank">my.harvard.edu</a><br />
            <br />
            <br />
            <br />
          </div>
        </td>
        
      </tr>
      
      <tr>
      
        <td width="50px">
        </td>
        
        <td>
          <div class="infosmall">
            Made by Renata Cummins for her final project in Computer 
            Science 50.<br />
            Email her at rcummins @ fas.harvard.edu with comments or
            suggestions. 
          </div>
        </td>
      </tr>
    </table>

  </body>

</html>
