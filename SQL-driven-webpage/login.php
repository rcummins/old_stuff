<?php

/***********************************************************************
* login.php
*
* Renata Cummins
* Final Project
*
* Has a box that allows previous users to log in to access the site, 
* and also has a box that allows new users to link to the registration 
* page. Also explains a little about the purpose of the website. 
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
              <h2>For members:</h2>
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
            
            <!-- Info about the site -->
            Welcome to Harvard XC Tracker. <br />
            <br />
            The goal of this website is to provide members of the Harvard
            Women's Cross Country Team with a better way to make comparisions
            between the different races they run than is provided by Harvard's
            sports website, 
            <a href="http://www.gocrimson.com">gocrimson.com</a>. 
            There, cross country meet results are organized by meet, which 
            allows runners to see how they stacked up against their competitors
            on that given day, but makes any other type of comparison 
            difficult. As a runner, I like to be able to compare my 
            performances throughout my career, which used to involve going
            through each separate page of results for every meet I ran in. 
            This website facilitates that process. 
          </div>
        </td>
        
      </tr>
      
      <tr>
      
        <td width="300px">
          <!-- Gray box with link to registration page -->
          <div class="register">
            <p>
              <br />
                Or, members of the Harvard cross country team <br />
                are elegible to 
                <a href="register.php">register</a> for an account.<br />
              <br />
            </p>
          </div>
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
