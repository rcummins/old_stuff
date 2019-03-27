<?php

/***********************************************************************
* register.php
*
* Renata Cummins
* Final Project
*
* A form that allows new users to register for a username and password
* in order to access the site. 
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
    <title>Harvard XC Tracker: Register</title>
    <link rel="stylesheet" type="text/css" href="css/mystyle.css"  />
  </head>
                              
  <body>
    <div align="center">
      <!-- Banner image at top of page -->
      <a href="index.php"><img alt="Harvard XC Tracker" border="0" src="images/banner.jpg" /></a>
    </div>

    <br />

    <table align="center" width="600px">
      <tr>
        
        <!-- box with a form to fill in to register -->
        <td>    
          <div class="register" align="center">
              <br />
              <h3>For members of the Harvard team:</h3>
              <br />
            <form action="register2.php" method="post">
              <table align="center" border="0">
                <tr>
                  
                  <!-- Name used to search for the user's meet results -->
                  <td class="field">Your name, as it usually appears in meet results: </td>
                  <td><input name="name" type="text" size="30" /></td>
                </tr>
                <tr><td> </td></tr>
                <tr>
                  
                  <!-- password used the next time user wants to log in -->
                  <td class="field">Password (max 16 characters): </td>
                  <td><input name="password" type="password" size="30" /></td>
                </tr>
                <tr><td> </td></tr>
                <tr>
                  
                  <!-- checks for typos -->
                  <td class="field">Type your password again: </td>
                  <td><input name="password2" type="password" size="30" /></td>
                </tr>
                <tr><td> </td></tr>
                <tr>
                  
                  <!-- I want their email address to discourage stalkers -->
                  <!-- from using the website and to be able to email -->
                  <!-- anyone who I suspect is a stalker -->
                  <td class="field">Your Harvard email address: </td>
                  <td><input name="email" type="text" size="30" /></td>
                </tr>
              </table>
              <p>
                <br />
                  <input type="submit" value="Register" /><br />
                <br />
                  Oops, I already have an account. 
                  <a href="login.php">Login</a>.
                <br />
              </p>
            </form>
          </div>  
        </td>
      </tr>
    </table>

  </body>

</html>
