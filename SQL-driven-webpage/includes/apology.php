<?php

/***********************************************************************
* apology.php
*
* Renata Cummins
* Final Project
*
* Apologizes for any error, providing text explaining the error. 
**********************************************************************/

?>
                                       
<!DOCTYPE html 
     PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
     "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">

  <head>
    <link href="css/mystyle.css" rel="stylesheet" type="text/css" />
    <title>Harvard XC Tracker: Sorry!</title>
  </head>

  <body>

    <div align="center">
        <!-- Banner image at top of page -->
        <a href="index.php">
          <img alt="Harvard XC Tracker" border="0" src="images/banner.jpg" >
        </a>
    </div>

    <br />
    <br />

    <div align="center">
      <h3>Sorry!</h3><br />
    </div>

    <!-- explains the error -->
    <div align="center">
      <?php echo $message; ?>
    </div>

    <!-- link to the previous page -->
    <div align="center" style="margin: 20px;">
      <a href="javascript:history.go(-1);">Back</a>
    </div>

  </body>

</html>
