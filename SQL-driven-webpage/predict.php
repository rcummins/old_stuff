<?php

/***********************************************************************
* predict.php
*
* Renata Cummins
* Final Project
*
* Displays the prediction for the user's time at Heps based on their
* time at the meet they chose on the previous page. Also has a box 
* that allos the user to make another prediction, and has some text and
* a graph that explains a little about how the prediction was made.  
**********************************************************************/

    // require common code
    require_once("includes/common.php");
    
?>

<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                  
<html xmlns="http://www.w3.org/1999/xhtml">
                  
  <head>
    <title>Harvard XC Tracker: Predict</title>
    <link rel="stylesheet" type="text/css" href="css/mystyle.css"  />
  </head>
                              
  <body>
    <div align="center">
      <!-- Banner image at top of page -->
      <a href="index.php"><img alt="Harvard XC Tracker" border="0" src="images/banner.jpg" /></a>
    </div>
    
    <table width="100%">
      <tr>

        <td width="240px">
          <div class="register">
            
            <!-- box that allows user to make another Heps prediction -->
            <p>
              Make another Prediction!<br />
            </p>
            <form action="predict.php" method="post">
              <table border="0">
                <tr>
                  <td>Race you want to use to predict:</td>
                </tr>
                <tr>
                  <td><select name="race">
                      <option value=""></option>
                      
                      <? // prepare sql to get races run by the user
                         $query = sprintf("SELECT race_date_num FROM Results
                                        WHERE name='%s'",$_SESSION["name"]);
                         
                         // execute query
                         $result = mysql_query($query);
                         
                         // repeat for each of the rows in the table
                         while ($row = mysql_fetch_array($result))
                         {
                            // find the race_date_num
                            $n = $row["race_date_num"];
                            
                            // prepare sql to get meet name
                            $query2 = sprintf("SELECT meet,date FROM Courses
                                      WHERE race_date_num=%d",$n);
                            
                            // execute query
                            $result2 = mysql_query($query2);
                            
                            // fetch the correct array from the Courses table
                            $row2 = mysql_fetch_array($result2);
                            
                            // find the meet name and date
                            $meet = $row2["meet"];
                            $textm = sprintf("%s",$meet);
                            $date = $row2["date"];
                            $textd = sprintf("%s",$date);
                            
                      ?>
                      
                            <!-- make an option with the meet's name & date -->
                            <option value="<? print($row["race_date_num"]); ?>">
                            <? echo $textm ?> <? echo $textd ?></option>
                      <? } ?>
                    </select>
                  </td>
                </tr>
              </table>
              <p>
                <br />
                  <input type="submit" value="Predict" />
                <br />
              </p>
            </form>
          </div>
          <br />
          <div class="register">
            <p>
                <a href="index.php">Back to your Career Report</a><br />
                <br />
                <a href="logout.php">Log out</a>
              <br />
            </p>
          </div>
        </td>
        
        <td>
          <div class="report">
                   
                <? // prepare sql to match up race_date_num with meet info
                   $query3 = sprintf("SELECT cid, meet, date 
                                      FROM Courses WHERE race_date_num=%d"
                                      ,$_POST["race"]);
                   
                   // execute query
                   $result3 = mysql_query($query3);
                   
                   // grab row
                   $row3 = mysql_fetch_array($result3);
                   
                   // prepare sql to get time run for that meet
                   $query4 = sprintf("SELECT mins, secs FROM Results
                                      WHERE race_date_num=%d
                                      AND name='%s'"
                                      ,$_POST["race"],$_SESSION["name"]);
                   
                   // execute query
                   $result4 = mysql_query($query4);
                   
                   // grab row
                   $row4 = mysql_fetch_array($result4);
                   
                   // convert time into seconds
                   $time = $row4["mins"] * 60 + $row4["secs"];
                   
                   // prepare sql to grab appropriate slope and intercept
                   $query5 = sprintf("SELECT slope, intercept FROM Predict
                                      WHERE x_cid=%d AND x_meet='%s'"
                                      ,$row3["cid"],$row3["meet"]);

                   // execute query
                   $result5 = mysql_query($query5);

                   // grab returned row
                   $row5 = mysql_fetch_array($result5);

                   // convert time to equivalent Heps Time
                   $equiv = $time * $row5["slope"] + $row5["intercept"];
                   $secs = $equiv % 60;
                   $mins = ($equiv - $secs) / 60;

                ?> 
                   
            <!-- title: [Name]'s predicted time for Heps -->
            <h3><? print($_SESSION["name"]) ?>'s predicted time for Heps</h3>
              <br />
              <br />
            
            <!-- Since you ran [time] at [meet] on [date], -->
            <!-- you will probably run [time] at Heps! -->
            <p>
              Since you ran <? 
               $text = sprintf("%1.0f",$row4["mins"]);
               echo $text ?>:<?
               $text2 = sprintf("%02d",$row4["secs"]);
               echo $text2; ?> at <? print($row3["meet"]); ?> on <?
               print($row3["date"]); ?>,<br />
              You will probably run<br />
               <br />
               <div class="time"> <? 
               $text = sprintf("%1.0f",$mins);
               echo $text ?>:<?
               $text2 = sprintf("%02d",$secs);
               echo $text2; ?><br />
               </div><br /> at Heps!<br />
            </p>    
          </div>
        </td>
      </tr>
      <tr>
        <td></td>
        <td>
          <div class="report">
            <h4>About the prediction:</h4>
            <p style="font-size: 12px; font-weight: normal;">I used the 
                statistics technique of linear regression to find
                an approximate relationship between runners' times for the two 
                different meets. <br />
                Graphically, this can be illustrated by making a scatter
                plot, with the times at one meet on one axis and the times
                at the other meet on the other axis. The points represent
                runners who ran both meets in the same year. The equation
                of the line that best fits the scatter plot is the equation
                that relates the times for the two meets. <br />
                <br />
                Here is an example graph showing the plot for NCAA Regionals 
                at the Van Cortlandt Park on the x-axis and Heps on the
                y-axis: <br />
              </p>
              
              <!-- pic of one of the graphs of the data and best-fit line -->
              <img src="images/graph.jpg" alt="graph" border="0" />
          </div>
        </td>
      </tr>
    </table>
  
  </body>

</html>
