<?php

/***********************************************************************
* index.php
*
* Renata Cummins
* Final Project
*
* Has a bar off the the side with the Heps prediction tool, and links 
* to Harvard course PRs and to log out. Next to that, there is a graph 
* showing the user's race pace throughout her career. This page also 
* displays a table with all the races the user has run during her 
* career. If the user is not logged in, redirects to the login page.
**********************************************************************/

    // require common code
    require_once("includes/common.php");
    
?>

<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                  
<html xmlns="http://www.w3.org/1999/xhtml">
                  
  <head>
    <title>Harvard XC Tracker: Home</title>
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
            
            <!-- Heps prediction tool -->
            <p>
              Predict your Heps Time!<br />
            </p>
            <form action="predict.php" method="post">
              <table width="100%" border="0">
                <tr>
                  <td>Race to use to predict:</td>
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
                            <option value="<? print($row["race_date_num"]) ?>">
                            <? echo $textm ?> <? echo $textd ?></option>
                      <? } ?>
                  </td>
                </tr>
              </table>
              <p>
                  <input type="submit" value="Predict" />
                <br />
              </p>
            </form>
          </div>
          <br />
          <div class="register">
            
            <!-- links to the Harvard Course PRs -->
            <p>
                Harvard Course PRs:<br />
                  <div class="links">
                    <a href="course.php?course=15">
                      Van Cortlandt Park 5k</a><br />
                    <a href="course.php?course=16">
                      Van Cortlandt Park 6k</a><br />
                    <a href="course.php?course=25">
                      Franklin Park 5k</a><br />
                    <a href="course.php?course=26">
                      Franklin Park 6k</a><br />
                    <a href="course.php?course=36">
                      Terre Haute 6k</a><br />
                    <a href="course.php?course=45">
                      Albany 5k</a><br />
                    <a href="course.php?course=66">
                      Lehigh 6k</a><br />
                    <a href="course.php?course=55">
                      Yale 5k</a><br />
                    <a href="course.php?course=75">
                      Princeton 5.1k</a>
                  </div>
              <br />
            </p>
          </div>
          <br />
          <div class="register">
            <p>
                <a href="logout.php">Log out</a>
              <br />
            </p>
          </div>
        </td>
        
        <td>
          <div class="report">
            
            <!-- title: Career Report for [user's name] -->
            <h3>Career Report for <? print($_SESSION["name"]) ?></h3><br />
            <br />

             <? // prepare sql to get races run, times from Results table
                $query = sprintf("SELECT mins, secs, race_date_num
                                  FROM Results WHERE name='%s'"
                                  ,$_SESSION["name"]);
                 
                // execute query
                $result = mysql_query($query);
                
                // grab every row returned
                while ($row = mysql_fetch_array($result))
                {
                  // prepare sql to get meet info from Courses table
                  $query2 = sprintf("SELECT length 
                                     FROM Courses WHERE race_date_num=%d"
                                     ,$row["race_date_num"]);
                  
                  // execute query
                  $result2 = mysql_query($query2);
                  
                  // grab the row that is returned
                  $row2 = mysql_fetch_array($result2); 
                  
                  // calculate the pace
                  $seconds = $row["mins"]*60 + $row["secs"];
                  $miles = $row2["length"]/1609;
                  $pace_in_secs = $seconds/$miles;
                  $pace = $pace_in_secs/60;
                  
                  // add pace to array
                  $array[] = $pace;
                  
                } 

                // graph showing mile pace over career ?>
                <img src="http://chart.apis.google.com/chart?<?
                
                // chart type = line graph 
                ?>cht=lxy<?
                
                // chart size in pixels
                ?>&chs=500x300<?
                
                // space date points evenly along x-axis
                ?>&chd=t:-1|<?
                
                // print each value in the pace array, comma separated
                for($i=0; $i < count($array) - 1; $i++)
                {
                    $pace_text = sprintf("%1.2f",$array[$i]);
                    print($pace_text);
                    print ",";
                }
                $n = count($array) -1;
                $pace_text = sprintf("%1.2f",$array[$n]);
                print($pace_text);
                
                // find next .25 number below minimum pace
                $minimum = min($array);
                $minimum0 = $minimum * 1000000;
                $min0 = $minimum0 - ($minimum0%250000);
                $min = $min0 / 1000000;
                $min_text = sprintf("%1.2f",$min);
                
                // find next .25 number above maximum pace
                $maximum = max($array);
                $maximum0 = $maximum * 1000000;
                $max0 = $maximum0 - ($maximum0%250000);
                $max25 = $max0 / 1000000;
                $max = $max25 + 0.25;
                $max_text = sprintf("%1.2f",$max);
                    
                // date range
                ?>&chds=<? print($min_text); ?>,<? print($max_text);
                
                // chart title = Your Mile Pace Over Your Career
                ?>&chtt=Your+Mile+Pace+Over+Your+Career<?
                
                // chart title style = black and size 14 font
                ?>&chts=000000,14<?
                
                // label y-axis with Minutes values
                ?>&chxt=y&chxl=0:<?
                for ($i = $min; $i <= $max + 0.1; $i += 0.25)
                {
                    print "|";
                    $mins_text = sprintf("%1d",$i);
                    print ($mins_text);
                    print "+";
                    $time100 = $i * 100;
                    $mins100 = $time100 - ($time100%100);
                    $secs100 = $time100 - $mins100;
                    $secs = $secs100 / 100; 
                    $secs60 = $secs * 60;
                    $secs_text = sprintf("%02d",$secs60);
                    print ($secs_text);
                    print "+";
                }
                
                // Put the labels in the right place
                ?>&chxp=0<?
                
                // print the percentage point of position 
                $count = ($max - $min) / 0.25;
                $increment = 100 / $count;
                $percent = 0;
                for ($i=0; $i <= $count + 0.1; $i++)
                {
                    print ",";
                    $percent_text = sprintf("%03d",$percent);
                    print($percent_text);
                    $percent = $percent + $increment;
                }
                print "|";
            
                // line color = crimson
                ?>&chco=770000<?
                
                // background color = gray 
                ?>&chf=bg,s,EEEEEE<?
                
                // line style = solid
                ?>&chls=2,1,0<?
                
                // margins around the chart
                ?>&chma=40,20,40,20<?
                
                // each data point is a black diamond
                ?>&chm=d,000000,0,-1,8,1" />

            <br />
          </div>
        </td>
      </tr>
      <tr>
        <td><br />
        </td>
        <td>
          <div class="report">
            <table border="1" cellpadding="5px" rules="rows" >
             
              <div align="center">
              <!-- print table headers for table with user's career history -->
              <tr>
                <th class="report">Date</th>
                <th class="report">Meet</th>
                <th class="report">Course</th>
                <th class="report">Distance</th>
                <th class="report">Time</th>
                <th class="report">Pace</th>
              </tr>
             
              <? // prepare sql to get races run, times from Results table
                 $query = sprintf("SELECT mins, secs, race_date_num
                                   FROM Results WHERE name='%s'"
                                   ,$_SESSION["name"]);
                
                 // execute query
                 $result = mysql_query($query);
                 
                 // grab every row returned
                 while ($row = mysql_fetch_array($result))
                 {
                   // prepare sql to get meet info from Courses table
                   $query2 = sprintf("SELECT date, course_name, cid, length, 
                                      event_name, meet
                                      FROM Courses WHERE race_date_num=%d"
                                      ,$row["race_date_num"]);
                   
                   // execute query
                   $result2 = mysql_query($query2);
                   
                   // grab the row that is returned
                   $row2 = mysql_fetch_array($result2); 
                   
                   // calculate the pace
                   $seconds = $row["mins"]*60 + $row["secs"];
                   $miles = $row2["length"]/1609;
                   $pace_in_secs = $seconds/$miles;
                   $pace_secs = $pace_in_secs%60;
                   $pace_mins = ($pace_in_secs - $pace_secs)/60;
                   
                   ?>
                   
                   <!-- print rows for table with user's career history -->
                   <tr>
                     
                     <!-- Date -->
                     <td class="report">
                       <a href="meet_results.php?race=<?
                         print($row["race_date_num"]); ?>">
                         <? print($row2["date"]); ?></a></td>
                     
                     <!-- Meet -->
                     <td class="report"><? print($row2["event_name"]); ?></td>
                     
                     <!-- Course -->
                     <td class="report"><? print($row2["course_name"]); ?></td>
                     
                     <!-- Distance -->
                     <td class="report"><? $distance = $row2["length"]/1000; 
                            echo $distance ?>k</td>
                     
                     <!-- Time -->
                     <td class="report"><? print($row["mins"]); ?>:<?
                            $text1 = sprintf("%02d",$row["secs"]);
                            echo $text1; ?></td>
                     
                     <!-- Pace -->
                     <td class="report"><? 
                            $text2 = sprintf("%1.0f",$pace_mins); 
                            echo $text2; ?>:<?
                            $text3 = sprintf("%02d",$pace_secs);
                            echo $text3; ?></td>
                     
                   </tr>
              <? } ?>
            
              </div>
            </table> 
          </div>
        </td>
      </tr>
    </table>
  
  </body>

</html>
