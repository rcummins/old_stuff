<?php

/***********************************************************************
* course.php
*
* Renata Cummins
* Final Project
*
* Displays all of the user's times on the spcified course, and then below
* that, all Harvard runners' fastest times (PRs) on the specified
* course. 
**********************************************************************/

    // require common code
    require_once("includes/common.php");
    
?>

<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                  
<html xmlns="http://www.w3.org/1999/xhtml">
                  
  <head>
    <title>Harvard XC Tracker: Course PRs</title>
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
            
            <!-- links to other Harvard Course PRs -->
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
                <a href="index.php">Back to your Career Report</a><br />
              <br />
                <a href="logout.php">Log out</a>
              <br />
            </p>
          </div>
        </td>

        <td>
          <div class="report">
            <h3><? 
              
              // prepare sql to get course name and length from Courses table
              $query3 = sprintf("SELECT course_name, length FROM Courses
                                WHERE cid=%d",$_GET["course"]);
                                
              // execute query
              $result3 = mysql_query($query3);
              
              // grab row returned
              $row3 = mysql_fetch_array($result3);
              
              // convert length to km
              $length = $row3["length"]/1000;
              
              // print course name
              print($row3["course_name"]); ?> <?
              
              // print length in km
              echo $length; ?>k Course</h3><br />
            <br />
            <br />
            
            <!-- title: yours races on the [course name and length] -->
            <h4>Your races on the <? 
                print($row3["course_name"]); ?> <?
                echo $length; ?>k Course</h4>
            <br />
            <table border="1" cellpadding="5px" rules="rows">
            
              <!-- print table headers for table with user's times -->
              <tr>
                <th class="report">Date</th>
                <th class="report">Meet</th>
                <th class="report">Time</th>
                <th class="report">Pace</th>
              </tr>
             
              <? // prepare sql to get time and race_date_num from Results
                 $query = sprintf("SELECT mins, secs, race_date_num 
                                   FROM Results WHERE name='%s' AND cid=%d"
                                   ,$_SESSION["name"],$_GET["course"]);
                
                 // execute query
                 $result = mysql_query($query);
                 
                 // grab every row returned
                 while ($row = mysql_fetch_array($result))
                 {
                   // prepare sql to get meet info from Courses table
                   $query2 = sprintf("SELECT date, event_name, length
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
                   
                   <!-- print table rows in table with user's times -->
                   <tr>
                     
                     <!-- Date -->
                     <td class="report">
                       <a href="meet_results.php?race=<?
                         print($row["race_date_num"]); ?>">
                         <? print($row2["date"]); ?></a></td>
                     
                     <!-- Meet -->
                     <td class="report"><? print($row2["event_name"]); ?></td>
                     
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
      <tr>
        <td></td>
        <td>
          <div class="report">
            
            <!-- title: Harvard PRs on the [course name and length] -->
            <h4>Harvard PRs on the <? 
                print($row3["course_name"]); ?> <?
                echo $length; ?>k Course (from results since 1998)</h4>
            <br />
            <table border="1" cellpadding="5px" rules="rows" >
             
              <div align="center">
              <!-- print table headers for table with Harvard PRs -->
              <tr>
                <th class="report">Runner</th>
                <th class="report">Time</th>
                <th class="report">Pace</th>
                <th class="report">Date</th>
                <th class="report">Meet</th>
              </tr>
             
              <? // prepare sql to get names of people who ran at the course
                 $query4 = sprintf("SELECT DISTINCT name FROM ResultsTime
                                   WHERE school='Harvard' AND cid=%d"
                                   ,$_GET["course"]);
                
                 // execute query
                 $result4 = mysql_query($query4);
                 
                 // grab every row returned
                 while ($row4 = mysql_fetch_array($result4))
                 {
                   // prepare sql to get time and race_date_num from Results
                   $query = sprintf("SELECT mins, secs, race_date_num 
                                     FROM ResultsTime WHERE name='%s' 
                                     AND cid=%d LIMIT 1"
                                     ,$row4["name"],$_GET["course"]);
                                     
                   // execute query
                   $result = mysql_query($query);
                   
                   // grab row returned
                   $row = mysql_fetch_array($result);
                   
                   // prepare sql to get meet info from Courses table
                   $query2 = sprintf("SELECT date, event_name
                                      FROM Courses WHERE race_date_num=%d"
                                      ,$row["race_date_num"]);
                   
                   // execute query
                   $result2 = mysql_query($query2);
                   
                   // grab the row that is returned
                   $row2 = mysql_fetch_array($result2); 
                   
                   // calculate the pace
                   $seconds = $row["mins"]*60 + $row["secs"];
                   $miles = $row3["length"]/1609;
                   $pace_in_secs = $seconds/$miles;
                   $pace_secs = $pace_in_secs%60;
                   $pace_mins = ($pace_in_secs - $pace_secs)/60;
                   
                   // calculate the Equivalent VCP 5k Time
                   
                   // TO DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   
                   $equiv_mins = $row["mins"];
                   $equiv_secs = $row["secs"];
                   
                   ?>
                   
                   <!-- print table rows for table with Harvard PRs -->
                   <tr>
                     
                     <!-- Runner -->
                     <td class="report">
                       <a href="indiv_career.php?name=<? 
                         print($row4["name"]); ?>">
                         <? print($row4["name"]); ?></a></td>
                     
                     <!-- Time -->
                     <td class="report"><? print($row["mins"]); ?>:<?
                            $text1 = sprintf("%02d",$row["secs"]);
                            echo $text1; ?></td>
                     
                     <!-- Pace -->
                     <td class="report"><? $text2 = sprintf("%1.0f",$pace_mins); 
                            echo $text2; ?>:<?
                            $text3 = sprintf("%02d",$pace_secs);
                            echo $text3; ?></td>
                     
                     <!-- Date -->
                     <td class="report">
                       <a href="meet_results.php?race=<?
                         print($row["race_date_num"]); ?>">
                         <? print($row2["date"]); ?></a></td>
                     
                     <!-- Meet -->
                     <td class="report"><? print($row2["event_name"]); ?></td>
                     
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
