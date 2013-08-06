<?php

/***********************************************************************
* meet_results.php
*
* Renata Cummins
* Final Project
*
* Displays the results of the meet whose date was clicked in a previous 
* page. Displays the Harvard results at the top, and if I have the full
* results uploaded in my database, the full results under that. 
**********************************************************************/

    // require common code
    require_once("includes/common.php");
    
?>

<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
                  
<html xmlns="http://www.w3.org/1999/xhtml">
                  
  <head>
    <title>Harvard XC Tracker: Meet Results</title>
    <link rel="stylesheet" type="text/css" href="css/mystyle.css"  />
  </head>
                              
  <body>
    <div align="center">
      <!-- Banner image at top of page -->
      <a href="index.php"><img alt="Harvard XC Tracker" border="0" src="images/banner.jpg" /></a>
    </div>
    
    <br />
    
    <table width="100%">
      <tr>

        <td width="240px">
          <div class="register">
            
            <!-- links to pages with Harvard Course PRs -->
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
                <a href="logout.php">Log out</a><br />
            </p>
          </div>
        </td>
        
        <td>
          <div class="report">
            
            <!-- title: Results From [meet] at [course] on [date], [length] -->
            <h3>Results from <? 
                
                // get $race from url
                $race = $_GET["race"];
                
                // sql to get meet, course, date, length from Courses table
                $query2 = sprintf("SELECT meet,course_name,date,length 
                                  FROM Courses WHERE race_date_num=%d",$race);
                                  
                // execute query
                $result2 = mysql_query($query2);
                
                // grab row
                $row2 = mysql_fetch_array($result2);
                
                // print meet name
                print($row2["meet"]); ?> at <?
                
                // print course name
                print($row2["course_name"]); ?> on <?
                
                // print meet date
                print($row2["date"]); ?> (<?
                
                // convert length to km
                $length = $row2["length"]/1000;
                echo $length; ?>k)</h3><br />
            <br />
            <h4>Harvard Results</h4>
            <br />
            <table border="1" cellpadding="5px" rules="rows" >
             
              <!-- print table headers for Harvard results from the meet -->
              <tr>
                <th class="report">Place</th>
                <th class="report">Runner</th>
                <th class="report">Time</th>
                <th class="report">Pace</th>
              </tr>
             
              <? // prepare sql to get races run, times from Results table
                 $query = sprintf("SELECT place, name, mins, secs 
                                   FROM Results 
                                   WHERE race_date_num=%d AND school='Harvard'"
                                   ,$race);
                
                 // execute query
                 $result = mysql_query($query);
                 
                 // grab every row returned
                 while ($row = mysql_fetch_array($result))
                 {
                   // calculate the pace
                   $seconds = $row["mins"]*60 + $row["secs"];
                   $miles = $row2["length"]/1609;
                   $pace_in_secs = $seconds/$miles;
                   $pace_secs = $pace_in_secs%60;
                   $pace_mins = ($pace_in_secs - $pace_secs)/60;
                   
                   ?>
                   
                   <!-- print table rows for table with Harvard results -->
                   <tr>
                     
                     <!-- Place -->
                     <td class="report"><? print($row["place"]); ?></td>
                     
                     <!-- Runner -->
                     <td class="report">
                       <a href="indiv_career.php?name=<? 
                         print($row["name"]); ?>">
                         <? print($row["name"]); ?></a></td>
                     
                     <!-- Time -->
                     <td class="report"><? print($row["mins"]); ?>:<?
                            $text1 = sprintf("%02d",$row["secs"]);
                            echo $text1; ?></td>
                     
                     <!-- Pace -->
                     <td class="report"><? $text2 = sprintf("%1.0f",$pace_mins); 
                            echo $text2; ?>:<?
                            $text3 = sprintf("%02d",$pace_secs);
                            echo $text3; ?></td>
                     
                   </tr>
              <? } ?>
            
            </table> 
          </div>
        </td>
      </tr>
      <tr>
        <td><br /></td>
        <td>
          <div class="report">
          <h4>Full Results</h4>
          <br />
          <table border="1" cellpadding="5px" rules="groups" >
        
          <!-- print table headers for table with full meet results -->
          <tr>
            <th class="report">Place</th>
            <th class="report">Runner</th>
            <th class="report">School</th>
            <th class="report">Time</th>
            <th class="report">Pace</th>
          </tr>
        
              <? // prepare sql to get races run, times from Results table
              $query3 = sprintf("SELECT place, name, school, mins, secs
                                FROM Results WHERE race_date_num=%d"
                                ,$race);
            
              // execute query
              $result3 = mysql_query($query3);
        
              // grab every row returned
              while ($row3 = mysql_fetch_array($result3))
              {
                // calculate the pace
                $seconds = $row3["mins"]*60 + $row3["secs"];
                $miles = $row2["length"]/1609;
                $pace_in_secs = $seconds/$miles;
                $pace_secs = $pace_in_secs%60;
                $pace_mins = ($pace_in_secs - $pace_secs)/60;
                
                ?>
            
                <!-- print table rows for table with full meet results -->
                <tr>
            
                    <!-- Place -->
                    <td class="report"><? print($row3["place"]); ?></td>
            
                    <!-- Runner -->
                    <td class="report"><? print($row3["name"]); ?></td>
            
                    <!-- School -->
                    <td class="report"><? print($row3["school"]); ?></td>
                
                    <!-- Time -->
                    <td class="report"><? print($row3["mins"]); ?>:<?
                                      $text1 = sprintf("%02d",$row3["secs"]);
                                      echo $text1; ?></td>
            
                    <!-- Pace -->
                    <td class="report"><? $text2 = sprintf("%1.0f",$pace_mins);
                                          echo $text2; ?>:<?
                                          $text3 = sprintf("%02d",$pace_secs);
                                          echo $text3; ?></td>
            
                </tr>
            <? } ?>
            
            </table>
          </div>
        </td>
      </tr>
        
    </table>
  
  </body>

</html>
