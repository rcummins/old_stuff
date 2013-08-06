<?php

/***********************************************************************
* login2.php
*
* Renata Cummins
* Final Project
*
* Checks whether the user's entered name and password are stored in 
* the Users table in the database and logs the user in if they are. 
**********************************************************************/

    // require common code
    require_once("includes/common.php"); 

    // escape username and password for safety
    $name = mysql_real_escape_string($_POST["name"]);
    $password = mysql_real_escape_string($_POST["password"]);

    // prepare SQL
    $sql = sprintf("SELECT name FROM Users WHERE name='%s' AND password='%s'",
                   $name, $password);

    // execute query
    $result = mysql_query($sql);

    // if we found a row, remember user and redirect to index
    if (mysql_num_rows($result) == 1)
    {
        // cache name in session
        $_SESSION["name"] = $name;

        // redirect to portfolio
        redirect("index.php");
    }

    // else report error
    else
        apologize("Invalid username and/or password.");

?>
