<?php

/***********************************************************************
* register2.php
*
* Renata Cummins
* Final Project
*
* Checks to make sure the user entered everything correctly in the form
* and apologizes if there are any errors. If everything goes well, logs 
* the user in and redirects them to the index page. 
**********************************************************************/

    // require common code
    require_once("includes/common.php"); 

    // escape entered values for safety
    $name = mysql_real_escape_string($_POST["name"]);
    $email = mysql_real_escape_string($_POST["email"]);
    $password = mysql_real_escape_string($_POST["password"]);
    $password2 = mysql_real_escape_string($_POST["password2"]);
    
    // check if name is blank
    if ($_POST["name"] == NULL)
        apologize("Please enter your name.");
        
    // check if email is blank
    if ($_POST["email"] == NULL)
        apologize("Please enter your email address.");
        
    // check if password or password2 is blank
    if ($_POST["password"] == NULL || $_POST["password2"] == NULL)
        apologize("Please enter your password in both password fields.");
        
    // check if password and password2 are the same
    if ($password != $password2)
        apologize("Your two passwords don't match.");

    // prepare SQL
    $sql = sprintf("INSERT INTO Users (name, password, email) 
                    VALUES ('%s', '%s', '%s')", $name, $password, $email);

    // execute query
    $result = mysql_query($sql);

    // if INSERT failed, apologize and ask user to enter different username
    if ($result == FALSE)
        apologize("This name is already registered. <br />
        If you have forgtten your password, please email Renata at 
        rcummins @ fas.harvard.edu and she will send you an email reminding 
        you of your password. <br />
        If you suspect that someone else has tried to sign up as you, or are 
        having any other problems, please also email Renata.");

    // if INSERT was successful, log the user in
    else
    {
        // cache name in session
        $_SESSION["name"] = $name;
        
        // redirect to index.php
        redirect("index.php");
    }
?>
