<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Query Genotype Data: Query Error</title>
</head>
<body>
<h1>Query Error</h1>
<p>Whoops, Looks like something went wrong when running your query.<br>
    Please check that you have submitted stocks for the right dataset. You can check under "View Sources" <a href="http://snpversity.maizegdb.org/Diversity/home.php">here</a>.</p>
<p> If you are manually selecting positions, please also verify that your selected end position is greater than the start site of the selected chromosome in the dataset.</p>
<?php
$uId = $_GET["query"];
echo "<p>" . "This error has been logged. Your query ID is: <strong>" . $uId . "</strong></p>";
?>
</body>
</html>
