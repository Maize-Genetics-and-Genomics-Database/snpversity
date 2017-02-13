<?php
/*Used for running actual query on server, updating DB, and returning results to send.php
 *CAUTION* Please be very careful when editing this sensitive and crucial file.
 */
require_once("db_controller.php");
require_once("query_timer.php");
$db_handle = new DBController();

$jsonObj = $_POST['data'];
$data = json_decode($jsonObj, true);

//Decode query from JSON-object
$assembly = $data["assembly"];
$chromosome = $data["chromosome"];
$positions = $data["positions"];
$startPosition = $data["startPosition"];
$endPosition = $data["endPosition"];
$taxaArray = $data["taxaArray"];
$output_format = $data["outputFormat"];
$data_set = $data["dataSet"];
$results_max = $data["resultsPerPage"];
$uId = $data["query"];

//Begin adjusting parameters
$sourceFormat = "hdf5"; //TODO: sourceFormat always hdf5?
$outputFile = "tassel/output/O" . $uId;
$taxaFile = "tassel/input/Q" . $uId . ".taxainfo";
$tassel = "tassel/tassel-gt-server-david-SPLIT.jar";
//Set sourcefile
$sourceFile = "tassel/".$data_set.".h5";

//Set output extension. Currently .hmp.txt or .vcf or .json
if ($output_format == "hapmap") {
    $output_extension = ".hmp.txt";
} elseif ($output_format == "vcf") {
    $output_extension = ".vcf";
} else {
    $output_extension = ".json";
}

/* Begin creating myfile, writing one taxa/line into $myfile.
  This file gets saved in input-folder. (See $taxafile). */
$myfile = fopen($taxaFile, "w") or die("Unable to open file!");
foreach ($taxaArray as &$txt) {
    fwrite($myfile, $txt);
    fwrite($myfile, "\n");
}
fclose($myfile);

$query_tassel = "java -jar " . $tassel . " sliceSplitFile -bf -sf " . $sourceFile . " -st " . $sourceFormat . " -df " . $outputFile . " -dt " . $output_format . " -tf " . $taxaFile . " -ch " . $chromosome . " -split " . $results_max;

if ($positions == "range") {
    $query_tassel = $query_tassel . "  -start " . $startPosition . "  -end " . $endPosition;
}
$query_log_to_DB = "INSERT INTO genotype_queries (token, tstamp, extension, files, format, resultspp, query, version) "
            . "VALUES ('" . $uId . "', " . $_SERVER["REQUEST_TIME"] . ", '" . $output_extension
            . "', 1, '" . $output_format . "', null, '" . $query_tassel . "', '".$assembly."');";
$db_handle->runQuery($query_log_to_DB);

//fwrite($myfile, $query);  //for debugging
//fclose($myfile);          //for debugging

// Log time
$time_pre_query = microtime(true);
$file_count  = exec($query_tassel);
$time_post_query = microtime(true);
error_log("\n{TIME} " . $uId. " Tassel time: " . ($time_post_query - $time_pre_query) . " (s)\n");
addDataPoint($data_set,count($taxaArray),$startPosition,$endPosition,($time_post_query - $time_pre_query));

//Split large file:
//$first_file = "http://david1.usda.iastate.edu/Diversity/tassel/output/1_O".$uId.$outputExtension;
if ($output_extension == ".json") {
  /*  $query_split_files = "cd ./tassel/output/ && ./SplitFile.py O" . $uId . $output_extension . " " . $results_max. " T"; //filepath, 200rows/page, delete original after split = true
    $time_pre_file_split = microtime(true);
    $file_count = exec($query_split_files);  //returns # of generated files
    $time_post_file_split = microtime(true);
    error_log("\n{TIME} " . $uId. " File Split time: " . ($time_post_file_split - $time_pre_file_split) . " (s)\n");*/
    $update_log_to_DB = "UPDATE genotype_queries SET resultspp=".$results_max.", files=".$file_count." WHERE token='".$uId."';";
    $db_handle->runQuery($update_log_to_DB);
    echo $uId . ',' . $file_count . ',' . $results_max; //uid,#files,#rows/page
} else {
    //not json-file
    echo "http://david1.usda.iastate.edu/Diversity/" . $outputFile . $output_extension;
}
?>