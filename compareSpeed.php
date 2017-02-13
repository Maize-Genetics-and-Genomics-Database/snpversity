<?php
//Begin adjusting parameters
$sourceFormat = "hdf5"; //TODO: sourceFormat always hdf5?
$taxaFile = "/root/Documents/test/2500.taxainfo";
$chromosome = 1;
$startPosition = "2918";
$endPosition = "101231039";
$results_max = 50;
//Set sourcefile
$sourceFile = "tassel/ZeaGBSv27publicImputed20150114.h5";
$output_extension = ".json";
$output_format = "json";
$validate_files = false;
echo (string)($endPosition - $startPosition);
// Define different files
$output_GSON = "/root/Documents/test/GSON_Output" . $output_extension;
$tassel_GSON = "tassel/tassel-gt-server-david-SPLIT.jar";
$output_BOON = "/root/Documents/test/BOON_Output" . $output_extension;
$tassel_BOON = "tassel/tassel-gt-server-david-BOON.jar";

// Print out Runtimes:
# run_time($output_BOON,$tassel_BOON,"BOON");
 run_time($output_GSON,$tassel_GSON,"GSON");
 if ($validate_files){
 // Open files, check for valid JSON
$handle_GSON = fopen("/root/Documents/test/1_GSON_Output". $output_extension, "r") or die("ERROR");
$handle_BOON = fopen("/root/Documents/test/1_BOON_Output". $output_extension, "r") or die("ERROR");
echo "\n{GSON-PARSED} \n";
for($x = 0; $x < 4; $x++){
var_dump(getRowFromJson($handle_GSON));
}
echo "\n{BOON-PARSED} \n";
for($x = 0; $x < 4; $x++){
    var_dump(getRowFromJson($handle_BOON));
}
// Close files
fclose($handle_BOON);
fclose($handle_GSON);
 }
// Functions
function run_time($output_file,$jar_file, $name){
    global $sourceFile,$sourceFormat,$output_format, $taxaFile, $chromosome, $results_max, $startPosition, $endPosition;
$query = "java -jar " . $jar_file . " sliceSplitFile -bf -sf " . $sourceFile . " -st " . $sourceFormat . " -df " . $output_file . " -dt ". $output_format . " -tf " . $taxaFile . " -ch " . $chromosome . " -split " . $results_max . "  -start " . $startPosition . "  -end " . $endPosition;
$time_pre_query = microtime(true);
echo "\n{".$name."-QUERY} " . $query;
$file_count  = exec($query);
echo "\n File Count: " . $file_count;
$time_post_query = microtime(true);
$report_string_gson = "\n{".$name."-TIME} " . " Tassel time: " . ($time_post_query - $time_pre_query) . " (s)\n";
echo $report_string_gson;
}
    //Converts next json-object (corresponding to 1 row in results table) to an associative array
function getRowFromJson($handle) {
    $json_str = "";
    while (($ch = fgetc($handle)) !== false) {
        if ($ch == '{') { // begin
            $json_str = $json_str . $ch;
            while (($ch_next = fgetc($handle)) !== '}') {
                $json_str = $json_str . $ch_next;
                //process JSON
            }
            return json_decode($json_str . '}', true, 8192);
        }
    }
}
?>
