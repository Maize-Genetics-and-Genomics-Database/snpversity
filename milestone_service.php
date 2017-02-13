<?php

//echo $_POST['data'];
$jsonObj=$_POST['data'];
$data=json_decode($jsonObj,true);

$chromosome = $data["chromosome"];
$position = $data["position"]; 
$taxa = $data["taxa"];
$assembly = $data["dataSet"];
$version = $data["version"];

$uId=microtime(true)*10000;
if($assembly=="amesTest")
	{$sourceFile = "tassel/amestest.hdf5";}
elseif($assembly=="AllZeaGBSv27public20140528")
	{$sourceFile = "tassel/AllZeaGBSv27public20140528.h5";}
elseif($assembly=="ZeaGBSv27publicImputed20150114")
	{$sourceFile = "tassel/ZeaGBSv27publicImputed20150114.h5";}

if($version=="stable")
	{$tassel="tassel/sNuke.jar";}
elseif($version=="stable-latest")
	{$tassel="tassel/tassel-gt-server-new.jar";}
elseif($version=="alpha")
	{$tassel="tassel/tassel-gt-server-new.jar";}

$output_format="hapmap";
$sourceFormat = "hdf5";
$outputFile = "tassel/output/O".$uId;
$taxaFile = "tassel/input/Q".$uId.".taxainfo";

$myfile = fopen($taxaFile, "w") or die("Unable to open file!");
fwrite($myfile, $taxa);
fwrite($myfile,"\n");
fclose($myfile);

$query_log_to_DB = "java -Xms8G -Xmx8G -jar " . $tassel . " sliceNuke -bf -sf " . $sourceFile . " -st " . $sourceFormat . " -df " . $outputFile . " -dt " . $output_format . " -tf ". $taxaFile . " -ch " . $chromosome . "  -start " .  $position . "  -end " . $position;
exec($query_log_to_DB);
echo "http://david1.usda.iastate.edu/Diversity/".$outputFile.".txt";

?>
