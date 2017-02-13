<?php
require_once("db_controller.php");
require_once("TaxaExtractor.php");
/**
 * Created by PhpStorm.
 * User: root
 * Date: 11/21/16
 * Time: 11:12 AM
 */


if (isset($_POST["dataSet"]) && isset($_POST["startPosition"]) && isset($_POST["endPosition"]) && (isset($_POST["taxa"]) || isset($_FILES["stockFile"]))) {
    // Calculate range
    $range = $_POST["endPosition"] - $_POST["startPosition"];
    // Pass file as arg
    if (isset($_FILES["stockFile"])){
        getTime($_POST["dataSet"],$range, $_POST["taxa"],$_FILES["stockFile"]);
    }
    else{
        getTime($_POST["dataSet"],$range, $_POST["taxa"]);
    }
}

function getTime($dataset, $range, $taxaArray, $file=''){
    $taxaExtrct = new TaxaExtractor($dataset);
    $taxaExtrct->setTaxaArray($taxaArray);
    // Include files
    if ($file != ''){
        $taxaExtrct->setTaxaDataWithFiles($taxaArray, $file['name'], $file['tmp_name']);
    }
    // Count Stocks
    $taxaArray = $taxaExtrct->extract();
    $count = count($taxaArray);
    // Define, execute, and return query as JSON-format
    $time_query = sprintf("cd time_estimate/ && /usr/local/bin/python2.7 fetch_time.py %s %d %s", $dataset, $count, $range);
    $results = array(
        "time" => intval(exec($time_query)),
        "stocks" => $count,
        "range" => $range
        );
    echo json_encode($results);
}