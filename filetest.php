<?php
/**
 * Created by PhpStorm.
 * User: root
 * Date: 9/19/16
 * Time: 12:20 PM
 */
require_once("TaxaExtractor.php");
$taxaExtrct = new TaxaExtractor("v2");
$taxaArray = array(
    1 => "12E:250032344",
    2 => "user-select2"
);
//var_dump($taxaExtrct->checkIfInbred("807",false));
$taxaFile = "/root/Documents/2010AmesTest.stockinfo";
$taxaFileName = "taxafile.stockinfo";
$taxaExtrct->setTaxaDataWithFiles($taxaArray,$taxaFileName,$taxaFile);
$taxaArray = $taxaExtrct->extract();

var_dump($taxaArray);