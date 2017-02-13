<?php
/**
 * Created by PhpStorm.
 * User: root
 * Date: 11/21/16
 * Time: 1:44 PM
 */

$time_estimate_dir = "time_estimate/";

function addDataPoint($dataset,$stocks_count,$start,$end,$time){
    global $time_estimate_dir;
    $range = $end - $start;
    $handle = fopen($time_estimate_dir . $dataset . ".csv", "a");
    $line = array($stocks_count,$range,$time);
    fputcsv($handle,$line);
    fclose($handle);
}