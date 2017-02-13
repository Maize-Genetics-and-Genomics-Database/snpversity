<?php
require_once("db_controller.php");
$projects = json_decode($_POST["projects"]);
$db_handle = new DBController();

foreach ($projects as $key => $proj) {
    if ($proj == "RandD"){
        $proj = "R&D";
    }
    if ($proj == "All"){
        $query_log_to_DB = "SELECT DISTINCT ON (dna_sample,project) dna_sample, lib_prep_id, inbred FROM allzeagbsv27";
    }
    else{
    $query_log_to_DB = "SELECT DISTINCT ON (dna_sample) dna_sample, lib_prep_id, inbred FROM allzeagbsv27 WHERE project = '" . $proj . "' ORDER BY dna_sample, lib_prep_id";
    }
    $results = $db_handle->runQuery($query_log_to_DB);
    echo '<optgroup label="' . $proj . '">';
    if ($proj != "All"){
    echo '<option value="'.$proj.'">All '.$proj.'</option>';
    }
    foreach ($results as $row) {
        if ($row["inbred"] == '' || ($row["inbred"] == $row["dna_sample"] )){
        echo '<option value="' . $row["dna_sample"] . ":" . $row["lib_prep_id"] . '">' . $row["dna_sample"] . '</option>';
        }
        else{
            echo '<option value="' . $row["dna_sample"] . ":" . $row["lib_prep_id"] . '">' . $row["inbred"] . " = ". $row["dna_sample"] . '</option>';
        }
    }
    echo '</optgroup>';
}
?>