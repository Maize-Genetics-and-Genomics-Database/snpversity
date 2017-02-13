<?php
require_once("db_controller.php");
$db_handle = new DBController();
if (isset($_POST["action"])) {
    if ($_POST["action"] == 'getGeneModels') {
        getGeneModels($db_handle,$_POST["assembly"],$_POST["input"]);
    }
    else if ($_POST["action"] == 'getRange' && isset($_POST["model"])) {
        getRangeOfModel($db_handle,$_POST["assembly"],$_POST["model"]);
    }
    else{

    }
}

function printGeneModels($db_handle, $assembly) {
    $results = $db_handle->runQuery("SELECT DISTINCT model FROM b73".$assembly."ranges LIMIT 500;");
    foreach ($results as $row) {
        echo '<option value="' . $row["model"] . '">' . $row["model"] . '</option>';
    }
}
function getGeneModels($db_handle, $assembly, $input){
    $results = $db_handle->runQuery("SELECT DISTINCT model FROM b73".$assembly."ranges WHERE model ~* '.*".$input.".*' LIMIT 50;");
    echo json_encode($results);
}
function getRangeOfModel($db_handle, $assembly,$model){
    $results = $db_handle->runQuery("SELECT chr, min(pos), max(ends) FROM b73".$assembly."ranges WHERE model='".$model."' GROUP BY chr;");
    echo json_encode($results);
}