<?php

require_once("db_controller.php");
$db_handle = new DBController();
if (isset($_POST["action"])) {
    if ($_POST["action"] == 'getHapMapLines') {
        printHapMapLines($db_handle);
    }
}

function printHapMapLines($db_handle) {
    $results = $db_handle->runQuery("SELECT * FROM hapmapv3;");
    foreach ($results as $row) {
        echo '<option value="' . $row["taxon"] . '">' . $row["taxon"] . '</option>';
    }
}
