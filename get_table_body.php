<?php

error_reporting(0);
require_once('db_controller.php');
$db_handle = new DBController();
$page_absolute = json_decode($_POST["file"]);
$assembly = json_decode($_POST["assembly"]);
$page_relative = substr($page_absolute, 41);
$handle = fopen($page_absolute, "r") or die("ERROR");
$sites_total = countObjectsFast($page_relative);
$table_name = "b73" . $assembly . "ranges";
$gene_models_all = array();
$gene_models_end = array();
if ($handle) {
    generateTableBody($handle);
}
fclose($handle);

/* Functions */

//Used to determine how many rows are needed in the table [+debugging]
function countObjectsFast($file) {
    return exec('fgrep -o "]}" ' . $file . ' | wc -l');
}

/* Prints entire table body */
function generateTableBody($file_pointer) {
    $site_no = 0; //counts current row
    global $sites_total, $gene_models_all;
    while ($site_no < $sites_total) {
        $site_no++;
        $row_result = getRowFromJson($file_pointer);
        if ($site_no == 1) {
            $gene_models_all = getAllDBVals($row_result);
        }
        printTableRow($row_result);
        unset($row_result);
    }
}

/* Prints (color-coded) row of main results table */
function printTableRow($row_result) {
    global $assembly, $gene_models_all, $gene_models_end;
    $gene_url = generateGbrowseURL($row_result, $assembly);
    echo $gene_url;

    while ($row_result["chrom_pos"] > $gene_models_all[0]["ends"] && isset($gene_models_all)) {
        if (empty($gene_models_all[0])){
            $gene_models_end = getLastValOfChrom($row_result);
            $gene_models_all = $gene_models_end;
        }
        elseif ($gene_models_all == $gene_models_end){
            break;
        }
        else{
        array_shift(&$gene_models_all);  // Remove elements that are not in range
        }
    }
    echo '<tr>
                <td style="background-color:white;">' . $row_result["rs#"] . '</td>
                <td style="background-color:white;">' . $row_result["allele"] . '</td>
                <td style="background-color:white;">' . $row_result["chrom_name"] . '</td>
                <td style="background-color:white;">' . $row_result["chrom_pos"] . '</td>';
        // Deal with IGR's
        if ($gene_models_all[0]["type"] == "IGR") {
        echo ' <td style="background-color:white;"><a href="' . $gene_url . '" title="View in Gbrowse" target="_blank">View</a></td>'
        . '<td style="background-color:white;">IGR</td>';
    }
    // Deal with other Gene Models
    elseif (!empty($gene_models_all[0]["model"])) {  // => not IGR
        $i = 0;
        $html_tds = "";
        try {
            $html_tds = $html_tds . ' <td style="background-color:white;">';
            while ($gene_models_all[$i]["pos"] <= $row_result["chrom_pos"]) {
                $html_tds = $html_tds . '<a href="' . $gene_url . '" title="View in Gbrowse" target="_blank">' . substr($gene_models_all[$i]["model"], 0, -4) . '</a><br>';
                $i++;
                if (isset($gene_models_end)){break;}  // Reached final result already (when close to end of chromosome).
                if ($i > 10){
                    throw new RangeException("Counted " . strval($i) . " different types for position " . strval($row_result["chrom_pos"]));
                }
            }
            $html_tds = $html_tds . '</td><td style="background-color:white;">';
            for ($j = 0; $j < $i; $j++) {
                $html_tds = $html_tds . $gene_models_all[$j]["type"] . '<br>';
            }
            // Check if empty or not
            if ($html_tds == ' <td style="background-color:white;"></td><td style="background-color:white;">'){
                throw new RangeException("No matches found");
            }
            echo $html_tds . '</td>';
        }
        catch(RangeException $e) {
            echo ' <td style="background-color:white;"><a href="' . $gene_url . '" title="View in Gbrowse" target="_blank">View</a></td>'
                . '<td style="background-color:white;">?</td>';
        }
    }
    // Else, not found in ranges DB
    else{
        echo ' <td style="background-color:white;"><a href="' . $gene_url . '" title="View in Gbrowse" target="_blank">View</a></td>'
        . '<td style="background-color:white;">?</td>'; 
    }
    foreach ($row_result["results"] as $allele) {
        $pos = strpos($row_result["allele"], $allele);
        printCell($pos, $allele);
    }
    echo '</tr>';
}

function generateGbrowseURL($row_result, $gbrowser_version) {
    $common_url = "http://www.maizegdb.org/gbrowse/maize_" . $gbrowser_version;
    if ($row_result["chrom_name"] == '0') {
        $chr_full = "UNMAPPED";
    } else {
        $chr_full = "Chr" . $row_result["chrom_name"];
    }

    $gene_url = $common_url . "/?start=" . strval((intval($row_result["chrom_pos"]) - 500))
            . ";stop=" . strval((intval($row_result["chrom_pos"]) + 500))
            . ";width=1000;ref=" . $chr_full . ";h_region=" .
            $row_result["chrom_pos"] . ".." . $row_result["chrom_pos"] . "@red;";
    return $gene_url;
}

/* Contains decision logic for highlighting a cell. */

function printCell($pos, $allele) {
    //Color-coding for different/same alleles as $result["allele"]
    if ($pos === 0 && $allele !== 'N') {
        echo '<td class="j">' . $allele . '</td>'; //major
    } else if ($pos === 2) {
        echo '<td class="n">' . $allele . '</td>';  //minor
    } else { // there could be '-' or '+' on position 4, which we want to treat in switch
        switch ($allele) {
            case 'A':
                echo '<td class="A">' . $allele . '</td>';
                break;
            case 'C':
                echo '<td class="C">' . $allele . '</td>';
                break;
            case 'G':
                echo '<td class="G">' . $allele . '</td>';
                break;
            case 'T':
                echo '<td class="T">' . $allele . '</td>';
                break;
            case 'R':
                echo '<td class="R">' . $allele . '</td>';
                break;
            case 'Y':
                echo '<td class="Y">' . $allele . '</td>';
                break;
            case 'S':
                echo '<td class="S">' . $allele . '</td>';
                break;
            case 'W':
                echo '<td class="W">' . $allele . '</td>';
                break;
            case 'K':
                echo '<td class="K">' . $allele . '</td>';
                break;
            case 'M':
                echo '<td class="M">' . $allele . '</td>';
                break;
            case 'B':
                echo '<td class="B">' . $allele . '</td>';
                break;
            case 'D':
                echo '<td class="D">' . $allele . '</td>';
                break;
            case 'H':
                echo '<td class="H">' . $allele . '</td>';
                break;
            case 'V':
                echo '<td class="V">' . $allele . '</td>';
                break;
            case '+':
                echo '<td class="ins">' . $allele . '</td>';
                break;
            case '-':
                echo '<td class="del">' . $allele . '</td>';
                break;
            case '0':
                echo '<td class="zero">' . $allele . '</td>';
                break;
            case 'N':
                echo '<td class="N">' . $allele . '</td>';
                break;
            default:
                echo '<td class="noMatch">' . $allele . '</td>';
        }
    }
}

//returns value $col_name from DB for a given $row_result
function getDBVal($row_result, $col_name, $table_name) {
    global $db_handle;
    $query = "SELECT " . $col_name . " FROM " . $table_name
            . " WHERE chr=" . $row_result["chrom_name"] . " AND "
            . "pos=" . $row_result["chrom_pos"] . ";";
    $gene_model = $db_handle->runQuery($query);
    return $gene_model;
}

function getLastValOfChrom($row_result){
    global $db_handle, $gene_models_all, $table_name;
    $query_max_pos = "SELECT max(pos) FROM " . $table_name
        . " WHERE chr=" . $row_result["chrom_name"];
    $max_pos = $db_handle->runQuery($query_max_pos);

    $query_max_pos_row = "SELECT * FROM " . $table_name. " WHERE chr=" . $row_result["chrom_name"] .
        " AND " . "pos=" . $max_pos[0]["max"] . " ORDER BY chr,pos LIMIT 1";
    $final_row = $db_handle->runQuery($query_max_pos_row);
    return $final_row;
}

//Similar to getDBVals, but connects to DB only once.
function getAllDBVals($row_result) {
    global $db_handle, $sites_total, $table_name;
    $max_vals = ceil($sites_total * 1.5);
    $query = "SELECT * FROM " . $table_name
            . " WHERE chr=" . $row_result["chrom_name"] . " AND "
            . "ends>=" . $row_result["chrom_pos"] . " ORDER BY chr,pos LIMIT " . $max_vals . ";";
    $gene_models = $db_handle->runQuery($query);
    return $gene_models;
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

/* DEPRECATED: Used to determine how many rows are needed in the table (slower, but works on arbitrary large file)
  function countObjects($file, $echo_on) {
  $count = 0;
  while (($ch = fgetc($file)) !== false) {
  if ($ch === '}') { //end of JSON-Object (1 row in table)
  $count++;
  }
  }
  return $count;
  } */
?>