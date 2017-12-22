<!DOCTYPE html>
<html lang="en">
<?php
error_reporting(0);
$parent = $_GET["parent"];
echo '<div id="parent-div" style="display:none;">' . $parent . '</div>';
require_once("db_controller.php");
require_once("TaxaExtractor.php");
$db_handle = new DBController();
$dir_www = "http://snpversity.maizegdb.org/Diversity/";
$dir_out = "tassel/output/";
$dir_in = "tassel/input/Q";
$uId = $_GET["query"];
$input_file_absolute = $dir_www . $dir_in . $uId . ".taxainfo";
$projects = array();
// query was run previously, skip processing
$exists_query_string = "SELECT EXISTS (SELECT 1 FROM genotype_queries WHERE token='" . $uId . "');";
$exists_query = $db_handle->runQuery($exists_query_string);

if ($exists_query[0]["exists"] == 't') {
    //reload values from previous query in DB genotype_queries
    $retrieve_query_string = "(SELECT * FROM genotype_queries WHERE token='" . $uId . "');";
    $result = $db_handle->runQuery($retrieve_query_string);
    $output_extension = $result[0]["extension"];
    $pages = $result[0]["files"];
    $output_format = $result[0]["format"];
    $results_max = $result[0]["resultspp"];  //results per page
    $assembly = $result[0]["version"];
    $taxaExtrct = new TaxaExtractor($assembly);
    //$taxaArray = array();  //need to read in taxa-array from previous file
    //echo $outputFormat. " <- format, resultsPP -> ".$resultsPerPage;
    $taxaArray = $taxaExtrct->getTaxaFromTaxaInfoFile($input_file_absolute);  //get taxa from file
    if ($output_extension == ".json") {
        $curl_result = $uId . ',' . $pages . ',' . $results_max;
    } else {
        $curl_result = $dir_www . $dir_out . 'O' . $uId . $output_extension;
    }
    goto serviceComplete;
}

// Get data from Form
$assembly = $_POST["assembly"];
$chromosome = $_POST["chromosome"];
$positions = $_POST["positions"];
$startPosition = $_POST["startPosition"];
$endPosition = $_POST["endPosition"];
$taxaArray = $_POST["taxa"];
$output_format = $_POST["outputFormat"];
$results_max = $_POST["resultsMax"];
$data_set = $_POST["dataSet"];

$taxaFileTemp = $_FILES["stockFile"]['tmp_name'];
$taxaFileName = $_FILES["stockFile"]['name'];
$taxaExtrct = new TaxaExtractor($assembly);
$taxaExtrct->setTaxaDataWithFiles($taxaArray, $taxaFileName, $taxaFileTemp);

$taxaArray = $taxaExtrct->extract();


$Obj = array(
    "assembly" => $assembly,
    "chromosome" => $chromosome,
    "positions" => $positions,
    "startPosition" => $startPosition,
    "endPosition" => $endPosition,
    "taxaArray" => $taxaArray,
    "outputFormat" => $output_format,
    "dataSet" => $data_set,
    "resultsPerPage" => $results_max,
    "query" => $uId
);

$jsonObj = json_encode($Obj);
$url = $dir_www . "service.php";

// Initialise libcurl
$curl = curl_init();
// Set Curl options:
curl_setopt_array($curl, array(
    CURLOPT_RETURNTRANSFER => 1, // return result as string value
    CURLOPT_URL => $url, //send request to $url
    CURLOPT_USERAGENT => 'Codular Sample cURL Request',
    CURLOPT_POST => 1, //enable POST-usage
    CURLOPT_POSTFIELDS => array(//send list of variables in form key => value
        "data" => $jsonObj
    )
));
$curl_result = curl_exec($curl);
// Close and store Curl-response
$status = curl_getinfo($curl, CURLINFO_HTTP_CODE);
curl_close($curl);

serviceComplete:
if ($output_format === 'json') {
    // Execute Curl-query, store results
    $service_results = explode(",", $curl_result);  // $service_results =
    $first_page_absolute = $dir_www . $dir_out . "1_O" . $uId . "." . $output_format;
    // $first_page_relative = substr($first_page_absolute, 41); // skips "http://david1.usda.iastate.edu/Diversity/ (used for counting obj)
    $pages = $service_results[1];
} else {
    echo "<a href=" . $curl_result . ">View File</a>";
    goto skipTable;
}
?>
<?php
// Skip page if parent is not set or not true
if (!isset($parent) || $parent !== "true") {
    goto skipTable;
}
?>
<head>
    <title>Diversity Viewer</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="css/send.css?version=7.9.3">
    <link rel="stylesheet" href="css/nucleotide_colors.css?version=0.6">
    <link rel="stylesheet" href="css/taxa_colors.css?version=0.2">
    <link rel="stylesheet" href="css/modal.css?version=0.4">
    <link rel="stylesheet" href="css/jquery.toolbar.css">
    <link rel="stylesheet" href="css/font-awesome.min.css">
    <link rel="stylesheet" href="//code.jquery.com/ui/1.11.4/themes/smoothness/jquery-ui.css">
    <!--<link rel="stylesheet" type="text/css"
          href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.7.1/themes/base/jquery-ui.css"/>-->
</head>
<body>
<!-- Results Table--->
<div id="wrapper">
    <div id='content'>
        <table cellspacing="0" id="tableOfResults">
            <thead>
            <tr>
                <th class="shortHead" id='siteHeader'>
                    <div class="applyHeaderFormat">Site</div>
                </th>
                <th class="shortHead" id='alleleHeader'>
                    <div class="applyHeaderFormat">Allele</div>
                </th>
                <th class="shortHead" id='chrHeader'>
                    <div class="applyHeaderFormat">Chr</div>
                </th>
                <th class="shortHead" id='positionHeader'>
                    <div class="applyHeaderFormat">Position</div>
                </th>
                <th class="shortHead" id='geneModelHeader'>
                    <div class="applyHeaderFormat">Gene (T01)</div>
                </th>
                <th class="shortHead" id='geneTypeHeader'>
                    <div class="applyHeaderFormat">Type</div>
                </th>
                <?php
                if ($assembly == "v3") {  //Should only enter this loop for HapMap-lines!
                    if (count($taxaArray) > 20) {
                        foreach ($taxaArray as $taxon) {
                            echo '<th class="rotate"><div><span title="' . $taxon . '">' . $taxon . '</span></div></th>';
                        }
                    } else {
                        foreach ($taxaArray as $taxon) {
                            echo '<th><div><span>' . $taxon . '</span></div></th>';
                        }
                    }
                } else {  //display only first part of taxon
                    if (count($taxaArray) > 20) {
                        foreach ($taxaArray as $taxon) {
                            $cssProjCode = getCSSProjectEncoding($taxon);  //used for coloring in table header
                            $taxon = taxonToInbred($taxon);
                            echo '<th class="rotate ' . $cssProjCode . '"><div><span title="' . $taxon . '">' . $taxon . '</span></div></th>';
                        }
                    } else {
                        foreach ($taxaArray as $taxon) {
                            $cssProjCode = getCSSProjectEncoding($taxon);  //used for coloring in table header
                            echo '<th class="' . $cssProjCode . '"><div><span>' . taxonToInbred($taxon) . '</span></div></th>';
                        }
                    }
                }
                ?>
            </tr>
            </thead>
            <tbody>
            </tbody>
        </table>
    </div>
    <div id="footer">
        <div id='leftFooter'>
            <strong>Zoom:</strong>
            <input type="text" id="amount" readonly
                   style="background-color: #61A961; color: white; border: none; display: inline">
            <div id="slider">
            </div>
        </div>
        <div id='rightFooter'>
            <strong>Color Labels:</strong>
            <div>
                <?php if ($assembly == "v2" || $assembly == "v4"){echo '<button style="display:inline;" type="button" id="showCodeT" onclick="showTaxaCodes();">Show Stock Box</button>';}?>

                <button style="display:inline;" type="button" id="showCodeN" onclick="showNucleotideCodes();">Show Nucleotide Box</button>
            </div>
        </div>
        <div id='centerFooter' align="center">
            <div id="descriptionFooter">
                <strong>Assembly: <?php echo "B73RefGen" . ucfirst($assembly) ?></strong>
            </div>
            <a href="#helpModal" class="tool-item btn-toolbar-dark"><i class="fa fa-question fa-lg"></i></a>
            <a href="<?php echo $input_file_absolute; ?>" target="_blank" class="tool-item btn-toolbar-dark"><i
                    class="fa fa-file-text-o fa-lg"></i></a>
            <a href="" class="tool-item btn-toolbar-dark convertCSV"><i class="fa fa-download fa-lg"></i></a>
            <div style="width: 50%;">
                <?php
                //Select page (if query is big)
                if ($pages > 1) {
                    echo '<button type="button" style="display:inline; vertical-align:bottom;" onclick="prevPg();">Prev</button>';
                    echo '<select id="pages" name="pages" style="display:inline; vertical-align:bottom;" onchange="rePopulate();" required>';
                    foreach (range(1, $pages) as $no) {
                        $file = $dir_www . $dir_out . $no . '_O' . $uId . '.' . $output_format;
                        $next_file = $dir_www . $dir_out . strval(((int)$no) + 1) . '_O' . $uId . '.' . $output_format;

                        $start_site = getFirstSiteFromJson($file);
                        $end_site = getFirstSiteFromJson($next_file);
                        if (is_numeric($end_site)) {
                            $end_site = (string)((int)$end_site - 1);
                        }
                        $range_str = (string)$start_site . ' - ' . $end_site;

                        echo '<option value="' . $file . '">' . $range_str . '</option>';
                    }
                    echo '</select><button type="button" style="display:inline;vertical-align:bottom;" onclick="nextPg();">Next</button>';
                }
                ?>
            </div>
        </div>
    </div>
</div>
<!-- divs hidden from main page -->
<div id="taxaTable">
</div>

<div id="nucleotideTable">
</div>
<div id="nucleotideModal" class="modalDialog">
</div>
<div id="geneModelsModal" class="modalDialog">
</div>
<div id="taxaCodesModal" class="modalDialog">
</div>
<div id="helpModal" class="modalDialog">
</div>
<div id="fileModal" class="modalDialog">
</div>
<div id="taxaSelectModal" class="modalDialog">
</div>
<div id="helpModal" class="modalDialog">
</div>
<div id="allzeagbsModal" class="modalDialog">
</div>
<div id="hapmapV3Modal" class="modalDialog">
</div>
<div id="geneModelsModal" class="modalDialog">
</div>
<div id="datasetsModal" class="modalDialog">
</div>
<div id="first-page" class="hidden"><?php echo $first_page_absolute; ?></div>
<div id="version" class="hidden"><?php echo $assembly; ?></div>
<div id="query-id" class="hidden"><?php echo $uId; ?></div>
<div id="loading">
    <h3>Constructing table...</h3>
    <img src="img/ajax-loader.gif" alt=""/>
</div>
<?php
skipTable:
?>
<script type="text/javascript" src="js/jquery.js"></script>
<script type="text/javascript" src="js/jquery.toolbar.js"></script>
<script type="text/javascript" src="js/blockUI.js"></script>
<script type="text/javascript" src="//code.jquery.com/ui/1.11.4/jquery-ui.js"></script>
<script type="text/javascript" src="js/modalLoader.js"></script>
<script type="text/javascript" src="js/stickytableheaders.min.js"></script>
<script type="text/javascript" src="js/send.js?version=1.2.9"></script>
</body>
</html>
<?php
/* Functions */

//Converts next json-object (corresponding to 1 row/SNIP in results table) to an associative array. Used for populating dropdown-menu of positions
function getFirstSiteFromJson($file)
{
    $file_relative = substr($file, 41);  // skips "http://david1.usda.iastate.edu/Diversity/ (dont want to use full url)

    if (file_exists($file_relative)) {
        $handle = fopen($file_relative, "r") or die("Could not open page based on file " . $file);
        $json_str = "";
        while (($ch = fgetc($handle)) !== false) {
            if ($ch == '{') { // begin
                $json_str = $json_str . $ch;
                while (($ch_next = fgetc($handle)) !== '}') {
                    $json_str = $json_str . $ch_next;
                    //process JSON
                }
                $json_arr = json_decode($json_str . '}', true, 8192);
                return strval($json_arr["chrom_pos"]);
            }
        }
    } else {
        return "End";
    }
}

function getCSSProjectEncoding($taxon)
{
    global $db_handle;
    $taxonClass = "";
    $taxon_id = substr($taxon, -9);
    $projOfTaxon = $db_handle->runQuery(
        "SELECT project "
        . "FROM allzeagbsv27 "
        . "WHERE lib_prep_id='" . $taxon_id . "';");
    switch ($projOfTaxon[0]['project']) {
        case "NAM":
            $taxonClass = "NAM";
            break;
        case "IBM":
            $taxonClass = "IBM";
            break;
        case "ApeKI 384-plex":
            $taxonClass = "ApeKI";
            break;
        case "Imputation Test":
            $taxonClass = "Imputation";
            break;
        case "2010 Ames Lines":
            $taxonClass = "Ames2010";
            break;
        case "R&D":
            $taxonClass = "RandD";
            break;
        case "Maize-BREAD":
            $taxonClass = "BREAD";
            break;
        case "AMES Inbreds":
            $taxonClass = "Inbreds";
            break;
        case "Ames282":
            $taxonClass = "Ames282";
            break;
        case "Old Maize Diversity":
            $taxonClass = "oldDiversity";
            break;
        default :
            break;
    }
    return $taxonClass;
}
function taxonToInbred($taxon){
    global $db_handle;
    $taxon_lib_prep_id = substr($taxon, -9);
    $inbreds = $db_handle->runQuery(
        "SELECT inbred "
        . "FROM allzeagbsv27 "
        . "WHERE lib_prep_id='" . $taxon_lib_prep_id . "';");
    if ($inbreds[0]["inbred"] != ""){
        return $inbreds[0]["inbred"];
    }
    return substr($taxon, 0, -10);
}

?>
