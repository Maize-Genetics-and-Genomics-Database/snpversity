<!doctype html public "-//w3c//dtd html 3.2//en">
<html lang="en">
<head>
    <title>Query Genotype Data</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="css/chosen.css?version=0.2">
    <link rel="stylesheet" href="css/nucleotide_colors.css?version=1.1">
    <link rel="stylesheet" href="css/taxa_colors.css?version=0.1">
    <link rel="stylesheet" href="css/home.css?version=0.0.2">
    <link rel="stylesheet" href="css/modal.css?version=1.9.1">
    <link rel="stylesheet" href="//code.jquery.com/ui/1.11.4/themes/smoothness/jquery-ui.css">
    <link rel="stylesheet" href="css/font-awesome.min.css">
</head>
<body>
<div class="notify-box"></div>
<?php
error_reporting(0);
require_once("db_controller.php");
$db_handle = new DBController();
$id = md5(uniqid(rand(), true));
$uId = "send.php?query=" . $id; //if combine POST and GET change form-action to: action="' . $uId . '"  // Add data to URL (GET)
echo '<form id="query" action="' . $uId . '" onsubmit="return validateForm(true);" method="post" enctype="multipart/form-data">';  // POST actual Form (want unlimited stocks)
echo '<input name="query" value="' . $id . '" style="display: none;" type="hidden"/>';
?>
<table class="table table-bordered" align="left">
    <thead>
    <tr>
        <th colspan="2">
            <div><strong>Query Genotype Data</strong></div>
            <a href="#helpModal">Help</a></th>
    </tr>
    </thead>
    <tbody>
    <tr>
        <td> Quick Select / Example:
        </td>
        <td>
            <div class="btn-group" role="group">
                <button type="button" id="select_small"
                        onclick="selectSmallQuery(document.getElementById('select_small').value)"
                        value="Imputation Test">Small
                </button>
                <button type="button" id="select_medium"
                        onclick="selectMediumQuery(document.getElementById('select_medium').value)" value="Ames282">
                    Medium
                </button>
                <input type="button" class="btn btn-default" value="2010 Ames Lines" id="select_2010AmesLines"
                       onclick="selectMediumQuery(document.getElementById('select_2010AmesLines').value)"/>
                <input type="button" class="btn btn-default" value="NAM" id="select_NAM"
                       onclick="selectMediumQuery(document.getElementById('select_NAM').value)"/>
            </div>
        </td>
    </tr>
    <tr>
        <td> Assembly and Data Set:
        </td>
        <td>
            <select id="assembly" name="assembly" onchange="onAssemblyChange();">
                <option value="v2">B73 RefGen_v2</option>
                <option value="v3">B73 RefGen_v3</option>
            </select>

            <select id="dataSet" onChange="" name="dataSet" required>
                <option value="ZeaGBSv27publicImputed20150114">ZeaGBSv27publicImputed20150114</option>
                <option value="AllZeaGBSv27public20140528">AllZeaGBSv27public20140528</option>
            </select>
        </td>
    </tr>
    <tr>
        <td>
            Select Project
            <br>
            and Stock:
        </td>
        <td>
            <span id="projects"><select data-placeholder="Select Project" id="project" name="project"
                                        class="chosen-select" multiple onChange="getGBSTaxa();">
                    <option value=""></option>
                    <?php
                    $query_log_to_DB = "SELECT DISTINCT project FROM allzeagbsv27;";
                    $results = $db_handle->runQuery($query_log_to_DB);

                    foreach ($results as $row) {
                        echo '<option value="' . $row["project"] . '">' . $row["project"] . '</option>';
                    }
                    ?>
                    <option value="All">All</option>
                </select></span>
            <select data-placeholder="Select Stock" name="taxa[]" id="taxa" class="chosen-select" multiple>
                <option value=""></option>
            </select>
            <a href="#taxaSelectModal" title="Download .stockinfo files to build custom files">View Sources</a>
        </td>
    </tr>
    <tr>
        <td>Custom Upload:</td>
        <td>
            <div>
                <input type="file" value="From File" id="stockFile" name="stockFile"/>
            </div>
            <div>
                <a href="example.csv" title="Comma separated value example file (with comments)" download>Download .CSV
                    File</a>
            </div>
        </td>
    </tr>
    <tr>
        <td>Select by Gene Model or Position: </td>
        <td>
            <input type="radio" name="select-region-type" value="site" checked onChange="onRegionTypeChange(this);"> Position
            <br>
            <input type="radio" name="select-region-type" value="model" onChange="onRegionTypeChange(this);"> Gene Model
        </td>
    </tr>
    <tr id="tr-model-select">
        <td> Select Gene Model and Offset: </td>
        <td>
            <div class="ui-widget">
                <input type="text" placeholder="Enter Gene Model" id="model-select" tabindex="7" oninput="populateGeneModelTags(this.value);">
                <div class="verify-box">
                </div>
            </div>
            <div>
                <label for="offset">+/- (bp): </label>
                <input type="text" value="1000" id="offset" oninput="getRangeOfGeneModel($('#model-select').val(),this.value)">
            </div>
        </td>
    </tr>
    <tr>
        <td> Chromosome:
        </td>
        <td>
            <select id="chromosome" name="chromosome" required onchange="onChromosomeChange();">
                <option value="0">0</option>
                <option value="1">1</option>
                <option value="2">2</option>
                <option value="3">3</option>
                <option value="4">4</option>
                <option value="5">5</option>
                <option value="6">6</option>
                <option value="7">7</option>
                <option value="8">8</option>
                <option value="9">9</option>
                <option value="10">10</option>
            </select>
        </td>
    </tr>
    <tr>
        <td>Positions (bp):
        </td>
        <td>
            <select id="positions" name="positions" onchange="onPositionsChange()">
                <option value="range">range</option>
                <option value="all">all</option>
            </select>
            <label for="startPosition" id="startPositionLabel"> Start: </label>
            <input id="startPosition" name="startPosition" type="number" min="0" step="1" required/>
            <label for="endPosition" id="endPositionLabel"> End: </label>
            <input id="endPosition" name="endPosition" type="number" min="5518" step="1" required/>
        </td>
    </tr>
    <tr>
        <td> Output Format:
        </td>
        <td>
            <select id="outputFormat" name="outputFormat" onchange="onOutputFormatChange();" required>
                <option value="json">browser</option>
                <option value="hapmap">hapmap</option>
                <option value="vcf">vcf</option>
            </select>
            <label for="resultsMax" id="resultsMaxLabel">SNP's per page: </label>
            <input id="resultsMax" name="resultsMax" type="number" min="10" step="10" required/>
        </td>
    </tr>
    <tr>
        <td colspan="2"><input type="submit" value="Submit Query"/>
            <button type="button" onclick="estimate();" value="Estimate">Estimate Query Time</button>
            <div class="alert alert-warning" role="alert" name='errorMsg' id='errorMsg'></div>
        </td>
    </tr>
    </tbody>
</table>
</form>
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
<div id="aboutModal" class="modalDialog">
</div>
<div id="loading">
    <h3>This may take a few minutes...</h3>
    <img src="img/ajax-loader.gif" alt=""/>
</div>
<script src="https://code.jquery.com/jquery-1.12.4.js"></script>
<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
<script type="text/javascript" src="js/chosen.jquery.js"></script>
<script type="text/javascript" src="js/blockUI.js"></script>
<script type="text/javascript" src="js/jquery.noty.packaged.min.js"></script>
<script type="text/javascript" src="js/modalLoader.js?version=0.1.2"></script>
<script type="text/javascript" src="js/home.js?version=0.3.5"></script>
</body>
</html>

